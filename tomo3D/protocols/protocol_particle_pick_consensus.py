# **************************************************************************
# *
# * Authors:    Carlos Oscar Sorzano (coss@cnb.csic.es)
# *             Tomas Majtner (tmajtner@cnb.csic.es)  -- streaming version
# *             David Herreros Calero (dherreros@cnb.csic.es) -- Tomo version
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'coss@cnb.csic.es'
# *
# **************************************************************************
"""
Consensus picking protocol
"""

import os

from math import sqrt
import numpy as np

from pyworkflow.object import Set, Pointer
import pyworkflow.protocol.params as params
from pyworkflow.protocol.constants import *
from pyworkflow.utils import getFiles, removeBaseExt, moveFile, removeExt

from pwem.convert.transformations import quaternion_from_matrix, weighted_tensor,\
                                         mean_quaternion, quaternion_matrix

from tomo.protocols import ProtTomoPicking
from tomo.objects import SetOfCoordinates3D, Coordinate3D

PICK_MODE_LARGER = 0
PICK_MODE_EQUAL = 1

class ProtTomoConsensusPicking(ProtTomoPicking):
    """
    Protocol to estimate the agreement between different particle picking
    algorithms. The protocol takes several Sets of Coordinates calculated
    by different programs and/or different parameter settings. Let's say:
    we consider N independent pickings. Then, a coordinate is considered
    to be a correct particle if M pickers have selected the same particle
    (within a radius in pixels specified in the form).

    If you want to be very strict, then set M=N; that is, a coordinate
    represents a particle if it has been selected by all particles (this
    is the default behaviour). Then you may relax this condition by setting
    M=N-1, N-2, ...

    If you want to be very flexible, set M=1, in this way it suffices that
    1 picker has selected the coordinate to be considered as a particle. Note
    that in this way, the cleaning of the dataset has to be performed by other
    means (screen particles, 2D and 3D classification, ...).
    """

    _label = 'picking consensus'
    outputName = 'consensusCoordinates'
    FN_PREFIX = 'consensusCoords_'

    def __init__(self, **args):
        ProtTomoPicking.__init__(self, **args)
        self.stepsExecutionMode = STEPS_SERIAL

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.MultiPointerParam,
                      pointerClass='SetOfCoordinates3D',
                      label="Input 3D coordinates", important=True,
                      help='Select the set of 3D coordinates to compare')
        form.addParam('consensusRadius', params.IntParam, default=10,
                      label="Radius",
                      help="All coordinates within this radius (in pixels) "
                           "are presumed to correspond to the same particle")
        form.addParam('consensus', params.IntParam, default=-1,
                      label="Consensus",
                      help="This parameter can take values from 1 to total number of inputs (being "
                           "-1 a special case).\n"
                           "*Set to -1* to indicate that it needs to be selected "
                           "by all algorithms: *AND* operation.\n"
                           "*Set to 1* to indicate that it suffices that only "
                           "1 algorithm selects the particle: *OR* operation.\n"
                           "Any other value will determine how many times need a particle"
                           " to be selected to be considered as a consensus particle.")
        form.addParam('mode', params.EnumParam, label='Consensus mode',
                      choices=['>=', '='], default=PICK_MODE_LARGER,
                      expertLevel=LEVEL_ADVANCED,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='If the number of votes to progress to the output '
                           'must be either (=) strictly speaking equals to '
                           'the consensus number or (>=) at least equals.')

        # FIXME: It's not using more than one since
        #         self.stepsExecutionMode = STEPS_SERIAL
        # form.addParallelSection(threads=4, mpi=0)

#--------------------------- INSERT steps functions ---------------------------
    def _insertAllSteps(self):
        self.checkedTomos = set()   # those tomos ready to be processed (tomoId)
        self.processedTomos = set() # those tomos already processed (tomoId)
        self.sampligRates = []
        coorSteps = self.insertNewCoorsSteps([])
        self._insertFunctionStep('createOutputStep',
                                 prerequisites=coorSteps, wait=True)

    def createOutputStep(self):
        pass

    def _getFirstJoinStepName(self):
        # This function will be used for streaming, to check which is
        # the first function that need to wait for all mics
        # to have completed, this can be overriden in subclasses
        # (e.g., in Xmipp 'sortPSDStep')
        return 'createOutputStep'

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == self._getFirstJoinStepName():
                return s
        return None

    def insertNewCoorsSteps(self, tomos):
        deps = []
        for tomo in tomos:
            stepId = self._insertFunctionStep("calculateConsensusStep",
                                              tomo.getObjId(),
                                              tomo.getFileName(),
                                              prerequisites=[])
            deps.append(stepId)
        return deps

    def _stepsCheck(self):
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        # If continue from an stopped run, don't repeat what is done
        if not self.checkedTomos:
            for fn in getFiles(self._getExtraPath()):
                fn = removeBaseExt(fn)
                if fn.startswith(self.FN_PREFIX):
                    self.checkedTomos.update([self.getTomoId(fn)])
                    self.processedTomos.update([self.getTomoId(fn)])

        streamClosed = []
        readyTomos = None
        allTomos = set()
        for coordSet in self.inputCoordinates:
            currentPickTomos, isSetClosed = getReadyTomos(coordSet.get())
            streamClosed.append(isSetClosed)
            if not readyTomos:  # first time
                readyTomos = currentPickTomos
            else:  # available mics are those ready for all pickers
                readyTomos.intersection_update(currentPickTomos)
            allTomos = allTomos.union(currentPickTomos)

        self.streamClosed = all(streamClosed)
        if self.streamClosed:
            # for non streaming do all and in the last iteration of streaming do the rest
            newTomoIds = allTomos.difference(self.checkedTomos)
        else:  # for streaming processing, only go for the ready mics in all pickers
            newTomoIds = readyTomos.difference(self.checkedTomos)

        if newTomoIds:
            self.checkedTomos.update(newTomoIds)

            inTomos = self.getMainInput().getPrecedents()
            newTomos = [inTomos[tomoId].clone() for tomoId in newTomoIds]

            fDeps = self.insertNewCoorsSteps(newTomos)
            outputStep = self._getFirstJoinStep()
            if outputStep is not None:
                outputStep.addPrerequisites(*fDeps)
            self.updateSteps()

    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return
        self.finished = self.streamClosed and self.checkedTomos == self.processedTomos
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        newFiles = getFiles(self._getTmpPath())
        if newFiles or self.finished:  # when finished to close the output set
            outSet = self._loadOutputSet(SetOfCoordinates3D, 'coordinates.sqlite')

            for fnTmp in newFiles:
                if '_vIds' in fnTmp or '_trMats' in fnTmp:
                    continue
                coords = np.loadtxt(fnTmp)
                vsIds = np.loadtxt(removeExt(fnTmp) + '_vIds.txt')
                with open(removeExt(fnTmp) + '_trMats.txt') as outFile:
                    shape_mat = outFile.readline()
                    shape_mat = eval('(' + shape_mat.split('(')[1])
                    trMats = np.loadtxt(removeExt(fnTmp) + '_trMats.txt').reshape(shape_mat)
                moveFile(fnTmp, self._getExtraPath())
                moveFile(removeExt(fnTmp) + '_vIds.txt', self._getExtraPath())
                for idv in np.unique(vsIds):
                    vesicle_coords = coords[np.where(vsIds == idv)]
                    vesicle_tr = trMats[np.where(vsIds == idv)]
                    if vesicle_coords.size == 3:  # special case with only one coordinate
                        vesicle_coords = [vesicle_coords]
                    for idx, coord in enumerate(vesicle_coords):
                        newCoord = Coordinate3D()
                        tomograms = self.getMainInput().getPrecedents()
                        newCoord.setVolume(tomograms[self.getTomoId(fnTmp)])
                        newCoord.setPosition(coord[0], coord[1], coord[2])
                        newCoord.setGroupId(idv)
                        matrix = vesicle_tr[idx]
                        if isinstance(self.inputCoordinates, list):
                            newCoord.setMatrix(matrix)
                        else:
                            newCoord.setMatrix(matrix)
                        outSet.append(newCoord)

            firstTime = not self.hasAttribute(self.outputName)
            self._updateOutputSet(self.outputName, outSet, streamMode)
            if firstTime:
                self.defineRelations(outSet)
            outSet.close()

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(STATUS_NEW)

    def defineRelations(self, outputSet):
        for inCorrds in self.inputCoordinates:
            self._defineTransformRelation(inCorrds, outputSet)

    def _loadOutputSet(self, SetClass, baseName):
        setFile = self._getPath(baseName)
        if os.path.exists(setFile):
            outputSet = SetClass(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

        #inMicsPointer = Pointer(self.getMapper().getParent(
        #                                    self.getMainInput().getMicrographs()),
        #                                    extended='outputMicrographs')
        outputSet.setBoxSize(self.getMainInput().getBoxSize())
        outputSet.setSamplingRate(self.getMainInput().getSamplingRate())
        inTomosPointer = Pointer(self.getMainInput().getPrecedents())
        outputSet.setPrecedents(inTomosPointer)

        return outputSet

    def calculateConsensusStep(self, tomoId, tomoName):

        print("Consensus calculation for Tomogram %d: '%s'"
              % (tomoId, tomoName))

        # Take the sampling rates just once
        if not self.sampligRates:
            for coordinates in self.inputCoordinates:
                tomograms = coordinates.get().getPrecedents()
                self.sampligRates.append(tomograms.getSamplingRate())

        # Get all coordinates for this tomogram
        coords = []
        vesicles = []
        transformations = []
        for idx, coordinates in enumerate(self.inputCoordinates):
            coordArray = np.asarray([x.getPosition() for x in
                                     coordinates.get().iterCoordinates(tomoId)],
                                     dtype=float)
            vIds = np.asarray([x.getGroupId() for x in
                               coordinates.get().iterCoordinates(tomoId)])
            trMat = np.asarray([x.getMatrix() for x in
                                coordinates.get().iterCoordinates(tomoId)])
            coordArray *= float(self.sampligRates[idx]) / float(self.sampligRates[0])
            coords.append(np.asarray(coordArray, dtype=int))
            vesicles.append(np.asarray(vIds, dtype=int))
            transformations.append(trMat)

        consensusWorker(coords, vesicles, transformations, self.consensus.get(),
                        self.consensusRadius.get(),
                        self._getTmpPath('%s%s.txt' % (self.FN_PREFIX, tomoId)),
                        self._getExtraPath('jaccard.txt'), self.mode.get())

        self.processedTomos.update([tomoId])

    def _validate(self):

        errors = []

        # Only for Scipion 2.0, next versions should have the default
        # PointerList validation and this can be removed
        if len(self.inputCoordinates) == 0:
                errors.append('inputCoordinates cannot be EMPTY.')
        # Consider empty pointers:
        else:
            for pointer in self.inputCoordinates:
                obj = pointer.get()
                if obj is None:
                    errors.append('%s is empty.' % obj)

        return errors

    def _summary(self):
        message = []
        for i, coordinates in enumerate(self.inputCoordinates):
            protocol = self.getMapper().getParent(coordinates.get())
            message.append("Method %d %s" % (i + 1, protocol.getClassLabel()))
        message.append("Radius = %d" % self.consensusRadius)
        message.append("Consensus = %d" % self.consensus)
        return message

    def _methods(self):
        return []

    def getMainInput(self):
        return self.inputCoordinates[0].get()

    @classmethod
    def getTomoId(self, fn):
        return int(removeBaseExt(fn).lstrip(self.FN_PREFIX))


def consensusWorker(coords, vesicles, trMats, consensus, consensusRadius, posFn, jaccFn=None,
                    mode=PICK_MODE_LARGER):
    """ Worker for calculate the consensus of N picking algorithms of
          M_n coordinates each one.

        coords: Array of N numpy arrays of M_n coordinates each one.
        consensus: Minimum number of votes to get a consensus coordinate
        consensusRadius: Tolerance to see two coordinates as the same (in pixels)
        posFn: Where to write the consensus coordinates
        jaccFN: Where to write the Jaccard index per micrograph
    """
    if len(coords) == 1:  # self consensus (remove duplicates)
        N0 = 0
        firstInput = 0
    else:  # in regular consensus all first coords are directly added to allCoords
        N0 = coords[0].shape[0]
        firstInput = 1


    # initializing arrays
    Ninputs = len(coords)
    Ncoords = sum([x.shape[0] for x in coords])
    allCoords = np.zeros([Ncoords, 3])
    allVesicles = np.zeros([Ncoords])
    allQuaternions = np.zeros([Ncoords, 4])
    votes = np.zeros(Ncoords)

    # Convert transformation matrices (rotations) to quaternions
    quaternions = []
    for n in range(len(trMats)):
        tomoMats = trMats[n]
        quatertions_tomo = np.asarray([quaternion_from_matrix(tomoMats[idq]) for idq in range(len(tomoMats))])
        quaternions.append(quatertions_tomo)

    inAllMicrographs = consensus <= 0 or consensus >= Ninputs

    # if nothing in the first and it should be in all, nothing to do
    if (not all([coords[idx].shape[0] for idx in range(Ninputs)])
            and inAllMicrographs):
        print("Returning from worker: doing AND consensus and, at least, one "
              "picker is empty for this micrograph (%s)." % posFn)
        return

    # Add all the first coordinates to 'allCoords' and 'votes' lists
    if N0 > 0:
        allCoords[0:N0, :] = coords[0]
        allVesicles[0:N0] = vesicles[0]
        allQuaternions[0:N0] = quaternions[0]
        votes[0:N0] = 1

    # Add the rest of coordinates to 'allCoords' and 'votes' lists
    Ncurrent = N0
    for n in range(firstInput, Ninputs):
        for idv in np.unique(vesicles[n]):
            coords_vesicle = coords[n][np.where(vesicles[n] == idv)]
            q_vesicle = quaternions[n][np.where(vesicles[n] == idv)]
            for coord, q in zip(coords_vesicle, q_vesicle):
                if Ncurrent > 0:
                    dist = np.sum((coord - allCoords[0:Ncurrent]) ** 2, axis=1)
                    imin = np.argmin(dist)
                    if sqrt(dist[imin]) < consensusRadius:
                        newCoord = (votes[imin] * allCoords[imin,] + coord) / (
                                    votes[imin] + 1)
                        allCoords[imin,] = newCoord
                        tuple_q = (allQuaternions[imin,], q)
                        T = weighted_tensor(tuple_q)
                        q_bar, _ = mean_quaternion(T)
                        allQuaternions[imin,] = q_bar
                        allVesicles[imin] = idv
                        votes[imin] += 1
                    else:
                        allCoords[Ncurrent, :] = coord
                        allVesicles[Ncurrent] = idv
                        allQuaternions[Ncurrent] = q
                        votes[Ncurrent] = 1
                        Ncurrent += 1
                else:
                    allCoords[Ncurrent, :] = coord
                    allVesicles[Ncurrent] = idv
                    allQuaternions[Ncurrent] = q
                    votes[Ncurrent] = 1
                    Ncurrent += 1

    # Select those in the consensus
    if consensus <= 0 or consensus > Ninputs:
        consensus = Ninputs
    elif not isinstance(consensus, int):
        consensus = consensus.get()

    if mode==PICK_MODE_LARGER:
        consensusCoords = allCoords[votes >= consensus, :]
        qConsensus = allQuaternions[votes >= consensus]
        vesiclesConsensus = allVesicles[votes >= consensus]
    else:
        consensusCoords = allCoords[votes == consensus, :]
        qConsensus = allQuaternions[votes == consensus]
        vesiclesConsensus = allVesicles[votes == consensus]

    # Convert consensus quaternions back to consensus transformations
    trConsensus = np.zeros([len(qConsensus), 4, 4])
    for idq, q in enumerate(qConsensus):
        tr = quaternion_matrix(q)
        trConsensus[idq,] = tr

    try:
        if jaccFn:
            jaccardIdx = float(len(consensusCoords)) / (
                    float(len(allCoords)) / Ninputs)
            # COSS: Possible problem with concurrent writes
            with open(jaccFn, "a") as fhJaccard:
                fhJaccard.write("%s \t %f\n" % (posFn, jaccardIdx))
    except Exception as exc:
        print("Some error occurred during Jaccard index calculation or "
              "writing it's file. Maybe a concurrence issue:\n%s" % exc)
    # Write the consensus file only if there
    # are some coordinates (size > 0)
    if consensusCoords.size:
        np.savetxt(posFn, consensusCoords)
        np.savetxt(removeExt(posFn) + '_vIds.txt', vesiclesConsensus)
        with open(removeExt(posFn) + '_trMats.txt', 'w') as outfile:
            outfile.write('# Array shape: {0}\n'.format(trConsensus.shape))
            for tr_slice in trConsensus:
                np.savetxt(outfile, tr_slice, fmt='%-7.2f')
                outfile.write('# New Transformation Matrix\n')


def getReadyTomos(coordSet):
    coorSet = SetOfCoordinates3D(filename=coordSet.getFileName())
    coorSet.loadAllProperties()
    setClosed = coorSet.isStreamClosed()
    coorSet.close()
    currentPickTomos = {tomoAgg["_volId"] for tomoAgg in
                        coordSet.aggregate(["MAX"], "_volId", ["_volId"])}
    return currentPickTomos, setClosed
