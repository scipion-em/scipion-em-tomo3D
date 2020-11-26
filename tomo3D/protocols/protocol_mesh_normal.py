# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
# *
# *  BCU, Centro Nacional de Biotecnologia, CSIC
# *
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from os import path
import numpy as np
from pyworkflow.protocol.params import PointerParam, FloatParam, LEVEL_ADVANCED, BooleanParam, IntParam, EnumParam
import pyworkflow.utils as pwutlis
import pwem.convert.transformations as tfs
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
from tomo.utils import normalFromMatrix
from ..utils import delaunayTriangulation, computeNormals


class XmippProtFilterbyNormal(EMProtocol, ProtTomoBase):
    """ This protocol takes surfaces or ROIs (SetOfMeshes) and a SetOfSubtomograms and filters them by different
    criteria related with the normal direction."""

    _label = 'filter subtomos by normal'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSubtomos', PointerParam, pointerClass="SetOfSubTomograms",
                      label='Subtomograms', help='SetOfSubtomograms to filter.')
        form.addParam('normalDir', BooleanParam, default=True,
                      label='Filter subtomograms by normal',
                      help='Remove the subtomograms that have a normal direction not equal to the normal direction of '
                           'the vesicle in the coordinate of the particle.')
        form.addParam('inputMeshes', PointerParam, label="Vesicles", pointerClass='SetOfMeshes', condition='normalDir',
                      help='Select the vesicles in which the subtomograms are.')
        form.addParam('tol', FloatParam, default=5, label='Tolerance in degrees',
                      condition='normalDir', expertLevel=LEVEL_ADVANCED,
                      help='Tolerance (in degrees) when comparing between particle and mesh normal directions.')
        form.addParam('tilt', BooleanParam, default=False,
                      label='Filter subtomograms by tilt angle',
                      help='Remove subtomograms depending on their tilt angle.')
        form.addParam('maxtilt', IntParam, default=150, label='Maximum allowed tilt', condition='tilt',
                      help='Remove the subtomograms that have a tilt angle bigger than the one specified in here, '
                           'considering tilt angle between 0 and 180 degrees.')
        form.addParam('mintilt', IntParam, default=30, label='Minimum allowed tilt', condition='tilt',
                      help='Remove the subtomograms that have a tilt angle smaller than the one specified in here, '
                           'considering tilt angle between 0 and 180 degrees.')

        # form.addParam('topBottom', BooleanParam, default=True,
        #               label='Remove subtomograms in the top and bottom of the vesicle',
        #               help='Remove the subtomograms that have been picked from the top and bottom parts because they '
        #                    'had a different view.')
        # form.addParam('mwDir', BooleanParam, default=True,
        #               label='Remove subtomograms in the missing wedge direction',
        #               help='Remove the subtomograms that are in the missing wedge direction because they are highly '
        #                    'affected by the missing wedge.')
        # form.addParam('mwDir', EnumParam, default=True, label='Missing wedge direction',
        #               help='Missing wedge direction of the tomograms.')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('computeNormalStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -------------------------------
    def computeNormalStep(self):
        inSet = self.inputSubtomos.get()
        self.outSet = self._createSetOfSubTomograms()
        self.outSet.copyInfo(inSet)

        if self.normalDir.get():
            tol = self.tol.get()*np.pi/180
            for subtomo in inSet:
                for mesh in self.inputMeshes.get().iterItems():
                    pathV = pwutlis.removeBaseExt(path.basename(mesh.getPath())).split('_vesicle_')
                    if pwutlis.removeBaseExt(path.basename(subtomo.getVolName())) == pathV[0]:
                        if str(self._getVesicleId(subtomo)) == pathV[1]:
                            normalsList = self._getNormalVesicleList(mesh)
                            normSubtomo, normVesicle = self._getNormalVesicle(normalsList, subtomo)
                            if abs(normSubtomo[0]-normVesicle[0]) < tol and abs(normSubtomo[1]-normVesicle[1]) < tol \
                                    and abs(normSubtomo[2]-normVesicle[2]) < tol:
                                if self.tilt.get():
                                    tilt = self._getTiltSubtomo(subtomo)
                                    if self.maxtilt.get() > tilt > self.mintilt.get():
                                        self.outSet.append(subtomo)
                                else:
                                    self.outSet.append(subtomo)

        if self.tilt.get() and not self.normalDir.get():
            for subtomo in inSet:
                tilt = self._getTiltSubtomo(subtomo)
                if self.maxtilt.get() > tilt > self.mintilt.get():
                    self.outSet.append(subtomo)

    def createOutputStep(self):
        self._defineOutputs(outputset=self.outSet)
        self._defineSourceRelation(self.inputSubtomos.get(), self.outSet)

    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        validateMsgs = []
        if not self.normalDir.get() and not self.tilt.get():
            validateMsgs.append('Some filter should be switched to "Yes"')
        return validateMsgs

    def _summary(self):
        summary = []
        if not self.isFinished():
            summary.append("Output subtomograms not ready yet.")
        else:
            if self.normalDir:
                summary.append("Remove subtomograms by normal direction (tolerance of %0.2f degrees)" %
                               self.tol.get())
            if self.tilt:
                summary.append("Remove subtomograms by tilt angle (max allowed tilt: %d, min allowed tilt: %d)" %
                               (self.maxtilt.get(), self.mintilt.get()))
        return summary

    def _methods(self):
        methods = []
        if not self.isFinished():
            methods.append("Output subtomograms not ready yet.")
        else:
            if self.normalDir:
                methods.append("Subtomograms that are not perpendicular to the membrane have been removed.")
            if self.tilt:
                methods.append("Subtomograms with tilt angle bigger than %d or smaller than %d have been removed." %
                               (self.maxtilt.get(), self.mintilt.get()))
            # if self.topBottom:
            #     methods.append("Particles in the top and bottom parts of the vesicles have been removed.")
            # if self.mwDir:
            #     methods.append("Particles in the missing wedge direction have been removed.")
        return methods

    # --------------------------- UTILS functions --------------------------------------------
    def _getVesicleId(self, subtomo):
        c = subtomo.getCoordinate3D()
        if c.hasGroupId():
            vesicleId = c.getGroupId()
        else:  # For now it works with several vesicles in the same tomo just for Pyseg subtomos
            vesicleId = 1
        return vesicleId

    def _getNormalVesicleList(self, mesh):
        triangulation = delaunayTriangulation(mesh.getMesh())
        normalsList = computeNormals(triangulation, associateCoords=True)
        return normalsList

    def _getNormalVesicle(self, normalsList, subtomo):
        normSubtomo = normalFromMatrix(subtomo.getTransform().getMatrix())
        coord = subtomo.getCoordinate3D()
        coors = np.asarray([coord.getX(), coord.getY(), coord.getZ()])
        points, normals = zip(*normalsList)
        points = np.asarray(points)
        idx = np.argmin(np.sum((points - coors) ** 2, axis=1))
        return normSubtomo, normals[idx]

    def _getTiltSubtomo(self, subtomo):
        _, tilt, _ = tfs.euler_from_matrix(subtomo.getTransform().getMatrix(), axes='szyz')
        tilt = -np.rad2deg(tilt)
        return tilt
