# **************************************************************************
# *
# * Authors:    David Herreros Calero (dherreros@cnb.csic.es) -- Tomo version
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


import numpy as np

import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from tomo.objects import Coordinate3D

from tomo.protocols import ProtTomoPicking
from tomo.utils import initDictVesicles, extractVesicles

from ..utils import delaunayTriangulation, computeNormals, rotation_matrix_from_vectors


class ProtTomoPickingReadjustNormals(ProtTomoPicking):
    """
    This protocol can be used to readjust the normals (Transformation Matrices) associated to a 
    SetOfCoordinates3D or to add them to vectorize a picking
    """

    _label = 'readjust normlas'
    outputName = 'outputCoordinates'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      label="Input coordinates", important=True,
                      help='Select the set of 3D coordinates to perform the normal adjustemnt')

#--------------------------- INSERT steps functions ----------------------------
    def _insertAllSteps(self):
        coordinates = self.inputCoordinates.get()
        self._insertFunctionStep('computeParams', coordinates)
        self._insertFunctionStep('normalAdjustment')
        self._insertFunctionStep('createOutputStep', coordinates)
        
    def computeParams(self, coordinates):
        self.tomo_vesicles, self.tomoNames = initDictVesicles(coordinates)
        for tomoName in self.tomoNames:
            self.tomo_vesicles = extractVesicles(coordinates, self.tomo_vesicles, tomoName)
    
    def normalAdjustment(self):
        for tomoName in self.tomoNames:
            for idn in  len(self.tomo_vesicles[tomoName]['vesicles']):
                shell = delaunayTriangulation(self.tomo_vesicles[tomoName]['vesicles'][idn])
                self.tomo_vesicles[tomoName]['normals'][idn] = computeNormals(shell)
    
    def createOutputStep(self, coordinates):
        outSet = self._createSetOfCoordinates3D(coordinates)
        outSet.setPrecedents(coordinates.getPrecedents())
        outSet.setBoxSize(coordinates.getBoxSize())
        outSet.setSamplingRate(coordinates.getSamplingRate())
        volIds = coordinates.aggregate(["MAX"], "_volId", ["_volId"])
        volIds = [d['_volId'] for d in volIds]
        for tomoName in self.tomoNames:
            for idn in  len(self.tomo_vesicles[tomoName]['vesicles']):
                outCoord = coordinates[self.tomo_vesicles[tomoName]['ids'][idn]].clone()
                outCoord.setBoxSize(outSet.getBoxSize())
                tr_Mat = rotation_matrix_from_vectors(self.tomo_vesicles[tomoName]['normals'][idn], 
                                                      np.array([0, 0, 1]))
                outCoord.setMatrix(tr_Mat)
                outSet.append(outCoord)
