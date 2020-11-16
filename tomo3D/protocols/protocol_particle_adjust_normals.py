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


class ProtTomoPickingAdjustNormals(ProtTomoPicking):
    """
    This protocol can be used to readjust the normals (Transformation Matrices) associated to a 
    SetOfCoordinates3D or to add them to vectorize a picking
    """

    _label = 'adjust normals'
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
            for idn in range(len(self.tomo_vesicles[tomoName]['vesicles'])):
                shell = delaunayTriangulation(self.tomo_vesicles[tomoName]['vesicles'][idn])
                self.tomo_vesicles[tomoName]['normals'][idn] = computeNormals(shell)
    
    def createOutputStep(self, coordinates):
        outSet = self._createSetOfCoordinates3D(coordinates)
        outSet.setPrecedents(coordinates.getPrecedents())
        outSet.setBoxSize(coordinates.getBoxSize())
        outSet.setSamplingRate(coordinates.getSamplingRate())
        for tomoName in self.tomoNames:
            ids = self.tomo_vesicles[tomoName]['ids']
            normals = self.tomo_vesicles[tomoName]['normals']
            for idv in range(len(ids)):
                for idp in range(len(ids[idv])):
                    outCoord = coordinates[int(ids[idv][idp])].clone()
                    outCoord.setBoxSize(outSet.getBoxSize())
                    tr_Mat = rotation_matrix_from_vectors(normals[idv][idp], np.array([0, 0, 1]))
                    outCoord.setMatrix(tr_Mat)
                    outSet.append(outCoord)
        self._defineOutputs(outputCoordinates=outSet)
        self._defineSourceRelation(coordinates, outSet)
