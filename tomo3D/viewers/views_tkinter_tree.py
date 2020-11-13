# **************************************************************************
# *
# * Authors:     David Herreros (dherreros@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from multiprocessing import Process

from pyworkflow import utils as pwutils
from pyworkflow.gui.dialog import ToolbarListDialog
from pyworkflow.gui.tree import TreeProvider

from .viewer_triangulations import TriangulationPlot
from ..utils import delaunayTriangulation

from tomo.utils import extractVesicles, initDictVesicles

class Tomo3DTreeProvider(TreeProvider):
    """ Populate Tree from SetOfTomograms. """

    def __init__(self, tomoList):
        TreeProvider.__init__(self)
        self.tomoList = tomoList

    def getColumns(self):
        return [('Tomogram', 300)]

    def getObjectInfo(self, tomo):
        tomogramName = pwutils.removeBaseExt(tomo.getFileName())

        return {'key': tomogramName, 'parent': None,
                'text': tomogramName,
                'tags': ("done")}

    def getObjectPreview(self, obj):
        return (None, None)

    def getObjectActions(self, obj):
        return []

    def _getObjectList(self):
        """Retrieve the object list"""
        return self.tomoList

    def getObjects(self):
        objList = self._getObjectList()
        return objList

    def configureTags(self, tree):
        tree.tag_configure("done", foreground="black")


class Tomo3DDialog(ToolbarListDialog):
    """
    This class extend from ListDialog to allow calling
    a pyvista viewer subprocess from a list of Tomograms.
    """

    def __init__(self, parent, coordinates, **kwargs):
        self.dictVesicles, _ = initDictVesicles(coordinates)
        self.coordinates = coordinates
        self.provider = kwargs.get("provider", None)
        ToolbarListDialog.__init__(self, parent,
                                   "Tomogram List",
                                   allowsEmptySelection=False,
                                   itemDoubleClick=self.doubleClickOnTomogram,
                                   **kwargs)

    def doubleClickOnTomogram(self, e=None):
        self.tomo = e
        tomoName = pwutils.removeBaseExt(self.tomo.getFileName())
        self.dictVesicles = extractVesicles(self.coordinates, self.dictVesicles, tomoName)
        self.proc = Process(target=createViewer, args=(self.tomo, self.dictVesicles))
        self.proc.start()

def createViewer(tomo, vesicles_dict):
    tomoName = pwutils.removeBaseExt(tomo.getFileName())
    normals = vesicles_dict[tomoName]['normals']
    vesicles = vesicles_dict[tomoName]['vesicles']
    shells = []
    for vesicle in vesicles:
        shells.append(delaunayTriangulation(vesicle))
    plt = TriangulationPlot(shells, clouds=vesicles, extNormals_List=normals)
    plt.initializePlot()
