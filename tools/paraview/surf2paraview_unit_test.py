import unittest
import os
import vtk
from surf2paraview import read_surf_file

CIRCLE_DATA_W_TYPE = os.path.join(os.path.dirname(__file__),
    "data.circle_w_type")

CIRCLE_DATA_WO_TYPE = os.path.join(os.path.dirname(__file__),
    "../examples/circle/data.circle")

SPHERE_DATA_W_TYPE = os.path.join(os.path.dirname(__file__),
    "data.sphere_w_type")

SPHERE_DATA_WO_TYPE = os.path.join(os.path.dirname(__file__),
    "../examples/sphere/data.sphere")

class TestSurf2ParaView(unittest.TestCase):

    def testReadCircleDataNoType(self):
        sif = open(CIRCLE_DATA_WO_TYPE, "r")
        self.readFileAndTest(sif, 50, 50)

    def testReadCircleDataWithType(self):
        sif = open(CIRCLE_DATA_W_TYPE, "r")
        self.readFileAndTest(sif, 50, 100)

    def testReadSphereDataNoType(self):
        sif = open(SPHERE_DATA_WO_TYPE, "r")
        self.readFileAndTest(sif, 1200, 602)

    def testReadSphereDataWithType(self):
        sif = open(SPHERE_DATA_W_TYPE, "r")
        self.readFileAndTest(sif, 1200, 3600)

    def readFileAndTest(self, sif, num_cells, num_points):
        ug = vtk.vtkUnstructuredGrid()
        read_surf_file(sif, ug)
        sif.close()
        self.assertEqual(ug.GetNumberOfCells(), num_cells)
        self.assertEqual(ug.GetNumberOfPoints(), num_points)

if __name__ == '__main__':
    unittest.main()
