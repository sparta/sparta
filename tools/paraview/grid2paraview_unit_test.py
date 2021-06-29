import unittest
import os
from grid2paraview import SpartaGridFile

GRID_FILE_200 = os.path.join(os.path.dirname(__file__),
    "tests/input_files/grid.200")

CIRCLE_GRID_FILE = os.path.join(os.path.dirname(__file__),
    "tests/input_files/circle_grid.txt")

class TestSpartaGridFile(unittest.TestCase):

    def testSpartaGridFile200(self):
       sgf = SpartaGridFile(GRID_FILE_200)
       self.checkDashedIdConversion(sgf)

    def testSpartaCircleGridFile(self):
       cgf = SpartaGridFile(CIRCLE_GRID_FILE)
       self.checkDashedIdConversion(cgf)

    def checkDashedIdConversion(self, sparta_grid_file):
       sparta_grid_file.iterate_local_cell_ids = True
       for local_cell_id in sparta_grid_file:
           dashed_id = sparta_grid_file.create_dashed_id(local_cell_id)
           new_local_cell_id = sparta_grid_file\
               .get_local_cell_id_from_dashed_cell_id(dashed_id)
           if dashed_id:
             self.assertEqual(local_cell_id, str(new_local_cell_id))

if __name__ == '__main__':
    unittest.main()
