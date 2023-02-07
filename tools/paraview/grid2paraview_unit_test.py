import unittest
import os
from grid2paraview import SpartaGridFile
from functools import cmp_to_key

GRID_FILE_200 = os.path.join(os.path.dirname(__file__),
    "grid.200")

CIRCLE_GRID_FILE = os.path.join(os.path.dirname(__file__),
    "circle_grid.txt")

class TestSpartaGridFile(unittest.TestCase):

    def testSpartaGridFile200(self):
        sgf = SpartaGridFile(GRID_FILE_200)
        self.checkDashedIdConversion(sgf)

    def testSpartaCircleGridFile(self):
        cgf = SpartaGridFile(CIRCLE_GRID_FILE)
        self.checkDashedIdConversion(cgf)

    def testSpartaGridFile200IterationSkip(self):
        no_skip = self.readGridFileWithSkip(GRID_FILE_200, 0, 1)
        no_skip.sort(key=cmp_to_key(SpartaGridFile.compare_dashed_ids))
        start_zero = self.readGridFileWithSkip(GRID_FILE_200, 0, 2)
        start_one = self.readGridFileWithSkip(GRID_FILE_200, 1, 2)
        skipped = start_zero + start_one
        skipped.sort(key=cmp_to_key(SpartaGridFile.compare_dashed_ids))
        self.assertTrue(no_skip == skipped)

    def testSpartaCircleGridFileIterationSkip(self):
        no_skip = self.readGridFileWithSkip(CIRCLE_GRID_FILE, 0, 1)
        no_skip.sort(key=cmp_to_key(SpartaGridFile.compare_dashed_ids))
        start_zero = self.readGridFileWithSkip(CIRCLE_GRID_FILE, 0, 4)
        start_one = self.readGridFileWithSkip(CIRCLE_GRID_FILE, 1, 4)
        start_two = self.readGridFileWithSkip(CIRCLE_GRID_FILE, 2, 4)
        start_three = self.readGridFileWithSkip(CIRCLE_GRID_FILE, 3, 4)
        skipped = start_zero + start_one + start_two + start_three
        skipped.sort(key=cmp_to_key(SpartaGridFile.compare_dashed_ids))
        self.assertTrue(no_skip == skipped)

    def readGridFileWithSkip(self, grid_file, start, skip):
        sgf = SpartaGridFile(grid_file)
        sgf.set_iteration_start(start)
        sgf.set_iteration_skip(skip)
        return [line for line in sgf]

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
