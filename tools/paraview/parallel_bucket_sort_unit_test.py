import unittest
import os
import random
from mpi4py import MPI
from functools import cmp_to_key
from parallel_bucket_sort import parallel_sort, flatten_list, sort_list
from grid2paraview import SpartaGridFile

GRID_FILE_200 = os.path.join(os.path.dirname(__file__),
    "tests/input_files/grid.200")

CIRCLE_GRID_FILE = os.path.join(os.path.dirname(__file__),
    "tests/input_files/circle_grid.txt")

class TestParallelBucketSort(unittest.TestCase):

    COMM = MPI.COMM_WORLD
    RANK = COMM.Get_rank()
    NUM_RANKS = COMM.Get_size()
    ROOT = 0
    MIN_NUM = 1
    MAX_NUM = 10000000

    def testSortDataInput(self):
        self.COMM.Barrier()
        with self.assertRaises(ValueError):
            parallel_sort(5)

    def testSortCompareInput(self):
        self.COMM.Barrier()
        with self.assertRaises(ValueError):
            parallel_sort([], 8.2)

    def testSortSize0(self):
        self.COMM.Barrier()
        data = self.generateData(0)
        self.checkResult(data, parallel_sort(data))

    def testSortSize1(self):
        self.COMM.Barrier()
        data = self.generateData(1)
        self.checkResult(data, parallel_sort(data))

    def testSortSize2(self):
        self.COMM.Barrier()
        data = self.generateData(2)
        self.checkResult(data, parallel_sort(data))

    def testSortSize10(self):
        self.COMM.Barrier()
        data = self.generateData(10)
        self.checkResult(data, parallel_sort(data))

    def testSortSize100(self):
        self.COMM.Barrier()
        data = self.generateData(100)
        self.checkResult(data, parallel_sort(data))

    def testSortSize1000(self):
        self.COMM.Barrier()
        data = self.generateData(1000)
        self.checkResult(data, parallel_sort(data))

    def testSortSize10000(self):
        self.COMM.Barrier()
        data = self.generateData(10000)
        self.checkResult(data, parallel_sort(data))

    def testSortSize100000(self):
        self.COMM.Barrier()
        data = self.generateData(100000)
        self.checkResult(data, parallel_sort(data))

    def testSortSize1000000(self):
        self.COMM.Barrier()
        data = self.generateData(1000000)
        self.checkResult(data, parallel_sort(data))

    def testSpartaGridFile200(self):
        self.COMM.Barrier()
        self.sortGridFile(GRID_FILE_200)

    def testSpartaCircleGridFile(self):
        self.COMM.Barrier()
        self.sortGridFile(CIRCLE_GRID_FILE)

    def sortGridFile(self, grid_file):
        sgf = SpartaGridFile(grid_file)
        sgf.set_iteration_start(self.RANK)
        sgf.set_iteration_skip(self.NUM_RANKS)
        data = [line for line in sgf]
        self.checkResult(data, parallel_sort(data,
            SpartaGridFile.compare_dashed_ids),
                SpartaGridFile.compare_dashed_ids)

    def checkResult(self, data, sorted_data, compare=None):
        self.COMM.Barrier()
        all_data = self.COMM.gather(data, root = self.ROOT)
        all_sorted_data = self.COMM.gather(sorted_data, root = self.ROOT)
        if self.RANK == self.ROOT:
            all_data = flatten_list(all_data)
            sort_list(all_data, compare)
            all_sorted_data = flatten_list(all_sorted_data)
            areEqual = all_data == all_sorted_data
        else:
            areEqual = None
        self.assertTrue(self.COMM.bcast(areEqual, root = self.ROOT))

    def generateData(self, size):
        return [random.randint(self.MIN_NUM, self.MAX_NUM)
            for i in range(self.findSizeForRank(size))]

    def findSizeForRank(self, size):
        count = int(size / self.NUM_RANKS)
        remainder = size % self.NUM_RANKS
        start = self.RANK * count + min(self.RANK, remainder)
        stop = (self.RANK + 1) * count + min(self.RANK + 1, remainder)
        return int(stop - start)

if __name__ == '__main__':
    unittest.main()
