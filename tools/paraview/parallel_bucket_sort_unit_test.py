import unittest
import os
import random
from mpi4py import MPI
from functools import cmp_to_key
from parallel_bucket_sort import parallel_sort, \
    flatten_list, sort_list, parallel_sort_to_file_buckets, \
        remove_file_buckets, get_bucket_file_name_for_rank, \
            is_rank_zero, get_rank, get_size, barrier, get_comm_world
from grid2paraview import SpartaGridFile

GRID_FILE_200 = os.path.join(os.path.dirname(__file__),
    "grid.200")

CIRCLE_GRID_FILE = os.path.join(os.path.dirname(__file__),
    "circle_grid.txt")

class TestParallelBucketSort(unittest.TestCase):

    MIN_NUM = 1
    MAX_NUM = 100000

    @classmethod
    def tearDownClass(cls):
        barrier()
        MPI.Finalize()

    def testSortDataInput(self):
        barrier()
        with self.assertRaises(ValueError):
            parallel_sort(5)

    def testSortCompareInput(self):
        barrier()
        with self.assertRaises(ValueError):
            parallel_sort([], 8.2)

    def testSortSize0(self):
        barrier()
        data = self.generateData(0)
        self.checkResult(data, parallel_sort(data))

    def testSortSize1(self):
        barrier()
        data = self.generateData(1)
        self.checkResult(data, parallel_sort(data))

    def testSortSize2(self):
        barrier()
        data = self.generateData(2)
        self.checkResult(data, parallel_sort(data))

    def testSortSize10(self):
        barrier()
        data = self.generateData(10)
        self.checkResult(data, parallel_sort(data))

    def testSortSize100(self):
        barrier()
        data = self.generateData(100)
        self.checkResult(data, parallel_sort(data))

    def testSortSize1000(self):
        barrier()
        data = self.generateData(1000)
        self.checkResult(data, parallel_sort(data))

    def testSortSize10000(self):
        barrier()
        data = self.generateData(10000)
        self.checkResult(data, parallel_sort(data))

    def testSortSize100000(self):
        barrier()
        data = self.generateData(100000)
        self.checkResult(data, parallel_sort(data))

    def testSortSize1000000(self):
        barrier()
        data = self.generateData(1000000)
        self.checkResult(data, parallel_sort(data))

    def testSpartaGridFile200(self):
        barrier()
        self.sortGridFile(GRID_FILE_200)
        self.sortGridFile(GRID_FILE_200, use_file_buckets=True)

    def testSpartaCircleGridFile(self):
        barrier()
        self.sortGridFile(CIRCLE_GRID_FILE)
        self.sortGridFile(CIRCLE_GRID_FILE, use_file_buckets=True)

    def testSpartaGridFile200SortToFileBuckets(self):
        barrier()
        self.sortGridFileToFileBuckets(GRID_FILE_200)

    def testSpartaCircleGridFileSortToFileBuckets(self):
        barrier()
        self.sortGridFileToFileBuckets(CIRCLE_GRID_FILE)

    def sortGridFileToFileBuckets(self, grid_file):
        data = self.loadDataFromGridFile(grid_file)
        parallel_sort_to_file_buckets(data,
            SpartaGridFile.compare_dashed_ids, "test")
        filename = get_bucket_file_name_for_rank(get_rank(), "test")
        result = []
        with open(filename, 'r') as f:
            for line in f:
                result.append(line.strip())
        sort_list(result, SpartaGridFile.compare_dashed_ids)
        self.checkResult(data, result, SpartaGridFile.compare_dashed_ids)
        remove_file_buckets("test")

    def sortGridFile(self, grid_file, use_file_buckets=False):
        data = self.loadDataFromGridFile(grid_file)
        self.checkResult(data, parallel_sort(data,
            SpartaGridFile.compare_dashed_ids, use_file_buckets),
                SpartaGridFile.compare_dashed_ids)

    def loadDataFromGridFile(self, grid_file):
        sgf = SpartaGridFile(grid_file)
        sgf.set_iteration_start(get_rank())
        sgf.set_iteration_skip(get_size())
        return [line for line in sgf]

    def checkResult(self, data, sorted_data, compare=None):
        barrier()
        all_data = get_comm_world().gather(data, root = 0)
        all_sorted_data = get_comm_world().gather(sorted_data, root = 0)
        if is_rank_zero():
            all_data = flatten_list(all_data)
            sort_list(all_data, compare)
            all_sorted_data = flatten_list(all_sorted_data)
            areEqual = all_data == all_sorted_data
        else:
            areEqual = None
        areEqual = get_comm_world().bcast(areEqual, 0)
        self.assertTrue(areEqual)

    def generateData(self, size):
        return [random.randint(self.MIN_NUM, self.MAX_NUM)
            for i in range(self.findSizeForRank(size))]

    def findSizeForRank(self, size):
        count = int(size / get_size())
        remainder = size % get_size()
        start = get_rank() * count + min(get_rank(), remainder)
        stop = (get_rank() + 1) * count + min(get_rank() + 1, remainder)
        return int(stop - start)

if __name__ == '__main__':
    unittest.main()
