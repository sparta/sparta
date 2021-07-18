from __future__ import print_function

from mpi4py import MPI
import argparse
import sys
import os
from grid2paraview import SpartaGridFile
from parallel_bucket_sort import parallel_sort_to_file_buckets

def main():
    args = parse_command_line()
    check_grid_file(args)
    sort_grid_file_to_files(args)

def parse_command_line():
    if is_rank_zero():
        parser = argparse.ArgumentParser()
        parser.add_argument("sparta_grid_file_path", help="Sparta grid file path")
        args = parser.parse_args()
        return args
    else:
        return None

def check_grid_file(args):
    COMM = MPI.COMM_WORLD
    error_detected = False
    if is_rank_zero():
        try:
            sgf = SpartaGridFile(args.sparta_grid_file_path)
        except:
            error_detected = True
        error_detected = COMM.bcast(error_detected, root = 0)
    else:
        error_detected = COMM.bcast(error_detected, root = 0)
         
    if error_detected:
        sys.exit(1)

def sort_grid_file_to_files(args):
    grid_file_path = get_grid_file_path(args)
    cells = get_grid_file_cells(grid_file_path)
    prefix = get_grid_file_prefix(grid_file_path)
    if is_rank_zero():
        print("Started parallel sort")
    parallel_sort_to_file_buckets(cells,
        SpartaGridFile.compare_dashed_ids, prefix)
    if is_rank_zero():
        print("Finished parallel sort")

def get_grid_file_path(args):
    COMM = MPI.COMM_WORLD
    path = None
    if is_rank_zero():
        path = args.sparta_grid_file_path
        path = COMM.bcast(path, root = 0)
    else:
        path = COMM.bcast(path, root = 0)
    return path

def get_grid_file_cells(grid_file_path):
    COMM = MPI.COMM_WORLD
    sgf = SpartaGridFile(grid_file_path)
    sgf.set_iteration_start(COMM.Get_rank())
    sgf.set_iteration_skip(COMM.Get_size())
    if is_rank_zero():
        print("Reading Sparta grid file " + grid_file_path)
    cells = []
    count = 0
    for cell in sgf:
        cells.append(cell)
        count += 1
        if is_rank_zero() and count % 10000 == 0:
            print("Read " + str(count) + " cell(s) from grid file")
    if is_rank_zero():
        print("Finished grid file read")
    return cells

def get_grid_file_prefix(grid_file_path):
    no_suffix = os.path.splitext(grid_file_path)[0]
    return os.path.basename(no_suffix)

def is_rank_zero():
    COMM = MPI.COMM_WORLD
    return COMM.Get_rank() == 0

if __name__ == "__main__":
    main()
