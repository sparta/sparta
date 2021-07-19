from __future__ import print_function

from mpi4py import MPI
import argparse
import sys
import os
from grid2paraview import SpartaGridFile
from parallel_bucket_sort import parallel_sort_to_file_buckets, \
    is_rank_zero, get_comm_world, get_rank, get_size, error_found_on_rank_zero

def main():
    args = parse_command_line()
    check_grid_file(args)
    grid_file_path = get_grid_file_path(args)
    sort_grid_file_to_files(grid_file_path)

def parse_command_line():
    if is_rank_zero():
        parser = argparse.ArgumentParser()
        parser.add_argument("sparta_grid_file_path", help="Sparta grid file path")
        args = parser.parse_args()
        return args
    else:
        return None

def check_grid_file(args):
    error_flag = False
    if is_rank_zero():
        try:
            sgf = SpartaGridFile(args.sparta_grid_file_path)
        except:
            error_flag = True

    if error_found_on_rank_zero(error_flag):
        sys.exit(1)

def sort_grid_file_to_files(grid_file_path):
    cells = get_grid_file_cells(grid_file_path)
    prefix = get_grid_file_prefix(grid_file_path)
    if is_rank_zero():
        print("Started parallel sort")
    parallel_sort_to_file_buckets(cells,
        SpartaGridFile.compare_dashed_ids, prefix)
    if is_rank_zero():
        print("Finished parallel sort")

def get_grid_file_path(args):
    path = None
    if is_rank_zero():
        path = args.sparta_grid_file_path
        path = get_comm_world().bcast(path, root = 0)
    else:
        path = get_comm_world().bcast(path, root = 0)
    return path

def get_grid_file_cells(grid_file_path):
    sgf = SpartaGridFile(grid_file_path)
    sgf.set_iteration_start(get_rank())
    sgf.set_iteration_skip(get_size())
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

if __name__ == "__main__":
    main()
