from __future__ import print_function

from mpi4py import MPI
from grid2paraview import SpartaGridFile, read_grid_description_file
from parallel_bucket_sort import parallel_sort_to_file_buckets, \
        get_bucket_file_name_for_rank, is_rank_zero, get_comm_world, \
            get_rank, get_size, barrier, error_found_on_rank_zero
from sort_sparta_grid_file import get_grid_file_prefix, sort_grid_file_to_files
import argparse
import sys
import os

def main():
    args = parse_command_line()
    check_command_line(args)
    grid_desc = create_grid_description(args)
    program_data = distribute_program_data(args, grid_desc)
    cells = read_grid_file_cells(program_data)

    print("rank " + str(get_rank()) + " has " + str(len(cells)) + " cell(s)")

def parse_command_line():
    args = None
    if is_rank_zero():
        parser = argparse.ArgumentParser()
        parser.add_argument("sparta_grid_description_file",
            help="SPARTA grid description input file name")
        parser.add_argument("paraview_output_file",
            help="ParaView output file name")
        args = parser.parse_args()
    return args

def check_command_line(args):
    error_flag = False
    if is_rank_zero():
        if not os.path.isfile(args.sparta_grid_description_file):
            print("Unable to open SPARTA grid description file: ",
                args.sparta_grid_description_file)
            error_flag = True

        if os.path.isdir(args.paraview_output_file):
            print("ParaView output directory exists: ",
                args.paraview_output_file)
            error_flag = True

    if error_found_on_rank_zero(error_flag):
        sys.exit(1)

def create_grid_description(args):
    error_flag = False
    grid_desc = None
    if is_rank_zero():
        grid_desc = {}
        gdf = open(args.sparta_grid_description_file, "r")
        read_grid_description_file(gdf, grid_desc)
        gdf.close()

        if "dimension" not in grid_desc:
            print("Error: grid description file does not have a dimension statement: ",
                args.sparta_grid_description_file)
            error_flag = True
            
        if "create_box" not in grid_desc:
            print("Error: grid description file does not have a create_box statement: ",
                args.sparta_grid_description_file)
            error_flag = True

        if "read_grid" not in grid_desc:
            print("Error: grid description file does not have a read_grid statement: ",
                args.sparta_grid_description_file)
            error_flag = True

        if os.path.isfile(args.paraview_output_file + '.pvd'):
            print("ParaView output file exists: ",
                args.paraview_output_file + '.pvd')
            error_flag = True

        try:
            sgf = SpartaGridFile(grid_desc["read_grid"])
        except:
            error_flag = True

    if error_found_on_rank_zero(error_flag):
        sys.exit(1)

    if is_rank_zero():
        print("Processing " + str(sgf.number_of_cells) +\
            " cell(s) on " + str(get_size()) + " MPI rank(s)")
        os.mkdir(args.paraview_output_file)

    return grid_desc

def distribute_program_data(args, grid_desc):
    pd = {}
    if is_rank_zero():
        pd["paraview_output_file"] = args.paraview_output_file
        pd["grid_desc"] = grid_desc
    else:
        pd["paraview_output_file"] = None
        pd["grid_desc"] = None

    pd["paraview_output_file"] = get_comm_world().bcast(
        pd["paraview_output_file"], root = 0)
    pd["grid_desc"] = get_comm_world().bcast(pd["grid_desc"], root = 0)
    return pd

def read_grid_file_cells(program_data):
    grid_file = program_data["grid_desc"]["read_grid"]

    if not exist_grid_file_files(grid_file):
        sort_grid_file_to_files(grid_file)

    return read_grid_file_files(grid_file)

def exist_grid_file_files(grid_file):
    files_exist = True
    if is_rank_zero():
        prefix = get_grid_file_prefix(grid_file)
        for rank in range(get_size()):
            file_name = get_bucket_file_name_for_rank(rank, prefix)
            files_exist = files_exist and os.path.isfile(file_name)
    return get_comm_world().bcast(files_exist, root = 0)

def read_grid_file_files(grid_file):
    prefix = get_grid_file_prefix(grid_file)
    file_name = get_bucket_file_name_for_rank(get_rank(), prefix)
    cells = []
    count = 1
    if is_rank_zero():
        print("Started sorted grid file read")
    with open(file_name, 'r') as f:
        for cell in f:
            cells.append(cell.strip())
            if is_rank_zero() and count % 10000 == 0:
                print("Read " + str(count) + " cells")
            count += 1
    if is_rank_zero():
        print("Finished sorted grid file read")
    return cells

if __name__ == "__main__":
    main()
