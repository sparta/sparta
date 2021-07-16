from __future__ import print_function

from mpi4py import MPI
from grid2paraview import SpartaGridFile, read_grid_description_file
from parallel_bucket_sort import parallel_sort
import argparse
import sys
import os

def read_grid_file(pd):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    num_ranks = comm.Get_size()

    sgf = SpartaGridFile(pd["grid_desc"]["read_grid"])
    sgf.set_iteration_start(rank)
    sgf.set_iteration_skip(num_ranks)
    if rank == 0:
        print("Reading Sparta grid file " + pd["grid_desc"]["read_grid"])
    cells = []
    count = 0
    for cell in sgf:
        cells.append(cell)
        count += 1
        if rank == 0 and count % 10000 == 0:
            print("Read " + str(count) + " cell(s) from grid file")
    if rank == 0:
        print("Finished grid file read")
        print("Starting parallel sort")
    cells = parallel_sort(cells, SpartaGridFile.compare_dashed_ids,
        use_file_buckets=True)
    if rank == 0:
        print("Finished parallel sort")

    print("rank " + str(rank) + " has " + str(len(cells)) + " cell(s)")

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    num_ranks = comm.Get_size()

    if rank == 0:
        parser = argparse.ArgumentParser()
        parser.add_argument("sparta_grid_description_file", help="SPARTA grid description input file name")
        parser.add_argument("paraview_output_file", help="ParaView output file name")
        args = parser.parse_args()

        try:
            gdf = open(args.sparta_grid_description_file, "r")
        except IOError:
            print("Unable to open SPARTA surf input file: ", args.sparta_grid_description_file)
            sys.exit(1)

        if os.path.isdir(args.paraview_output_file):
            print("ParaView output directory exists: ", args.paraview_output_file)
            sys.exit(1)

        grid_desc = {}
        read_grid_description_file(gdf, grid_desc)
        gdf.close()

        if "dimension" not in grid_desc:
            print("Error: grid description file does not have a dimension statement: ",
                args.sparta_grid_description_file)
            sys.exit(1)

        if "create_box" not in grid_desc:
            print("Error: grid description file does not have a create_box statement: ",
                args.sparta_grid_description_file)
            sys.exit(1)

        if "read_grid" not in grid_desc:
            print("Error: grid description file does not have a read_grid statement: ",
                args.sparta_grid_description_file)
            sys.exit(1)

        if os.path.isfile(args.paraview_output_file + '.pvd'):
            print("ParaView output file exists: ",
                args.paraview_output_file + '.pvd')
            sys.exit(1)

        sgf = SpartaGridFile(grid_desc["read_grid"])
        print("Processing " + str(sgf.number_of_cells) +\
            " cell(s) on " + str(num_ranks) + " MPI rank(s)")
        os.mkdir(args.paraview_output_file)

    pd = {}
    if rank == 0:
        pd["paraview_output_file"] = args.paraview_output_file
        pd["grid_desc"] = grid_desc
    else:
        pd["paraview_output_file"] = None
        pd["grid_desc"] = None

    pd["paraview_output_file"] = comm.bcast(pd["paraview_output_file"], root=0)
    pd["grid_desc"] = comm.bcast(pd["grid_desc"], root=0)

    read_grid_file(pd)
