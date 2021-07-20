from __future__ import print_function
from __future__ import division

from mpi4py import MPI
from grid2paraview import SpartaGridFile, read_grid_description_file
from parallel_bucket_sort import parallel_sort_to_file_buckets, \
        get_bucket_file_name_for_rank, is_rank_zero, get_comm_world, \
            get_rank, get_size, barrier, error_found_on_rank_zero
from sort_sparta_grid_file import get_grid_file_prefix, sort_grid_file_to_files
import argparse
import sys
import os
import vtk

def main():
    args = parse_command_line()
    check_command_line(args)
    grid_desc = create_grid_description(args)
    program_data = distribute_program_data(args, grid_desc)
    unstructured_grid = create_grid(program_data)

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(unstructured_grid)
    writer.SetFileName("grid_" + str(get_rank()) + ".vtu")
    writer.Write()

    print("rank " + str(get_rank()) + " has " + \
        str(unstructured_grid.GetNumberOfCells()) + " cell(s)")

def parse_command_line():
    args = None
    if is_rank_zero():
        parser = argparse.ArgumentParser()
        parser.add_argument("sparta_grid_description_file",
            help="SPARTA grid description input file name")
        parser.add_argument("paraview_output_file",
            help="ParaView output file name")
        group = parser.add_mutually_exclusive_group()
        group.add_argument('-r', '--result',
            help="Optional list of SPARTA dump result files", nargs='+')
        group.add_argument('-rf', '--resultfile',
            help="Optional filename containing path names of SPARTA dump result files")
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

def create_grid(program_data):
    grid_file = program_data["grid_desc"]["read_grid"]
    if not exist_grid_file_files(program_data):
        sort_grid_file_to_files(grid_file)

    return create_grid_from_grid_files(program_data)

def exist_grid_file_files(program_data):
    files_exist = True
    if is_rank_zero():
        for rank in range(get_size()):
            files_exist = files_exist and \
                os.path.isfile(get_bucket_file_name(rank, program_data))
    return get_comm_world().bcast(files_exist, root = 0)

def create_grid_from_grid_files(program_data):
    grid_file = program_data["grid_desc"]["read_grid"]
    sgf = SpartaGridFile(grid_file)
    append = vtk.vtkAppendFilter()
    append.MergePointsOn()
    count = 1
    if is_rank_zero():
        print("Started sorted grid file read")
        print("Creating grid")
    with open(get_bucket_file_name(get_rank(), program_data), 'r') as f:
        for cell in f:
            ug_cell = create_unstructured_grid_cell(
                cell.strip(), sgf, program_data)
            append.AddInputData(ug_cell)
            if is_rank_zero() and count % 10000 == 0:
                print("Read " + str(count) + " cells")
            count += 1
    if is_rank_zero():
        print("Finished sorted grid file read")
        print("Finished creating grid")
    append.Update()
    return append.GetOutput()

def create_unstructured_grid_cell(cell_dashed_id, sgf, program_data):
    if program_data['grid_desc']['dimension'] == 3:
        return create_unstructured_grid_hex(cell_dashed_id, sgf, program_data)
    else:
        return create_unstructured_grid_quad(cell_dashed_id, sgf, program_data)

def create_unstructured_grid_hex(cell_dashed_id, sgf, program_data):
    sgc = SpartaGridCell(cell_dashed_id, sgf, program_data)
    ug = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    points.InsertNextPoint(sgc.origin_x, sgc.origin_y, sgc.origin_z)
    points.InsertNextPoint(sgc.origin_x + sgc.length_x,
        sgc.origin_y, sgc.origin_z)
    points.InsertNextPoint(sgc.origin_x + sgc.length_x,
        sgc.origin_y + sgc.length_y, sgc.origin_z)
    points.InsertNextPoint(sgc.origin_x,
        sgc.origin_y + sgc.length_y, sgc.origin_z)
    points.InsertNextPoint(sgc.origin_x,
        sgc.origin_y, sgc.origin_z + sgc.length_z)
    points.InsertNextPoint(sgc.origin_x + sgc.length_x,
        sgc.origin_y, sgc.origin_z + sgc.length_z)
    points.InsertNextPoint(sgc.origin_x + sgc.length_x,
        sgc.origin_y + sgc.length_y, sgc.origin_z + sgc.length_z)
    points.InsertNextPoint(sgc.origin_x,
        sgc.origin_y + sgc.length_y, sgc.origin_z + sgc.length_z)
    ug.SetPoints(points)

    hex = vtk.vtkHexahedron()
    for i in range(8):
        hex.GetPointIds().SetId(i, i)

    ug.InsertNextCell(hex.GetCellType(), hex.GetPointIds())
    add_global_id(cell_dashed_id, sgf, ug)
    return ug

def create_unstructured_grid_quad(cell_dashed_id, sgf, program_data):
    sgc = SpartaGridCell(cell_dashed_id, sgf, program_data)
    ug = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    points.InsertNextPoint(sgc.origin_x, sgc.origin_y, 0.0)
    points.InsertNextPoint(sgc.origin_x + sgc.length_x, sgc.origin_y, 0.0)
    points.InsertNextPoint(sgc.origin_x + sgc.length_x,
        sgc.origin_y + sgc.length_y, 0.0)
    points.InsertNextPoint(sgc.origin_x, sgc.origin_y + sgc.length_y, 0.0)
    ug.SetPoints(points)

    quad = vtk.vtkQuad()
    for i in range(4):
        quad.GetPointIds().SetId(i, i)

    ug.InsertNextCell(quad.GetCellType(), quad.GetPointIds())
    add_global_id(cell_dashed_id, sgf, ug)
    return ug

def add_global_id(cell_dashed_id, sgf, unstructured_grid):
    gids = vtk.vtkIdTypeArray()
    gids.SetName("GlobalIds")
    gids.InsertNextTuple1(sgf.get_local_cell_id_from_dashed_cell_id(
        cell_dashed_id))
    unstructured_grid.GetCellData().AddArray(gids)

def get_bucket_file_name(rank, program_data):
    grid_file = program_data["grid_desc"]["read_grid"]
    prefix = get_grid_file_prefix(grid_file)
    return get_bucket_file_name_for_rank(rank, prefix)

class SpartaGridCell:

    def __init__(self, dashed_id, sparta_grid_file, program_data):
        self.__origin_x = 0.0
        self.__origin_y = 0.0
        self.__origin_z = 0.0
        self.__length_x = 1.0
        self.__length_y = 1.0
        self.__length_z = 1.0
        self._calc_origin_and_lengths(dashed_id, sparta_grid_file, program_data)

    @property
    def origin_x(self):
        return self.__origin_x

    @property
    def origin_y(self):
        return self.__origin_y

    @property
    def origin_z(self):
        return self.__origin_z

    @property
    def length_x(self):
        return self.__length_x

    @property
    def length_y(self):
        return self.__length_y

    @property
    def length_z(self):
        return self.__length_z

    def _calc_origin_and_lengths(self, dashed_id,
        sparta_grid_file, program_data):

        top_level_box = program_data['grid_desc']['create_box']
        box = {}
        box["xlo"] = top_level_box['xlo']
        box["ylo"] = top_level_box['ylo']
        box["zlo"] = top_level_box['zlo']
        box["xhi"] = top_level_box['xhi']
        box["yhi"] = top_level_box['yhi']
        box["zhi"] = top_level_box['zhi']
        dimension = program_data['grid_desc']["dimension"]
        cells = SpartaGridFile.get_cells_in_dashed_id(dashed_id)
        cells.reverse()

        for idx, cell in enumerate(cells):
            lvl_dims = sparta_grid_file.get_level_dimensions(idx+1)
            spacing = self._get_level_spacing(box, lvl_dims)
            location = self._get_location_of_cell_id(cell, lvl_dims, dimension)
            box["xlo"] += location[0] * spacing[0]
            box["ylo"] += location[1] * spacing[1]
            box["zlo"] += location[2] * spacing[2]
            box["xhi"] = box["xlo"] + spacing[0]
            box["yhi"] = box["ylo"] + spacing[1]
            box["zhi"] = box["zlo"] + spacing[2]

        self.__origin_x = box["xlo"]
        self.__origin_y = box["ylo"]
        self.__origin_z = box["zlo"]
        self.__length_x = spacing[0]
        self.__length_y = spacing[1]
        self.__length_z = spacing[2]
        if dimension == 2:
            self.__origin_z = 0.0
            self.__length_z = 0.0

    def _get_level_spacing(self, box, level_dimensions):
        spacing_x = (box["xhi"] - box["xlo"])/level_dimensions["x"]
        spacing_y = (box["yhi"] - box["ylo"])/level_dimensions["y"]
        spacing_z = (box["zhi"] - box["zlo"])/level_dimensions["z"]
        return (spacing_x, spacing_y, spacing_z)

    def _get_location_of_cell_id(self, cell_id, level_dimensions, dimension):
        Dx = level_dimensions["x"]
        Dy = level_dimensions["y"]
        Dz = level_dimensions["z"]
        cell_id_start_zero = cell_id - 1
        zloc = 0
        if dimension == 3:
            zloc = cell_id_start_zero//(Dx*Dy)
        yloc = (cell_id_start_zero - zloc*Dx*Dy)//Dx
        xloc = cell_id_start_zero - yloc*Dx - zloc*Dx*Dy
        return (xloc, yloc, zloc)

if __name__ == "__main__":
    main()
