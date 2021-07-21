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
import glob

def main():
    args = parse_command_line()
    check_command_line(args)
    grid_desc = create_grid_description(args)
    time_steps = get_time_steps(args)
    program_data = distribute_program_data(args, grid_desc, time_steps)
    unstructured_grid = create_grid(program_data)
    write_grid(unstructured_grid, program_data)
    if is_rank_zero():
        write_pvd_file(program_data)
    barrier()
    if is_rank_zero():
        print("grid2paraview_cells finished")

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

def get_time_steps(args):
    error_flag = False
    time_steps_dict = None
    if is_rank_zero():
        result_file_list = get_time_steps_file_list(args)
        if result_file_list is not None:
            time_steps_dict = {}
            for f in result_file_list:
                try:
                    fh = open(f, "r")
                except IOError:
                    print("Unable to open SPARTA result file: ", f)
                    error_flag = True
                    break

                for line in fh:
                    s = clean_line(line)
                    if s.lower().replace(" ", "") == "item:timestep":
                        for line in fh:
                            time = int(line)
                            if time in time_steps_dict.keys():
                                time_steps_dict[time].append(f)
                            else:
                                time_steps_dict[time] = [f]
                            break
                    break
                fh.close()
        else:
            error_flag = True

    if error_found_on_rank_zero(error_flag):
        sys.exit(1)

    return time_steps_dict

def clean_line(line):
    line = line.partition('#')[0]
    return line.strip()

def get_time_steps_file_list(args):
    time_steps_file_list = []
    if args.result:
        for f in args.result:
            time_steps_file_list.extend(glob.glob(f))
    elif args.resultfile:
        try:
            rf = open(args.resultfile, "r")
            for name in rf:
                time_steps_file_list.append(name.rstrip())
            rf.close()
        except IOError:
            print("Unable to open SPARTA result file input list file: ",
                args.resultfile)
            time_steps_file_list = None
    return time_steps_file_list

def write_grid(unstructured_grid, program_data):
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(unstructured_grid)
    time_steps = program_data["time_steps"]
    output_prefix = program_data["paraview_output_file"]

    if is_rank_zero():
        print("Writing grid over " + str(len(time_steps.keys())) + " time step(s)")

    if time_steps:
        id_map = create_cell_global_id_to_local_id_map(unstructured_grid)
        for time in sorted(time_steps.keys()):
            read_time_step_data(time_steps[time][0], unstructured_grid, id_map)
            write_grid_to_file(writer, output_prefix, time)
    else:
        write_grid_to_file(writer, output_prefix, 0)

    if is_rank_zero():
        print("Finished writing grid")

def write_grid_to_file(writer, output_prefix, time):
    filepath = os.path.join(output_prefix, output_prefix + '_' + \
        str(get_rank()) + '_' + str(time) + '.vtu')
    writer.SetFileName(filepath)
    writer.Write()

def create_cell_global_id_to_local_id_map(unstructured_grid):
    id_map = {}
    gids = unstructured_grid.GetCellData().GetArray("GlobalIds")
    if gids:
        for i in range(gids.GetNumberOfTuples()):
            id_map[int(gids.GetTuple1(i))] = i
        unstructured_grid.GetCellData().RemoveArray("GlobalIds")
    return id_map

def read_time_step_data(file_name, unstructured_grid, id_hash):
    if is_rank_zero():
        print("Reading Sparta flow file " + file_name)

    fh = open(file_name, "r")
    array_names = get_array_names(fh)

    id_index = 0
    try:
        id_index = array_names.index('id')
    except ValueError:
        print("Error reading SPARTA result file: ", f)
        print("id column not given in file.")
        return

    if not unstructured_grid.GetCellData().GetNumberOfArrays():
        for name in array_names:
            array = vtk.vtkDoubleArray()
            array.SetName(name)
            array.SetNumberOfComponents(1)
            array.SetNumberOfTuples(unstructured_grid.GetNumberOfCells())
            array.FillComponent(0, 0.0)
            unstructured_grid.GetCellData().AddArray(array)

    if unstructured_grid.GetCellData().GetNumberOfArrays() != len(array_names):
        print("Error reading SPARTA result file: ", f)
        print("Expected data columns:  ",
            unstructured_grid.GetCellData().GetNumberOfArrays())
        print("Found data columns:  ", len(array_names))
        return

    arrays = []
    for val in array_names:
        arrays.append(unstructured_grid.GetCellData().GetArray(val))

    count = 0
    for line in fh:
        s = clean_line(line)
        sl = s.split()
        if len(sl) == len(array_names):
            index = int(sl[id_index])
            if index not in id_hash:
                continue
            for idx, val in enumerate(array_names):
                arrays[idx].SetValue(id_hash[index], float(sl[idx]))
        else:
            print("Error reading SPARTA result file: ", f)
            print("Flow data line cannot be processed:  ", line)
            return
        if is_rank_zero() and count % 100000 == 0:
            print("Read " + str(count) + " lines from flow file")
        count += 1

    fh.close()
    if is_rank_zero():
        print("Finished reading Sparta flow file " + file_name)

def write_pvd_file(program_data):
    if is_rank_zero():
        print("Writing pvd file")
    time_steps = program_data["time_steps"]
    output_prefix = program_data["paraview_output_file"]
    fh = open(output_prefix + ".pvd", "w")
    fh.write('<?xml version="1.0"?>\n')
    fh.write('<VTKFile type="Collection" version="0.1"\n')
    fh.write('               byte_order="LittleEndian"\n')
    fh.write('               compressor="vtkZLibDataCompressor">\n')
    fh.write('   <Collection>    \n')

    if time_steps:
        for time in sorted(time_steps.keys()):
            array_names = []
            file_list = time_steps[time]
            if file_list:
                try:
                    afh = open(file_list[0], "r")
                    array_names = get_array_names(afh)
                    afh.close()
                except IOError:
                    print("Unable to open SPARTA result file: ", f)
                    return
            write_pvd_time_step(fh, time, array_names, output_prefix)
    else:
        array_names = []
        time = 0
        write_pvd_time_step(fh, time, array_names, output_prefix)

    fh.write('   </Collection>    \n')
    fh.write('</VTKFile>    \n')
    fh.close()
    if is_rank_zero():
        print("Finished writing pvd file")

def write_pvd_time_step(file_handle, time, array_names, output_prefix):
    file_handle.write('    <DataSet timestep="' + str(time) + '" group="" part="0"   \n')
    filepath = os.path.join(output_prefix, output_prefix + '_'  + str(time) + '.pvtu')
    file_handle.write('             file="' + filepath + '"/>\n')
    write_pvtu_file(array_names, time, output_prefix)

def write_pvtu_file(array_names, time, output_prefix):
    filepath = os.path.join(output_prefix, output_prefix + '_'  +\
        str(time) + '.pvtu')
    fh = open(filepath, "w")
    fh.write('<?xml version="1.0"?>\n')
    fh.write('<VTKFile type="PUnstructuredGrid" version="0.1"\n')
    fh.write('               byte_order="LittleEndian"\n')
    fh.write('               compressor="vtkZLibDataCompressor">\n')
    fh.write('<PUnstructuredGrid GhostLevel="0">\n')
    fh.write('<PCellData>\n')

    for name in array_names:
        fh.write('<PDataArray type="Float64" Name="' + name + '"/>\n')

    fh.write('</PCellData>\n')
    fh.write('<PPoints>\n')
    fh.write('<PDataArray type="Float32" Name="Points" NumberOfComponents="3"/>\n')
    fh.write('</PPoints>\n')

    for rank in range(get_size()):
        fh.write('<Piece Source="' + output_prefix + '_' + str(rank) + '_' + \
            str(time) + '.vtu"/>\n')

    fh.write('</PUnstructuredGrid>\n')
    fh.write('</VTKFile>\n')

def get_array_names(file_handle):
    array_names = []
    for line in file_handle:
        s = clean_line(line)
        if s.lower().replace(" ", "")[:10] == "item:cells":
            for name in s.split()[2:]:
                array_names.append(name)
            break
    return array_names

def distribute_program_data(args, grid_desc, time_steps):
    pd = {}
    if is_rank_zero():
        pd["paraview_output_file"] = args.paraview_output_file
        pd["grid_desc"] = grid_desc
        pd["time_steps"] = time_steps
    else:
        pd["paraview_output_file"] = None
        pd["grid_desc"] = None
        pd["time_steps"] = None

    pd["paraview_output_file"] = get_comm_world().bcast(
        pd["paraview_output_file"], root = 0)
    pd["grid_desc"] = get_comm_world().bcast(pd["grid_desc"], root = 0)
    pd["time_steps"] = get_comm_world().bcast(pd["time_steps"], root = 0)
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
    merge_points = create_merge_points(program_data)
    unstructured_grid = vtk.vtkUnstructuredGrid()
    unstructured_grid.SetPoints(merge_points.GetPoints())
    global_ids = vtk.vtkIdTypeArray()
    global_ids.SetName("GlobalIds")
    count = 1
    if is_rank_zero():
        print("Started sorted grid file read")
        print("Creating grid")
    with open(get_bucket_file_name(get_rank(), program_data), 'r') as f:
        for cell_id in f:
            cell = create_grid_cell(cell_id.strip(), sgf,
                program_data, merge_points)
            unstructured_grid.InsertNextCell(cell.GetCellType(),
                cell.GetPointIds())
            global_ids.InsertNextTuple1(
                sgf.get_local_cell_id_from_dashed_cell_id(cell_id))
            if is_rank_zero() and count % 10000 == 0:
                print("Read " + str(count) + " cells")
            count += 1
    if is_rank_zero():
        print("Finished sorted grid file read")
        print("Finished creating grid")
    unstructured_grid.GetCellData().AddArray(global_ids)
    return unstructured_grid

def create_merge_points(program_data):
    merge_points = vtk.vtkMergePoints()
    points = vtk.vtkPoints()
    top_level_box = program_data['grid_desc']['create_box']
    bounds = []
    bounds.append(top_level_box['xlo'])
    bounds.append(top_level_box['xhi'])
    bounds.append(top_level_box['ylo'])
    bounds.append(top_level_box['yhi'])
    bounds.append(top_level_box['zlo'])
    bounds.append(top_level_box['zhi'])
    merge_points.InitPointInsertion(points, bounds)
    return merge_points

def create_grid_cell(cell_dashed_id, sgf, program_data, merge_points):

    if program_data['grid_desc']['dimension'] == 3:
        return create_hex_cell(cell_dashed_id, sgf, program_data, merge_points)
    else:
        return create_quad_cell(cell_dashed_id, sgf, program_data, merge_points)

def create_hex_cell(cell_dashed_id, sgf, program_data, merge_points):
    sgc = SpartaGridCell(cell_dashed_id, sgf, program_data)
    hex = vtk.vtkHexahedron()

    hex.GetPointIds().SetId(0, merge_points.InsertNextPoint(
        [sgc.origin_x, sgc.origin_y, sgc.origin_z]))
    hex.GetPointIds().SetId(1, merge_points.InsertNextPoint(
        [sgc.origin_x + sgc.length_x, sgc.origin_y, sgc.origin_z]))
    hex.GetPointIds().SetId(2, merge_points.InsertNextPoint(
        [sgc.origin_x + sgc.length_x, sgc.origin_y + sgc.length_y, sgc.origin_z]))
    hex.GetPointIds().SetId(3, merge_points.InsertNextPoint(
        [sgc.origin_x, sgc.origin_y + sgc.length_y, sgc.origin_z]))
    hex.GetPointIds().SetId(4, merge_points.InsertNextPoint(
        [sgc.origin_x, sgc.origin_y, sgc.origin_z + sgc.length_z]))
    hex.GetPointIds().SetId(5, merge_points.InsertNextPoint(
        [sgc.origin_x + sgc.length_x, sgc.origin_y, sgc.origin_z + sgc.length_z]))
    hex.GetPointIds().SetId(6, merge_points.InsertNextPoint(
        [sgc.origin_x + sgc.length_x, sgc.origin_y + sgc.length_y,
            sgc.origin_z + sgc.length_z]))
    hex.GetPointIds().SetId(7, merge_points.InsertNextPoint(
        [sgc.origin_x, sgc.origin_y + sgc.length_y, sgc.origin_z + sgc.length_z]))

    return hex

def create_quad_cell(cell_dashed_id, sgf, program_data, merge_points):
    sgc = SpartaGridCell(cell_dashed_id, sgf, program_data)
    quad = vtk.vtkQuad()

    quad.GetPointIds().SetId(0, merge_points.InsertNextPoint(
        [sgc.origin_x, sgc.origin_y, 0.0]))
    quad.GetPointIds().SetId(1, merge_points.InsertNextPoint(
        [sgc.origin_x + sgc.length_x, sgc.origin_y, 0.0]))
    quad.GetPointIds().SetId(2, merge_points.InsertNextPoint(
        [sgc.origin_x + sgc.length_x, sgc.origin_y + sgc.length_y, 0.0]))
    quad.GetPointIds().SetId(3, merge_points.InsertNextPoint(
        [sgc.origin_x, sgc.origin_y + sgc.length_y, 0.0]))

    return quad

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
