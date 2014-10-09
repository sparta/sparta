
#   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
#   http://sparta.sandia.gov
#   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov,
#   Thomas Otahal, tjotaha@sandia.gov
#   Sandia National Laboratories

#   Copyright (2014) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under 
#   the GNU General Public License.

#   See the README file in the top-level SPARTA directory.

import math
import argparse
import sys
import os
import vtk
import glob
import platform
import multiprocessing as mp

def get_chunk(dim, chunk_size):
  dmod = divmod(dim, chunk_size)
  c = []
  if dmod[0] == 0:
    c.append([1, dmod[1]])
  elif dmod[0] == 1 and dmod[1] == 0:
    c.append([1, chunk_size])
  else:
    for i in range(dmod[0]):
      c.append([i*chunk_size+1, (i+1)*chunk_size])
    if dmod[1] != 0:
      c.append([(dmod[0])*chunk_size+1, (dmod[0])*chunk_size + dmod[1]])
  return c

def find_chunking(chunks, grid_desc, args):
  Cx = grid_desc["create_grid"][1]["Cx"]
  Cy = grid_desc["create_grid"][1]["Cy"]
  xc = get_chunk(Cx, args.xchunk)
  yc = get_chunk(Cy, args.ychunk)

  if grid_desc["dimension"] == 3:
    Cz = grid_desc["create_grid"][1]["Cz"]
    zc = get_chunk(Cz, args.zchunk)

    for k in zc:
      for j in yc:
        for i in xc:
          chunks.append({"x": i, "y": j, "z": k}) 
  else:
    for j in yc:
      for i in xc:
        chunks.append({"x": i, "y": j, "z": [1,1]}) 

def process_grid_chunk(chunk_id, chunk_info, num_chunks, grid_desc, time_steps_dict, output_file):
  xi = grid_desc["create_box"]["xhi"] - grid_desc["create_box"]["xlo"]
  yi = grid_desc["create_box"]["yhi"] - grid_desc["create_box"]["ylo"]
  zi = grid_desc["create_box"]["zhi"] - grid_desc["create_box"]["zlo"]
  px = grid_desc["create_grid"][1]["Cx"]
  py = grid_desc["create_grid"][1]["Cy"]
  pz = grid_desc["create_grid"][1]["Cz"]
  ug = None
  spacing = [xi/float(px), yi/float(py), zi/float(pz)]
  origin = [grid_desc["create_box"]["xlo"] + spacing[0]*(chunk_info["x"][0]-1), \
            grid_desc["create_box"]["ylo"] + spacing[1]*(chunk_info["y"][0]-1), \
            grid_desc["create_box"]["zlo"] + spacing[2]*(chunk_info["z"][0]-1)]
  ndims = [chunk_info["x"][1] - chunk_info["x"][0] + 2, \
           chunk_info["y"][1] - chunk_info["y"][0] + 2, \
           chunk_info["z"][1] - chunk_info["z"][0] + 2]

  if grid_desc["dimension"] == 3:
    ug = create_3d_amr_grids(grid_desc, 1, 0, 0, origin, spacing, ndims, chunk_info)
  elif grid_desc["dimension"] == 2:
    ug = create_2d_amr_grids(grid_desc, 1, 0, 0, origin, spacing, ndims, chunk_info)

  id_hash = {}
  gids = ug.GetCellData().GetArray("GlobalIds")
  for i in range(gids.GetNumberOfTuples()):
    id_hash[int(gids.GetTuple1(i))] = i

  ug.GetCellData().RemoveArray("GlobalIds")

  writer = vtk.vtkXMLUnstructuredGridWriter()
  writer.SetInputData(ug)
  for time in sorted(time_steps_dict.keys()):
    read_time_step_data(time_steps_dict[time], ug, id_hash)
    filepath = os.path.join(output_file, output_file + '_' + str(chunk_id) + '_' + str(time) + '.vtu')
    writer.SetFileName(filepath)
    writer.Write()

  return chunk_id

def create_3d_amr_grids(grid_desc, level, parent_bit_mask, parent_id, origin, spacing, ndims, chunk_info):
  next_level = level+1
  if next_level in grid_desc["create_grid"]:
    if level == 1:
      Dx = grid_desc["create_grid"][1]["Cx"]
      Dy = grid_desc["create_grid"][1]["Cy"]
      Dz = grid_desc["create_grid"][1]["Cz"]
      level_one_bit_mask = int(math.floor(math.log(int(Dx*Dy*Dz),2)) + 1)
      xc = range(chunk_info["x"][0], chunk_info["x"][1]+1)
      yc = range(chunk_info["y"][0], chunk_info["y"][1]+1)
      zc = range(chunk_info["z"][0], chunk_info["z"][1]+1)
    Px = grid_desc["create_grid"][next_level]["Px"]
    Py = grid_desc["create_grid"][next_level]["Py"]
    Pz = grid_desc["create_grid"][next_level]["Pz"]
    Cx = grid_desc["create_grid"][next_level]["Cx"]
    Cy = grid_desc["create_grid"][next_level]["Cy"]
    Cz = grid_desc["create_grid"][next_level]["Cz"]
    xi = spacing[0]
    yi = spacing[1]
    zi = spacing[2]
    r_spacing = [xi/float(Cx), yi/float(Cy), zi/float(Cz)]
    n_spacing = [xi, yi, zi]
    r_ndims = [Cx + 1, Cy + 1, Cz + 1]
    n_ndims = [2, 2, 2]
    num_of_cells = int( (ndims[0]-1)*(ndims[1]-1)*(ndims[2]-1) )
    bit_mask = int(math.floor(math.log(num_of_cells,2)) + 1)
    k_append = vtk.vtkAppendFilter()
    k_append.MergePointsOn()
    cell_index = 1
    for k in range(ndims[2]-1):
      if level == 1:
        zindex = (zc[k]-1)*Dx*Dy
      zl = origin[2] + k*spacing[2]
      y_append = vtk.vtkAppendFilter()
      y_append.MergePointsOn()
      y_ug = None
      for j in range(ndims[1]-1):
        if level == 1:
          yindex = (yc[j]-1)*Dx
        yl = origin[1] + j*spacing[1]
        x_append = vtk.vtkAppendFilter()
        x_append.MergePointsOn()
        x_ug = None
        for i in range(ndims[0]-1):
          xl = origin[0] + i*spacing[0]
 
          if level == 1:
            bit_mask = level_one_bit_mask
            cell_index = xc[i] + yindex + zindex
            refine = xc[i] in Px and yc[j] in Py and zc[k] in Pz
          else:
            refine = i+1 in Px and j+1 in Py and k+1 in Pz

          if refine:
            refined_cell_index = cell_index*(2**parent_bit_mask) + parent_id
            r_ug = create_3d_amr_grids(grid_desc, next_level, bit_mask+parent_bit_mask, refined_cell_index, \
                                       [xl, yl, zl], r_spacing, r_ndims, chunk_info)
          else:
            r_ug = create_3d_amr_grids(grid_desc, next_level, parent_bit_mask, cell_index - 1, \
                                       [xl, yl, zl], n_spacing, n_ndims, chunk_info)
          x_append.AddInputData(r_ug)
          cell_index += 1
        x_append.Update()
        x_ug = x_append.GetOutput()
        y_append.AddInputData(x_ug)
      y_append.Update()
      y_ug = y_append.GetOutput()
      k_append.AddInputData(y_ug)
    k_append.Update()
    return k_append.GetOutput()
  else:
    ug = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    for k in range(ndims[2]):
      k_offset = k*spacing[2]
      for j in range(ndims[1]):
        j_offset = j*spacing[1]
        for i in range(ndims[0]):
          points.InsertNextPoint(origin[0] + i*spacing[0], \
                                 origin[1] + j_offset, \
                                 origin[2] + k_offset)
    ug.SetPoints(points)

    hex = vtk.vtkHexahedron()
    for k in range(ndims[2]-1):
      kl_offset = k*ndims[1]*ndims[0]
      ku_offset = (k+1)*ndims[1]*ndims[0]
      for j in range(ndims[1]-1):
        jl_offset = j*ndims[0]
        ju_offset = (j + 1)*ndims[0]
        ll = jl_offset + kl_offset
        uu = ju_offset + ku_offset
        lu = jl_offset + ku_offset
        ul = ju_offset + kl_offset
        for i in range(ndims[0]-1):
          hex.GetPointIds().SetId(0, i + ll)
          hex.GetPointIds().SetId(1, i + 1 + ll)
          hex.GetPointIds().SetId(2, i + 1 + ul)
          hex.GetPointIds().SetId(3, i + ul)
          hex.GetPointIds().SetId(4, i + lu)
          hex.GetPointIds().SetId(5, i + 1 + lu)
          hex.GetPointIds().SetId(6, i + 1 + uu)
          hex.GetPointIds().SetId(7, i + uu)

          ug.InsertNextCell(hex.GetCellType(), hex.GetPointIds())

    gids = vtk.vtkIdTypeArray()
    gids.SetName("GlobalIds")
    if level == 1:
      Dx = grid_desc["create_grid"][1]["Cx"]
      Dy = grid_desc["create_grid"][1]["Cy"]
      xc = range(chunk_info["x"][0], chunk_info["x"][1]+1)
      yc = range(chunk_info["y"][0], chunk_info["y"][1]+1)
      zc = range(chunk_info["z"][0], chunk_info["z"][1]+1)
      for k in range(ndims[2]-1):
        zindex = (zc[k]-1)*Dx*Dy
        for j in range(ndims[1]-1):
          yindex = (yc[j]-1)*Dx
          for i in range(ndims[0]-1):
            gids.InsertNextTuple1(xc[i]+yindex+zindex)
    else:
      for i in range(ug.GetNumberOfCells()):
        gids.InsertNextTuple1((i+1)*(2**parent_bit_mask) + parent_id)

    ug.GetCellData().AddArray(gids)
    ug.GetCellData().SetActiveGlobalIds("GlobalIds")

  return ug

def create_2d_amr_grids(grid_desc, level, parent_bit_mask, parent_id, origin, spacing, ndims, chunk_info):
  next_level = level+1
  if next_level in grid_desc["create_grid"]:
    if level == 1:
      Dx = grid_desc["create_grid"][1]["Cx"]
      Dy = grid_desc["create_grid"][1]["Cy"]
      level_one_bit_mask = int(math.floor(math.log(int(Dx*Dy),2)) + 1)
      xc = range(chunk_info["x"][0], chunk_info["x"][1]+1)
      yc = range(chunk_info["y"][0], chunk_info["y"][1]+1)
    Px = grid_desc["create_grid"][next_level]["Px"]
    Py = grid_desc["create_grid"][next_level]["Py"]
    Cx = grid_desc["create_grid"][next_level]["Cx"]
    Cy = grid_desc["create_grid"][next_level]["Cy"]
    num_of_cells = int( (ndims[0]-1)*(ndims[1]-1) )
    bit_mask = int(math.floor(math.log(num_of_cells,2)) + 1)
    zl = 0.0
    xi = spacing[0]
    yi = spacing[1]
    zi = 0.0
    r_spacing = [xi/float(Cx), yi/float(Cy), 0.0]
    n_spacing = [xi, yi, zi]
    r_ndims = [Cx + 1, Cy + 1, 1]
    n_ndims = [2, 2, 2]
    y_append = vtk.vtkAppendFilter()
    y_append.MergePointsOn()
    cell_index = 1
    for j in range(ndims[1]-1):
      if level == 1:
        yindex = (yc[j]-1)*Dx
      yl = origin[1] + j*spacing[1]
      x_append = vtk.vtkAppendFilter()
      x_append.MergePointsOn()
      x_ug = None
      for i in range(ndims[0]-1):
        xl = origin[0] + i*spacing[0]

        if level == 1:
          bit_mask = level_one_bit_mask
          cell_index = xc[i] + yindex
          refine = xc[i] in Px and yc[j] in Py
        else:
          refine = i+1 in Px and j+1 in Py

        if refine:
          r_ug = create_2d_amr_grids(grid_desc, next_level, bit_mask+parent_bit_mask, cell_index, \
                                     [xl, yl, zl], r_spacing, r_ndims, chunk_info)
        else:
          r_ug = create_2d_amr_grids(grid_desc, next_level, parent_bit_mask, cell_index - 1, \
                                     [xl, yl, zl], n_spacing, n_ndims, chunk_info)
        x_append.AddInputData(r_ug)
        cell_index += 1
      x_append.Update()
      x_ug = x_append.GetOutput()
      y_append.AddInputData(x_ug)
    y_append.Update()
    return y_append.GetOutput()
  else:
    ug = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    for j in range(ndims[1]):
      j_offset = j*spacing[1]
      for i in range(ndims[0]):
        points.InsertNextPoint(origin[0] + i*spacing[0], \
                               origin[1] + j_offset, \
                               0.0)
    ug.SetPoints(points)

    quad = vtk.vtkQuad()
    for j in range(ndims[1]-1):
      jl_offset = j*ndims[0]
      ju_offset = (j+1)*ndims[0]
      for i in range(ndims[0]-1):
        quad.GetPointIds().SetId(0, i + jl_offset)
        quad.GetPointIds().SetId(1, i + 1 + jl_offset)
        quad.GetPointIds().SetId(2, i + 1 + ju_offset)
        quad.GetPointIds().SetId(3, i + ju_offset)
        ug.InsertNextCell(quad.GetCellType(), quad.GetPointIds())

    gids = vtk.vtkIdTypeArray()
    gids.SetName("GlobalIds")
    if level == 1:
      Dx = grid_desc["create_grid"][1]["Cx"]
      xc = range(chunk_info["x"][0], chunk_info["x"][1]+1)
      yc = range(chunk_info["y"][0], chunk_info["y"][1]+1)
      for j in range(ndims[1]-1):
        yindex = (yc[j]-1)*Dx
        for i in range(ndims[0]-1):
          gids.InsertNextTuple1(xc[i]+yindex)
    else: 
      for i in range(ug.GetNumberOfCells()):
        gids.InsertNextTuple1((i+1)*(2**parent_bit_mask) + parent_id)

    ug.GetCellData().AddArray(gids)
    ug.GetCellData().SetActiveGlobalIds("GlobalIds")

  return ug
  
def clean_line(line):
  line = line.partition('#')[0]
  return line.strip()

def create_parent_ranges(item, level, index, parent, grid_desc):
  if item == "*":
    grid_desc["create_grid"][level][index] = range(1, parent + 1)
  elif item[0] == "*": 
    rb = int(item.split('*')[1])
    grid_desc["create_grid"][level][index] = range(1, rb + 1)
  elif item[-1:] == "*": 
    lb = int(item.split('*')[0])
    grid_desc["create_grid"][level][index] = range(lb, parent + 1)
  elif len(item.split('*')) == 1:
    grid_desc["create_grid"][level][index] = [int(item)]
  elif len(item.split('*')) == 2:
    b = item.split('*')
    grid_desc["create_grid"][level][index] = range(int(b[0]), int(b[1]) + 1)

def read_grid_levels(gl_array, grid_desc, level):
  if len(gl_array) % 8 or \
     gl_array[0].lower() != "level" or \
     int(gl_array[1]) != level: 
    print "Error reading SPARTA grid description file"
    print "create_grid specification is invalid: ", ' '.join(gl_array)
    sys.exit(1)

  grid_desc["create_grid"][level] = {}
  Px = grid_desc["create_grid"][level - 1]["Cx"]
  Py = grid_desc["create_grid"][level - 1]["Cy"]
  Pz = grid_desc["create_grid"][level - 1]["Cz"]

  create_parent_ranges(gl_array[2], level, "Px", Px, grid_desc)
  create_parent_ranges(gl_array[3], level, "Py", Py, grid_desc)
  create_parent_ranges(gl_array[4], level, "Pz", Pz, grid_desc)

  grid_desc["create_grid"][level]["Cx"] = int(gl_array[5])
  grid_desc["create_grid"][level]["Cy"] = int(gl_array[6])
  grid_desc["create_grid"][level]["Cz"] = int(gl_array[7])

  if len(gl_array) > 8:
    read_grid_levels(gl_array[8:], grid_desc, level + 1)

def read_grid_description_file(sif, grid_desc):
  for line in sif:
    s = clean_line(line)
    if s.lower()[:9] == "dimension" and len(s.split()) == 2:
      dimension = int(s.split()[1])
      if dimension != 2 and dimension != 3:
        print "Error reading SPARTA grid description file"
        print "dimension must be either 2 or 3: ", dimension
        sys.exit(1)
      else:
        grid_desc["dimension"] = dimension
    elif s.lower()[:10] == "create_box" and len(s.split()) == 7:
      grid_desc["create_box"] = {}
      if s.split()[1] < s.split()[2] and \
         s.split()[3] < s.split()[4] and \
         s.split()[5] < s.split()[6]:
         grid_desc["create_box"]["xlo"] = float(s.split()[1])
         grid_desc["create_box"]["xhi"] = float(s.split()[2])
         grid_desc["create_box"]["ylo"] = float(s.split()[3])
         grid_desc["create_box"]["yhi"] = float(s.split()[4])
         grid_desc["create_box"]["zlo"] = float(s.split()[5])
         grid_desc["create_box"]["zhi"] = float(s.split()[6])
      else:
        print "Error reading SPARTA grid description file"
        print "create_box specification is invalid: ", s
        sys.exit(1)
    elif s.lower()[:11] == "create_grid" and len(s.split()) > 3:
      grid_desc["create_grid"] = {}
      if int(s.split()[1]) > 0 and \
         int(s.split()[2]) > 0 and \
         int(s.split()[3]) > 0:
        grid_desc["create_grid"][1] = {}
        grid_desc["create_grid"][1]["Cx"] = int(s.split()[1])
        grid_desc["create_grid"][1]["Cy"] = int(s.split()[2])
        grid_desc["create_grid"][1]["Cz"] = int(s.split()[3])
        if len(s.split()) > 4:
          read_grid_levels(s.split()[4:], grid_desc, 2)
      else:
        print "Error reading SPARTA grid description file"
        print "create_grid specification is invalid: ", s
    elif len(s):
      print "Error reading SPARTA grid description file"
      print "File contains unrecognized keyword: ", s
      sys.exit(1)

def read_time_steps(result_file_list, time_steps_dict):
  for f in result_file_list:
    try:
      fh = open(f, "r")
    except IOError:
      print "Unable to open SPARTA result file: ", f
      sys.exit(1)

    for line in fh:
      s = clean_line(line)
      if s.lower().replace(" ", "") == "item:timestep":
        time = int(fh.next())
        if time in time_steps_dict.keys():
          time_steps_dict[time].append(f)
        else:
          time_steps_dict[time] = [f]
        break

    fh.close()

def read_time_step_data(time_step_file_list, ug, id_hash):
  for f in time_step_file_list:
    try:
      fh = open(f, "r")
    except IOError:
      print "Unable to open SPARTA result file: ", f
      return

    array_names = []
    for line in fh:
      s = clean_line(line)
      if s.lower().replace(" ", "")[:10] == "item:cells":
        for name in s.split()[2:]: 
          array_names.append(name)
        break

    id_index = 0
    try:
      id_index = array_names.index('id')
    except ValueError:
      print "Error reading SPARTA result file: ", f
      print "id column not given in file."
      return

    if not ug.GetCellData().GetNumberOfArrays():
      for name in array_names:
        array = vtk.vtkDoubleArray()
        array.SetName(name)
        array.SetNumberOfComponents(1)
        array.SetNumberOfTuples(ug.GetNumberOfCells())
        array.FillComponent(0, 0.0)
        ug.GetCellData().AddArray(array)

    if ug.GetCellData().GetNumberOfArrays() != len(array_names):
      print "Error reading SPARTA result file: ", f
      print "Expected data columns:  ", ug.GetCellData().GetNumberOfArrays()
      print "Found data columns:  ", len(array_names)
      return
   
    arrays = []
    for val in array_names:
      arrays.append(ug.GetCellData().GetArray(val))

    cells_read = 0
    for line in fh:
      s = clean_line(line)
      sl = s.split()
      if len(sl) == len(array_names):
        index = int(sl[id_index])
        if index not in id_hash:
          continue
        for idx, val in enumerate(array_names):
          arrays[idx].SetValue(id_hash[index], float(sl[idx]))
        cells_read += 1
        if cells_read == ug.GetNumberOfCells():
          break
      else:
        print "Error reading SPARTA result file: ", f
        print "Flow data line cannot be processed:  ", line
        return

    fh.close()

def write_pvd_file(time_steps, file_name, chunking):
  fh = open(file_name + ".pvd", "w")
  fh.write('<?xml version="1.0"?>\n')
  fh.write('<VTKFile type="Collection" version="0.1"\n')
  fh.write('               byte_order="LittleEndian"\n')
  fh.write('               compressor="vtkZLibDataCompressor">\n')
  fh.write('   <Collection>    \n')
  
  for time in time_steps:
    for chunk in range(len(chunking)):
      filepath = os.path.join(file_name, file_name + '_' + str(chunk) + '_' + str(time) + '.vtu') 
      fh.write('    <DataSet timestep="' + str(time) + '" group="" part="' + str(chunk) + '"   \n')
      fh.write('             file="' + filepath + '"/>\n')

  fh.write('   </Collection>    \n')
  fh.write('</VTKFile>    \n')
  fh.close()

def report_chunk_complete(chunk_id):
  report_chunk_complete.count += 1
  rem = report_chunk_complete.num_chunks - report_chunk_complete.count
  print "Completed grid chunk number: ", chunk_id, ",", rem, "chunk(s) remaining" 
report_chunk_complete.count = 0
report_chunk_complete.num_chunks = 0

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("sparta_grid_description_file", help="SPARTA grid description input file name")
  parser.add_argument("paraview_output_file", help="ParaView output file name")
  parser.add_argument('-r', '--result', help="Optional list of SPARTA dump result files", nargs='+')
  parser.add_argument('-xc', '--xchunk', \
                      help="Optional x grid chunk size (positive integer; default 100)", \
                      default=100, type=int )
  parser.add_argument('-yc', '--ychunk', \
                      help="Optional y grid chunk size (positive integer; default 100)", \
                      default=100, type=int )
  parser.add_argument('-zc', '--zchunk', \
                      help="Optional z grid chunk size (positive integer; default 100)", \
                      default=100, type=int )
  args = parser.parse_args()

  try:
    gdf = open(args.sparta_grid_description_file, "r")
  except IOError:
    print "Unable to open SPARTA surf input file: ", args.sparta_grid_description_file
    sys.exit(1)

  if os.path.isfile(args.paraview_output_file + '.pvd'):
    print "ParaView output file exists: ", args.paraview_output_file + '.pvd'
    sys.exit(1)

  if os.path.isdir(args.paraview_output_file):
    print "ParaView output directory exists: ", args.paraview_output_file
    sys.exit(1)
  else:
    os.mkdir(args.paraview_output_file)
 
  if args.xchunk < 1:
    print "Invalid xchunk size given: ", args.xchunk
    sys.exit(1)

  if args.ychunk < 1:
    print "Invalid ychunk size given: ", args.ychunk
    sys.exit(1)

  if args.zchunk < 1:
    print "Invalid zchunk size given: ", args.zchunk
    sys.exit(1)

  grid_desc = {}
  read_grid_description_file(gdf, grid_desc)
  gdf.close()

  time_steps_dict = {}
  time_steps_file_list = []
  if args.result:
    for f in args.result:
      time_steps_file_list.extend(glob.glob(f))
  else:
    time_steps_dict[0] = []

  read_time_steps(time_steps_file_list, time_steps_dict)

  chunking = []
  find_chunking(chunking, grid_desc, args)

  sys.stdin = open(os.devnull)

  print "Processing ", len(chunking), " grid chunk(s)."
  report_chunk_complete.num_chunks = len(chunking)

  import platform
  if platform.system() == 'Linux' or platform.system() == 'Darwin':
    import multiprocessing as mp

    pool = mp.Pool()
    for idx, chunk in enumerate(chunking):
      pool.apply_async(process_grid_chunk, \
                       args=(idx, chunk, len(chunking), grid_desc, time_steps_dict, args.paraview_output_file, ), \
                       callback = report_chunk_complete)
    pool.close()
    pool.join()
  else:
    for idx, chunk in enumerate(chunking):
      res = process_grid_chunk(idx, chunk, len(chunking), grid_desc, time_steps_dict, args.paraview_output_file)
      report_chunk_complete(res)

  write_pvd_file(sorted(time_steps_dict.keys()), args.paraview_output_file, chunking)

  print "Done."

