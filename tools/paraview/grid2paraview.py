
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
import platform
import time
import glob
from datetime import timedelta

def open_grid_file(filename):
  gf = None
  try:
    if filename.lower().endswith('.gz'):
      import gzip
      gf = gzip.open(filename, "r")
    else:
      gf = open(filename, "r")
  except IOError:
    print "Unable to open SPARTA grid file: ", filename
    sys.exit(1)
  return gf

def create_grid_from_grid_file(grid_desc):
  gf = open_grid_file(grid_desc["read_grid"])

  for line in gf:
    s = clean_line(line)
    if len(s.split()) == 5:
      if int(s.split()[0]) == 1 and \
         int(s.split()[1]) == 0:
        grid_desc["create_grid"] = {}
        grid_desc["create_grid"][1] = {}
        grid_desc["create_grid"][1]["Cx"] = int(s.split()[2])
        grid_desc["create_grid"][1]["Cy"] = int(s.split()[3])
        grid_desc["create_grid"][1]["Cz"] = int(s.split()[4])
        if grid_desc["create_grid"][1]["Cz"] == 1:
          grid_desc["dimension"] = 2
        else:
          grid_desc["dimension"] = 3
      else:
        print "Error reading SPARTA grid file"
        print "top level grid specification is invalid: ", s
        sys.exit(1)
      gf.close()
      return

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

def create_and_write_grid_chunk(chunk_id, chunk_info, num_chunks, \
                                grid_desc, time_steps_dict, output_file):
  ug = process_grid_chunk(chunk_id, chunk_info, num_chunks, \
    grid_desc, time_steps_dict, output_file)
  write_grid_chunk(ug, chunk_id, num_chunks, \
    grid_desc, time_steps_dict, output_file)

def process_grid_chunk(chunk_id, chunk_info, num_chunks, \
                       grid_desc, time_steps_dict, output_file):
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

  if "read_grid" in grid_desc:
    read_grid_file(grid_desc, chunk_info)

  if grid_desc["dimension"] == 3:
    ug = create_3d_amr_grids(grid_desc, 1, 0, 0, origin, \
      spacing, ndims, chunk_info, "")
  elif grid_desc["dimension"] == 2:
    ug = create_2d_amr_grids(grid_desc, 1, 0, 0, origin, \
      spacing, ndims, chunk_info, "")

  return ug

def create_cell_global_id_to_local_id_map(ug):
  id_map = {}
  gids = ug.GetCellData().GetArray("GlobalIds")
  if gids:
    for i in range(gids.GetNumberOfTuples()):
      id_map[int(gids.GetTuple1(i))] = i
    ug.GetCellData().RemoveArray("GlobalIds")
  return id_map

def write_grid_chunk(ug, chunk_id, num_chunks, grid_desc, time_steps_dict, \
  output_file):
  if not ug:
    return

  id_map = create_cell_global_id_to_local_id_map(ug)

  writer = vtk.vtkXMLUnstructuredGridWriter()
  writer.SetInputData(ug)
  for idx, time in enumerate(sorted(time_steps_dict.keys())):
    pt = ParallelTimer()
    read_time_step_data(time_steps_dict[time], ug, id_map)
    filepath = os.path.join(output_file, output_file + '_' + str(chunk_id) +
      '_' + str(time) + '.vtu')
    writer.SetFileName(filepath)
    writer.Write()
    pt.report_collective_time("grid file write time %s (wall clock time) for step "\
      + str(idx))

  cr = ChunkReport()
  cr.reportChunkComplete(chunk_id)

def create_3d_amr_grids(grid_desc, level, parent_bit_mask, parent_id, \
                        origin, spacing, ndims, chunk_info, dashed_id):
  if level_contains_refined_cells(level, grid_desc, dashed_id):
    if level == 1:
      Dx = grid_desc["create_grid"][1]["Cx"]
      Dy = grid_desc["create_grid"][1]["Cy"]
      Dz = grid_desc["create_grid"][1]["Cz"]
      level_one_bit_mask = int(math.floor(math.log(int(Dx*Dy*Dz),2)) + 1)
      xc = range(chunk_info["x"][0], chunk_info["x"][1]+1)
      yc = range(chunk_info["y"][0], chunk_info["y"][1]+1)
      zc = range(chunk_info["z"][0], chunk_info["z"][1]+1)
    cell_info = {'Cx' : 1, 'Cy' : 1, 'Cz' : 1}
    get_cell_size(level, grid_desc, cell_info)
    Cx = cell_info['Cx']
    Cy = cell_info['Cy']
    Cz = cell_info['Cz']
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
            refine = is_3d_cell_refined(level, xc[i], yc[j], zc[k], cell_index, \
                                        grid_desc, dashed_id, xi, yi, zi, r_spacing, r_ndims)
          else:
            refine = is_3d_cell_refined(level, i+1, j+1, k+1, cell_index, grid_desc, \
                                        dashed_id, xi, yi, zi, r_spacing, r_ndims)

          refined_cell_index = cell_index*(2**parent_bit_mask) + parent_id

          if refine:
            if not dashed_id:
              next_dashed_id = str(cell_index)
            else:
              next_dashed_id = dashed_id + "-" + str(cell_index)
            r_ug = create_3d_amr_grids(grid_desc, level+1, bit_mask+parent_bit_mask, refined_cell_index, \
                                       [xl, yl, zl], r_spacing, r_ndims, chunk_info, next_dashed_id)
          else:
            r_ug = build_3d_grid(-1, refined_cell_index, [xl, yl, zl], n_spacing, n_ndims, chunk_info, grid_desc)
          if r_ug:
            x_append.AddInputData(r_ug)
          cell_index += 1
        if x_append.GetInputList().GetNumberOfItems():
          x_append.Update()
          x_ug = x_append.GetOutput()
          y_append.AddInputData(x_ug)
      if y_append.GetInputList().GetNumberOfItems():
        y_append.Update()
        y_ug = y_append.GetOutput()
        k_append.AddInputData(y_ug)
    if k_append.GetInputList().GetNumberOfItems():
      k_append.Update()
      return k_append.GetOutput()
    else:
      return None
  else:
    return build_3d_grid(parent_bit_mask, parent_id, origin, spacing, ndims, chunk_info, grid_desc)

def read_grid_file(grid_desc, chunk_info):
  gf = open_grid_file(grid_desc["read_grid"])
  Dx = grid_desc["create_grid"][1]["Cx"]
  Dy = grid_desc["create_grid"][1]["Cy"]
  grid_desc["parent_grid"] = {}

  for line in gf:
    s = clean_line(line)
    if len(s.split()) == 5:
      id = s.split()[1].split('-')
      if len(id) and int(id[0]) != 0: 
        index = int(id[0]) - 1
        zloc = 0
        if grid_desc["dimension"] == 3:
          zloc = math.floor(index/(Dx*Dy))
        yloc = math.floor((index - zloc*Dx*Dy)/Dx)
        xloc = index - yloc*Dx - zloc*Dx*Dy
        xloc += 1
        yloc += 1
        zloc += 1
        if xloc >= chunk_info["x"][0] and xloc <= chunk_info["x"][1] and \
           yloc >= chunk_info["y"][0] and yloc <= chunk_info["y"][1] and \
           zloc >= chunk_info["z"][0] and zloc <= chunk_info["z"][1]:
          cld = grid_desc["parent_grid"]
          for pid in id:
            if int(pid) in cld:
              cld = cld[int(pid)]['np']
            else:
              cld[int(pid)] = {'px':int(s.split()[2]), 'py':int(s.split()[3]), 'pz':int(s.split()[4]), 'np':{}}
  gf.close()

def level_contains_refined_cells(level, grid_desc, dashed_id):
  if "parent_grid" in grid_desc:
    if level == 1: 
      return bool(grid_desc["parent_grid"])
    else:
      s = dashed_id.split('-')
      d = None 
      for id in s:
        if not d:
          if int(id) in  grid_desc["parent_grid"]:
            d = grid_desc["parent_grid"][int(id)]['np']
          else:
            return False
        else:
          if int(id) in d:
            d = d[int(id)]['np']
          else:
            return False
      return bool(d)
  else:
    return level+1 in grid_desc["create_grid"]

def is_3d_cell_refined(level, i, j, k, cell_index, grid_desc, dashed_id, xi, yi, zi, r_spacing, r_ndims):
  if "parent_grid" not in grid_desc:
    Px = grid_desc["create_grid"][level+1]["Px"]
    Py = grid_desc["create_grid"][level+1]["Py"]
    Pz = grid_desc["create_grid"][level+1]["Pz"]
    return (i in Px and j in Py and k in Pz)
  else:
    if level == 1:
      if cell_index in grid_desc["parent_grid"]:
        Cx = grid_desc["parent_grid"][cell_index]['px']
        Cy = grid_desc["parent_grid"][cell_index]['py']
        Cz = grid_desc["parent_grid"][cell_index]['pz']
        r_spacing[0] = xi/float(Cx)
        r_spacing[1] = yi/float(Cy)
        r_spacing[2] = zi/float(Cz)
        r_ndims[0] = Cx + 1
        r_ndims[1] = Cy + 1
        r_ndims[2] = Cz + 1
        return True
      else:
        return False
    else:
      s = dashed_id.split('-')
      d = None
      lc = 2
      for id in s:
        if not d:
          d = grid_desc["parent_grid"][int(id)]['np']
        else:
          d = d[int(id)]['np']
        if cell_index in d and lc == level:
          Cx = d[cell_index]['px']
          Cy = d[cell_index]['py']
          Cz = d[cell_index]['pz']
          r_spacing[0] = xi/float(Cx)
          r_spacing[1] = yi/float(Cy)
          r_spacing[2] = zi/float(Cz)
          r_ndims[0] = Cx + 1
          r_ndims[1] = Cy + 1
          r_ndims[2] = Cz + 1
          return True
        lc += 1
      return False

def is_2d_cell_refined(level, i, j, cell_index, grid_desc, dashed_id, xi, yi, r_spacing, r_ndims):
  if "parent_grid" not in grid_desc:
    Px = grid_desc["create_grid"][level+1]["Px"]
    Py = grid_desc["create_grid"][level+1]["Py"]
    return (i in Px and j in Py)
  else:
    if level == 1: 
      if cell_index in grid_desc["parent_grid"]:
        Cx = grid_desc["parent_grid"][cell_index]['px']
        Cy = grid_desc["parent_grid"][cell_index]['py']
        r_spacing[0] = xi/float(Cx)
        r_spacing[1] = yi/float(Cy)
        r_ndims[0] = Cx + 1
        r_ndims[1] = Cy + 1
        return True
      else:
        return False
    else:
      s = dashed_id.split('-')
      d = None 
      lc = 2
      for id in s:
        if not d:
          d = grid_desc["parent_grid"][int(id)]['np']
        else:
          d = d[int(id)]['np']
        if cell_index in d and lc == level:
          Cx = d[cell_index]['px']
          Cy = d[cell_index]['py']
          r_spacing[0] = xi/float(Cx)
          r_spacing[1] = yi/float(Cy)
          r_ndims[0] = Cx + 1
          r_ndims[1] = Cy + 1
          return True
        lc += 1
      return False
     
def get_cell_size(level, grid_desc, cell_info):
  if "parent_grid" not in grid_desc:
    cell_info['Cx'] = grid_desc["create_grid"][level+1]["Cx"]
    cell_info['Cy'] = grid_desc["create_grid"][level+1]["Cy"]
    if grid_desc["dimension"] == 3:
      cell_info['Cz'] = grid_desc["create_grid"][level+1]["Cz"]

def find_2d_intersected_cells(intersecting_planes, parent_bit_mask, parent_id, origin, \
                              spacing, ndims, chunk_info, grid_desc):
  append = vtk.vtkAppendFilter()
  append.MergePointsOn()
  a1 = spacing[0]/2.0
  a2 = spacing[1]/2.0
  index = 1
  for j in range(ndims[1]-1):
    p2 = a2 + origin[1] + j*spacing[1]
    for i in range(ndims[0]-1):
      p1 = a1 + origin[0] + i*spacing[0]
      for plane in intersecting_planes:
        n1 = plane["nx"]
        n2 = plane["ny"]
        n3 = plane["nz"]
        p01 = plane["px"]
        p02 = plane["py"]
        p03 = plane["pz"]
 
        d = math.fabs(n1*(p1-p01) + n2*(p2-p02))
        rhs = a1*n1 + a2*n2
        if d <= math.fabs(rhs):
          ug = vtk.vtkUnstructuredGrid()

          points = vtk.vtkPoints()
          for pj in range(2):
            j_offset = pj*spacing[1]
            for pi in range(2):
              points.InsertNextPoint(p1 - a1 + pi*spacing[0], \
                                     p2 - a2 + j_offset, \
                                     0.0)

          ug.SetPoints(points)

          quad = vtk.vtkQuad()

          quad.GetPointIds().SetId(0, 0)
          quad.GetPointIds().SetId(1, 1)
          quad.GetPointIds().SetId(2, 3)
          quad.GetPointIds().SetId(3, 2)
       
          ug.InsertNextCell(quad.GetCellType(), quad.GetPointIds())

          gids = vtk.vtkIdTypeArray()
          gids.SetName("GlobalIds")
          if parent_id == 0:
            Dx = grid_desc["create_grid"][1]["Cx"]
            xc = range(chunk_info["x"][0], chunk_info["x"][1]+1)
            yc = range(chunk_info["y"][0], chunk_info["y"][1]+1)
            yindex = (yc[j]-1)*Dx
            gids.InsertNextTuple1(xc[i]+yindex)
          else:
            if parent_bit_mask == -1:
              gids.InsertNextTuple1(parent_id)
            else:
              gids.InsertNextTuple1(index*(2**parent_bit_mask) + parent_id)

          ug.GetCellData().AddArray(gids)
          append.AddInputData(ug)
          break
      index += 1
  if append.GetInputList().GetNumberOfItems():
    append.Update()
    return append.GetOutput()
  else:
    return None

def find_3d_intersected_cells(intersecting_planes, parent_bit_mask, parent_id, origin, \
                              spacing, ndims, chunk_info, grid_desc):
  append = vtk.vtkAppendFilter()
  append.MergePointsOn()
  a1 = spacing[0]/2.0
  a2 = spacing[1]/2.0
  a3 = spacing[2]/2.0
  index = 1
  for k in range(ndims[2]-1):
    p3 = a3 + origin[2] + k*spacing[2]
    for j in range(ndims[1]-1):
      p2 = a2 + origin[1] + j*spacing[1]
      for i in range(ndims[0]-1):
        p1 = a1 + origin[0] + i*spacing[0]
        for plane in intersecting_planes:
          n1 = plane["nx"]
          n2 = plane["ny"]
          n3 = plane["nz"]
          p01 = plane["px"]
          p02 = plane["py"]
          p03 = plane["pz"]
          d = math.fabs(n1*(p1-p01) + n2*(p2-p02) + n3*(p3-p03))
          rhs = a1*n1 + a2*n2 + a3*n3
          if d <= math.fabs(rhs):
            ug = vtk.vtkUnstructuredGrid()
            points = vtk.vtkPoints()
            for pk in range(2):
              k_offset = pk*spacing[2]
              for pj in range(2):
                j_offset = pj*spacing[1]
                for pi in range(2):
                  points.InsertNextPoint(p1 - a1 + pi*spacing[0], \
                                         p2 - a2 + j_offset, \
                                         p3 - a3 + k_offset)

            ug.SetPoints(points)

            hex = vtk.vtkHexahedron()

            hex.GetPointIds().SetId(0, 0)
            hex.GetPointIds().SetId(1, 1)
            hex.GetPointIds().SetId(2, 3)
            hex.GetPointIds().SetId(3, 2)
            hex.GetPointIds().SetId(4, 4)
            hex.GetPointIds().SetId(5, 5)
            hex.GetPointIds().SetId(6, 7)
            hex.GetPointIds().SetId(7, 6)

            ug.InsertNextCell(hex.GetCellType(), hex.GetPointIds())

            gids = vtk.vtkIdTypeArray()
            gids.SetName("GlobalIds")
            if parent_id == 0:
              Dx = grid_desc["create_grid"][1]["Cx"]
              Dy = grid_desc["create_grid"][1]["Cy"]
              xc = range(chunk_info["x"][0], chunk_info["x"][1]+1)
              yc = range(chunk_info["y"][0], chunk_info["y"][1]+1)
              zc = range(chunk_info["z"][0], chunk_info["z"][1]+1)
              zindex = (zc[k]-1)*Dx*Dy
              yindex = (yc[j]-1)*Dx
              gids.InsertNextTuple1(xc[i]+yindex+zindex)
            else:
              if parent_bit_mask == -1:
                gids.InsertNextTuple1(parent_id)
              else:
                gids.InsertNextTuple1(index*(2**parent_bit_mask) + parent_id)

            ug.GetCellData().AddArray(gids)
            append.AddInputData(ug)
            break
        index += 1
  if append.GetInputList().GetNumberOfItems():
    append.Update()
    return append.GetOutput()
  else:
    return None

def cells_on_slice_planes(parent_bit_mask, parent_id, origin, spacing, ndims, chunk_info, grid_desc):
  intersecting_planes = []
  for plane in grid_desc["slice"]:
    intersecting_planes.append(plane)

  if grid_desc["dimension"] == 3:
    return find_3d_intersected_cells(intersecting_planes, parent_bit_mask, parent_id, origin, \
                                     spacing, ndims, chunk_info, grid_desc)
  else:
    return find_2d_intersected_cells(intersecting_planes, parent_bit_mask, parent_id, origin, \
                                     spacing, ndims, chunk_info, grid_desc)

def build_3d_grid(parent_bit_mask, parent_id, origin, spacing, ndims, chunk_info, grid_desc):
  if "slice" in grid_desc:
    return cells_on_slice_planes(parent_bit_mask, parent_id, origin, \
                                 spacing, ndims, chunk_info, grid_desc)

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
  if parent_id == 0:
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
      if parent_bit_mask == -1:
        gids.InsertNextTuple1(parent_id)
      else:
        gids.InsertNextTuple1((i+1)*(2**parent_bit_mask) + parent_id)

  ug.GetCellData().AddArray(gids)
  return ug

def build_2d_grid(parent_bit_mask, parent_id, origin, spacing, ndims, chunk_info, grid_desc):
  if "slice" in grid_desc:
    return cells_on_slice_planes(parent_bit_mask, parent_id, origin, \
                                 spacing, ndims, chunk_info, grid_desc)

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
  if parent_id == 0:
    Dx = grid_desc["create_grid"][1]["Cx"]
    xc = range(chunk_info["x"][0], chunk_info["x"][1]+1)
    yc = range(chunk_info["y"][0], chunk_info["y"][1]+1)
    for j in range(ndims[1]-1):
      yindex = (yc[j]-1)*Dx
      for i in range(ndims[0]-1):
        gids.InsertNextTuple1(xc[i]+yindex)
  else:
    for i in range(ug.GetNumberOfCells()):
      if parent_bit_mask == -1:
        gids.InsertNextTuple1(parent_id)
      else:
        gids.InsertNextTuple1((i+1)*(2**parent_bit_mask) + parent_id)

  ug.GetCellData().AddArray(gids)
  return ug

def create_2d_amr_grids(grid_desc, level, parent_bit_mask, parent_id, \
                        origin, spacing, ndims, chunk_info, dashed_id):
  if level_contains_refined_cells(level, grid_desc, dashed_id):
    if level == 1:
      Dx = grid_desc["create_grid"][1]["Cx"]
      Dy = grid_desc["create_grid"][1]["Cy"]
      level_one_bit_mask = int(math.floor(math.log(int(Dx*Dy),2)) + 1)
      xc = range(chunk_info["x"][0], chunk_info["x"][1]+1)
      yc = range(chunk_info["y"][0], chunk_info["y"][1]+1)
    cell_info = {'Cx' : 1, 'Cy' : 1}
    get_cell_size(level, grid_desc, cell_info)
    Cx = cell_info['Cx']
    Cy = cell_info['Cy']
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
          refine = is_2d_cell_refined(level, xc[i], yc[j], cell_index, \
                                      grid_desc, dashed_id, xi, yi, r_spacing, r_ndims)
        else:
          refine = is_2d_cell_refined(level, i+1, j+1, cell_index, grid_desc, \
                                      dashed_id, xi, yi, r_spacing, r_ndims)

        refined_cell_index = cell_index*(2**parent_bit_mask) + parent_id
        
        if refine:
          if not dashed_id:
            next_dashed_id = str(cell_index)
          else:
            next_dashed_id = dashed_id + "-" + str(cell_index)
          r_ug = create_2d_amr_grids(grid_desc, level+1, bit_mask+parent_bit_mask, refined_cell_index, \
                                     [xl, yl, zl], r_spacing, r_ndims, chunk_info, next_dashed_id)
        else:
          r_ug = build_2d_grid(-1, refined_cell_index, [xl, yl, zl], n_spacing, n_ndims, chunk_info, grid_desc)
        x_append.AddInputData(r_ug)
        cell_index += 1
      if x_append.GetInputList().GetNumberOfItems():
        x_append.Update()
        x_ug = x_append.GetOutput()
        y_append.AddInputData(x_ug)
    if y_append.GetInputList().GetNumberOfItems():
      y_append.Update()
      return y_append.GetOutput()
    else:
      return None
  else:
    return build_2d_grid(parent_bit_mask, parent_id, origin, spacing, ndims, chunk_info, grid_desc)
  
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
    elif s.lower()[:5] == "slice" and len(s.split()) == 7:
      p = {}
      p["nx"] = float(s.split()[1])
      p["ny"] = float(s.split()[2])
      p["nz"] = float(s.split()[3])
      norm = math.sqrt(math.pow(p["nx"],2) + math.pow(p["ny"],2) + math.pow(p["nz"],2))
      p["nx"] = p["nx"]/norm
      p["ny"] = p["ny"]/norm
      p["nz"] = p["nz"]/norm
      p["px"] = float(s.split()[4])
      p["py"] = float(s.split()[5])
      p["pz"] = float(s.split()[6])
      if "slice" not in grid_desc:
        grid_desc["slice"] = []
      grid_desc["slice"].append(p)
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
    elif s.lower()[:9] == "read_grid" and len(s.split()) == 2:
      filename = s.split()[1]
      if not os.path.isfile(filename):
        print "Error reading SPARTA grid description file"
        print "read_grid filename is not available: ", filename
        sys.exit(1)
      else:
        grid_desc["read_grid"] = filename
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

def read_array_names(fh, array_names):
  for line in fh:
    s = clean_line(line)
    if s.lower().replace(" ", "")[:10] == "item:cells":
      for name in s.split()[2:]:
        array_names.append(name)
      break

def read_time_step_data(time_step_file_list, ug, id_hash):
  for f in time_step_file_list:
    try:
      fh = open(f, "r")
    except IOError:
      print "Unable to open SPARTA result file: ", f
      return

    array_names = []
    read_array_names(fh, array_names)

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
        print "Error reading SPARTA result file: ", f
        print "Flow data line cannot be processed:  ", line
        return

    fh.close()

def write_pvd_file(time_steps_dict, file_name, num_chunks):
  fh = open(file_name + ".pvd", "w")
  fh.write('<?xml version="1.0"?>\n')
  fh.write('<VTKFile type="Collection" version="0.1"\n')
  fh.write('               byte_order="LittleEndian"\n')
  fh.write('               compressor="vtkZLibDataCompressor">\n')
  fh.write('   <Collection>    \n')
  for time in sorted(time_steps_dict.keys()):
      fh.write('    <DataSet timestep="' + str(time) + '" group="" part="0"   \n')
      filepath = os.path.join(file_name, file_name + '_'  + str(time) + '.pvtu')
      fh.write('             file="' + filepath + '"/>\n')
      file_list = time_steps_dict[time]
      if file_list:
        try:
          afh = open(file_list[0], "r")
        except IOError:
          print "Unable to open SPARTA result file: ", f
          return
      else:
          return
      array_names = []
      read_array_names(afh, array_names)
      afh.close()
      write_pvtu_file(array_names, file_name, num_chunks, time)
  fh.write('   </Collection>    \n')
  fh.write('</VTKFile>    \n')
  fh.close()

def write_pvtu_file(array_names, output_file, num_chunks, time):
  filepath = os.path.join(output_file, output_file + '_'  + str(time) + '.pvtu')
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

  for chunk_id in range(num_chunks):
    fh.write('<Piece Source="' + output_file + '_' + str(chunk_id) + '_' + \
      str(time) + '.vtu"/>\n')

  fh.write('</PUnstructuredGrid>\n')
  fh.write('</VTKFile>\n')

def write_slice_pvd_file(time_steps_dict, output_file):
  fh = open(output_file + ".pvd", "w")
  fh.write('<?xml version="1.0"?>\n')
  fh.write('<VTKFile type="Collection" version="0.1"\n')
  fh.write('               byte_order="LittleEndian"\n')
  fh.write('               compressor="vtkZLibDataCompressor">\n')
  fh.write('   <Collection>    \n')
  for time in sorted(time_steps_dict.keys()):
      fh.write('    <DataSet timestep="' + str(time) + '" group="" part="0"   \n')
      filepath = os.path.join(output_file, output_file + '_'  + \
        str(time) + '.pvtu')
      fh.write('             file="' + filepath + '"/>\n')
      file_list = time_steps_dict[time]
      if file_list:
        try:
          afh = open(file_list[0], "r")
        except IOError:
          print "Unable to open SPARTA result file: ", f
          return
      else:
          return
      array_names = []
      read_array_names(afh, array_names)
      afh.close()
      write_slice_pvtu_file(array_names, output_file, time)
  fh.write('   </Collection>    \n')
  fh.write('</VTKFile>    \n')
  fh.close()

def write_slice_pvtu_file(array_names, output_file, time):
  filepath = os.path.join(output_file, output_file + '_' + str(time) + '.pvtu')
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

  rexp_path = os.path.join(output_file, output_file + '_*_' + str(time) + '.vtu')
  file_list = glob.glob(rexp_path)

  for f in file_list:
    fh.write('<Piece Source="' + os.path.basename(f) +'"/>\n')

  fh.write('</PUnstructuredGrid>\n')
  fh.write('</VTKFile>\n')

class ChunkReport:
  def __init__(self):
    controller = vtk.vtkMultiProcessController.GetGlobalController()
    self.num_procs = controller.GetNumberOfProcesses()

  def reportChunkComplete(self, chunk_id):
    if self.num_procs == 1:
      print "Completed grid chunk number: ", chunk_id

class ParallelTimer:
  def __init__(self):
    from mpi4py import MPI
    if MPI.Is_initialized():
      self.comm = MPI.COMM_WORLD
      self.rank = self.comm.Get_rank()
      self.size = self.comm.Get_size()
      self.start_time = time.time()
    else:
      self.size = 1
      self.rank = 0

  def report_rank_zero_time(self, message):
    if self.rank == 0:
      time_secs = time.time() - self.start_time
      print ""
      print message % timedelta(seconds=round(time_secs))
      print ""

  def report_collective_time(self, message):
    if self.size == 1:
      return
    from mpi4py import MPI
    time_secs = time.time() - self.start_time
    total_time = self.comm.reduce(time_secs, op=MPI.SUM)
    max_time = self.comm.reduce(time_secs, op=MPI.MAX)
    min_time = self.comm.reduce(time_secs, op=MPI.MIN)
    if self.rank == 0:
      print ""
      print "Average " + message % timedelta(seconds=total_time/float(self.size))
      print "Maxiumum " + message % timedelta(seconds=max_time)
      print "Minimum " + message % timedelta(seconds=min_time)
      print ""

def report_collective_grid_sizes(ug):
  from mpi4py import MPI
  num_cells = ug.GetNumberOfCells()
  mem_size = ug.GetActualMemorySize()
  if MPI.Is_initialized():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    total_cell_count = comm.reduce(num_cells, op=MPI.SUM)
    max_cell_count = comm.reduce(num_cells, op=MPI.MAX)
    min_cell_count = comm.reduce(num_cells, op=MPI.MIN)
    total_mem_used = comm.reduce(mem_size, op=MPI.SUM)
    min_mem_used = comm.reduce(mem_size, op=MPI.MIN)
    max_mem_used = comm.reduce(mem_size, op=MPI.MAX)
    if rank == 0:
      print "Average grid cell count over MPI ranks: {:.1e}"\
        .format(total_cell_count/float(size))
      print "Minimum grid cell count over MPI ranks: {:.1e}"\
        .format(min_cell_count)
      print "Maximum grid cell count over MPI ranks: {:.1e}"\
        .format(max_cell_count)
      print "Average grid memory used over MPI ranks: {} MB"\
        .format((total_mem_used/float(size))/1000.0)
      print "Minimum grid memory used over MPI ranks: {} MB"\
        .format(min_mem_used/1000.0)
      print "Maximum grid memory used over MPI ranks: {} MB"\
        .format(max_mem_used/1000.0)

def setup_for_MPI(params_dict):
  from mpi4py import MPI
  pd = params_dict
  comm = MPI.COMM_WORLD
  pd["size"] = comm.Get_size()
  pd["rank"] = comm.Get_rank()

  pd["chunking"] = comm.bcast(pd["chunking"], root=0)
  pd["grid_desc"] = comm.bcast(pd["grid_desc"], root=0)
  pd["time_steps_dict"] = comm.bcast(pd["time_steps_dict"], root=0)
  pd["paraview_output_file"] = comm.bcast(pd["paraview_output_file"], root=0)
  pd["catalystscript"] = comm.bcast(pd["catalystscript"], root=0)

def run_pvbatch_output(params_dict):
  rank =  params_dict["rank"]
  size =  params_dict["size"]
  chunking =  params_dict["chunking"]
  grid_desc = params_dict["grid_desc"]
  time_steps_dict = params_dict["time_steps_dict"]
  paraview_output_file = params_dict["paraview_output_file"]
  catalystscript = params_dict["catalystscript"]

  if rank == 0:
    print "Processing grid chunk(s) on " + str(size) + " MPI ranks"

  c = len(chunking)/size
  r = len(chunking) % size
  if rank < r:
    start = rank * (c + 1)
    stop = start + c
  else:
    start = rank * c + r
    stop = start + (c - 1)

  append = vtk.vtkAppendFilter()
  append.MergePointsOn()
  for idx, chunk in enumerate(chunking[start:stop+1]):
    g = process_grid_chunk(idx, chunk, len(chunking), \
      grid_desc, time_steps_dict, paraview_output_file)
    if g:
      append.AddInputData(g)
  ug = None
  if append.GetInputList().GetNumberOfItems():
    append.Update()
    ug = append.GetOutput()

  if ug is None:
    ug = vtk.vtkUnstructuredGrid()

  report_collective_grid_sizes(ug)

  if catalystscript is not None:
    if rank == 0:
      print "Calling Catalyst over " + str(len(time_steps_dict)) + " time step(s) ..."
    import coprocessor
    coprocessor.initialize()
    coprocessor.addscript(catalystscript)
    id_map = create_cell_global_id_to_local_id_map(ug)
    for idx, time in enumerate(sorted(time_steps_dict.keys())):
      pt = ParallelTimer()
      read_time_step_data(time_steps_dict[time], ug, id_map)
      coprocessor.coprocess(time, idx, ug, paraview_output_file + '.pvd')
      pt.report_collective_time("catalyst output time %s (wall clock time) for step "\
        + str(idx))
    coprocessor.finalize()
  else:
    if rank == 0:
      print "Writing grid files over " + str(len(time_steps_dict)) + " time step(s) ..."
    write_grid_chunk(ug, rank, size, grid_desc, time_steps_dict, \
      paraview_output_file)

if __name__ == "__main__":
  controller = vtk.vtkMultiProcessController.GetGlobalController()
  num_procs = controller.GetNumberOfProcesses()
  local_proc_id = controller.GetLocalProcessId()

  if local_proc_id == 0:
    pt = ParallelTimer()
    parser = argparse.ArgumentParser()
    parser.add_argument("sparta_grid_description_file", help="SPARTA grid description input file name")
    parser.add_argument("paraview_output_file", help="ParaView output file name")
    parser.add_argument('-c', '--catalystscript', help="Run a ParaView Catalyst Python script with pvbatch")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-r', '--result', help="Optional list of SPARTA dump result files", nargs='+')
    group.add_argument('-rf', '--resultfile', \
                        help="Optional filename containing path names of SPARTA dump result files")
    parser.add_argument('-xc', '--xchunk', \
                        help="Optional x grid chunk size (positive integer; default 100)", \
                        default=100, type=int)
    parser.add_argument('-yc', '--ychunk', \
                        help="Optional y grid chunk size (positive integer; default 100)", \
                        default=100, type=int)
    parser.add_argument('-zc', '--zchunk', \
                        help="Optional z grid chunk size (positive integer; default 100)", \
                        default=100, type=int)
    args = parser.parse_args()

    try:
      gdf = open(args.sparta_grid_description_file, "r")
    except IOError:
      print "Unable to open SPARTA surf input file: ", args.sparta_grid_description_file
      sys.exit(1)

    if os.path.isdir(args.paraview_output_file) and not args.catalystscript:
      print "ParaView output directory exists: ", args.paraview_output_file
      sys.exit(1)
 
    if args.xchunk < 1:
      print "Invalid xchunk size given: ", args.xchunk
      sys.exit(1)

    if args.ychunk < 1:
      print "Invalid ychunk size given: ", args.ychunk
      sys.exit(1)

    if args.zchunk < 1:
      print "Invalid zchunk size given: ", args.zchunk
      sys.exit(1)

    if args.catalystscript:
      if num_procs == 1:
        print "Error: ParaView Catalyst only available with pvbatch and more than one MPI rank"
        sys.exit(1)
      try:
        cf = open(args.catalystscript, "r")
        cf.close()
      except IOError:
        print "Error: Unable to open Catalyst Python file: ", args.catalystscript
        sys.exit(1)

    grid_desc = {}
    read_grid_description_file(gdf, grid_desc)
    gdf.close()

    if "dimension" not in grid_desc:
      print "Error: grid description file does not have a dimension statement: ", \
        args.sparta_grid_description_file
      sys.exit(1)

    if "create_box" not in grid_desc:
      print "Error: grid description file does not have a create_box statement: ", \
        args.sparta_grid_description_file
      sys.exit(1)

    if "read_grid" not in grid_desc and "create_grid" not in grid_desc:
      print "Error: grid description file does not have a read_grid or a create_grid statement: ", \
        args.sparta_grid_description_file
      sys.exit(1)

    if "slice" in grid_desc and args.catalystscript:
      print "Error: Slice plane output not available with ParaView Catalyst"
      sys.exit(1)

    if "slice" not in grid_desc:
      if os.path.isfile(args.paraview_output_file + '.pvd') and not \
         args.catalystscript:
        print "ParaView output file exists: ", args.paraview_output_file + '.pvd'
        sys.exit(1)
    else:
      for idx, slice in enumerate(grid_desc["slice"]):
        file_name = args.paraview_output_file + '_slice' + str(idx) + "-" +\
                      str(round(slice["nx"],4)) + "_" + \
                      str(round(slice["ny"],4)) + "_" + \
                      str(round(slice["nz"],4))
        if os.path.isfile(file_name + '.pvd'):
          print "ParaView output file exists: ", file_name + '.pvd'
          sys.exit(1)

    if "read_grid" in grid_desc:
      create_grid_from_grid_file(grid_desc)

    time_steps_dict = {}
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
        print "Unable to open SPARTA result file input list file: ", args.result_file
        sys.exit(1)

    if not time_steps_file_list:
      time_steps_dict[0] = []

    read_time_steps(time_steps_file_list, time_steps_dict)

    chunking = []
    find_chunking(chunking, grid_desc, args)

    sys.stdin = open(os.devnull)

    print "Processing ", len(chunking), " grid chunk(s)."

    if not args.catalystscript:
      os.mkdir(args.paraview_output_file)

  import platform
  if platform.system() == 'Linux' or platform.system() == 'Darwin':
    if num_procs == 1:
      import multiprocessing as mp
      pool = mp.Pool()
      for idx, chunk in enumerate(chunking):
        pool.apply_async(create_and_write_grid_chunk, \
                         args=(idx, chunk, len(chunking), grid_desc, \
                               time_steps_dict, args.paraview_output_file, ))
      pool.close()
      pool.join()
    else:
      pd = {}
      if local_proc_id == 0:
        pd["catalystscript"] = args.catalystscript
        pd["paraview_output_file"] = args.paraview_output_file
        pd["chunking"] = chunking
        pd["grid_desc"] = grid_desc
        pd["time_steps_dict"] = time_steps_dict
      else:
        pd["catalystscript"] = None
        pd["paraview_output_file"] = None
        pd["chunking"] = None
        pd["grid_desc"] = None
        pd["time_steps_dict"] = None

      setup_for_MPI(pd)
      run_pvbatch_output(pd)
  else:
    for idx, chunk in enumerate(chunking):
      create_and_write_grid_chunk(idx, chunk, len(chunking), grid_desc, time_steps_dict,
                                  args.paraview_output_file)
  if local_proc_id == 0:
    if args.catalystscript is None:
      if "slice" not in grid_desc:
        if num_procs == 1:
          write_pvd_file(time_steps_dict, args.paraview_output_file, len(chunking))
        else:
          write_pvd_file(time_steps_dict, args.paraview_output_file, num_procs)
      else:
        write_slice_pvd_file(time_steps_dict, args.paraview_output_file)
          
    pt.report_rank_zero_time("Done in %s (wall clock time)")
