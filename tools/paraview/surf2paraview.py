from __future__ import print_function
#   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
#   http://sparta.sandia.gov
#   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov,
#   Thomas Otahal, tjotaha@sandia.gov
#   Sandia National Laboratories

#   Copyright (2014) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under 
#   the GNU General Public License.

#   See the README file in the top-level SPARTA directory.

import argparse
import sys
import os
import vtk
import glob
from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase

class UGridSource(VTKPythonAlgorithmBase):
  def __init__(self, ug, time_steps_dict):
    VTKPythonAlgorithmBase.__init__(self,
      nInputPorts=0, nOutputPorts=1, outputType='vtkUnstructuredGrid')
    self.__ug = ug
    self.__time_steps_dict = time_steps_dict

  def RequestInformation(self, request, inInfo, outInfo):
    info = outInfo.GetInformationObject(0)
    tsteps = sorted(self.__time_steps_dict.keys())
    if tsteps:
      info.Set(vtk.vtkStreamingDemandDrivenPipeline.TIME_STEPS(),
        tsteps, len(tsteps))
      info.Set(vtk.vtkStreamingDemandDrivenPipeline.TIME_RANGE(),
        [tsteps[0], tsteps[-1]], 2)
    return 1

  def RequestData(self, request, inInfo, outInfo):
    info = outInfo.GetInformationObject(0)
    output = vtk.vtkUnstructuredGrid.GetData(outInfo)
    tstep = info.Get(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP())
    if tstep in self.__time_steps_dict:
      read_time_step_data(self.__time_steps_dict[tstep], self.__ug)
    output.ShallowCopy(self.__ug)
    return 1

def clean_line(line):
  line = line.partition('#')[0]
  return line.strip()

def read_points(sif, num_points, num_elements, ug, three_d_file):
  num_items_per_line = 4
  if not three_d_file:
    num_items_per_line = 3

  if ug.GetPoints():
    print("Error reading SPARTA surf input file")
    print("Points section of file occurs more than once")
    sys.exit(1)
  
  points = vtk.vtkPoints()
  for line in sif:
    s = clean_line(line)
    if s and len(s.split()) == num_items_per_line:
      x = float(s.split()[1])
      y = float(s.split()[2])
      z = 0.0
      if three_d_file:
        z = float(s.split()[3])
      points.InsertNextPoint(x,y,z)
    elif s and len(s.split()) == 1:
      if s.split()[0].lower() == 'triangles' or \
         s.split()[0].lower() == 'lines':
        read_elements(sif, num_points, num_elements, ug, three_d_file)
      else:
        print("Error reading SPARTA surf input file")
        print("File contains bad section header: ", s)
        sys.exit(1)
    elif s:
      print("Error reading SPARTA surf input file")
      print("Points section contains bad line: ", s)
      sys.exit(1)

  if points.GetNumberOfPoints() != num_points:
    print("Error reading SPARTA surf input file")
    print("The number of points in point section: ", points.GetNumberOfPoints())
    print("Does not agree with number of points in header: ", num_points)
    sys.exit(1)
  else:
    ug.SetPoints(points)

def read_elements(sif, num_points, num_elements, ug, three_d_file):
  num_items_per_line = 4
  if not three_d_file:
    num_items_per_line = 3

  if ug.GetNumberOfCells():
    print("Error reading SPARTA surf input file")
    print("Trangles or Lines section of file occurs more than once")
    sys.exit(1)

  for line in sif:
    s = clean_line(line)
    if s and len(s.split()) == num_items_per_line:
      if three_d_file:
        tri = vtk.vtkTriangle()
        i = int(s.split()[1])
        j = int(s.split()[2])
        k = int(s.split()[3])
        tri.GetPointIds().SetId(0,i-1)
        tri.GetPointIds().SetId(1,j-1)
        tri.GetPointIds().SetId(2,k-1)
        ug.InsertNextCell(tri.GetCellType(), tri.GetPointIds());
      else:
        li = vtk.vtkLine()
        i = int(s.split()[1])
        j = int(s.split()[2])
        li.GetPointIds().SetId(0,i-1)
        li.GetPointIds().SetId(1,j-1)
        ug.InsertNextCell(li.GetCellType(), li.GetPointIds());
    elif s and len(s.split()) == 1:
      if s.split()[0].lower() == 'points':
        read_points(sif, num_points, num_elements, ug, three_d_file)
      else:
        print("Error reading SPARTA surf input file")
        print("File contains bad section header: ", s)
        sys.exit(1)
    elif s:
      print("Error reading SPARTA surf input file")
      print("Triangles or Lines section contains bad line: ", s)
      sys.exit(1)

  if ug.GetNumberOfCells() != num_elements:
    print("Error reading SPARTA surf input file")
    print("The number of elements in Triangles or Lines section: ", ug.GetNumberOfCells())
    print("Does not agree with number of elements in header: ", num_elements)
    sys.exit(1)

def read_surf_file(sif, ug):
  # Skip first line of file
  next(sif)

  num_points = 0
  num_elements = 0
  three_d_file = True
  for line in sif:
    s = clean_line(line)
    if s and len(s.split()) == 2:
      if s.split()[1].lower() == 'points':
        num_points = int(s.split()[0]) 
      elif s.split()[1].lower() == 'triangles':
        num_elements = int(s.split()[0]) 
      elif s.split()[1].lower() == 'lines':
        num_elements = int(s.split()[0]) 
        three_d_file = False
    elif s:
      if num_points <= 0:
        print("Error reading SPARTA surf input file")
        print("Number of points are: ", num_points)
        sys.exit(1)

      if num_elements <= 0:
        print("Error reading SPARTA surf input file")
        print("Number of elements are: ", num_points)
        sys.exit(1)
 
      if len(s.split()) != 1:
        print("Error reading SPARTA surf input file")
        print("Missing section header: Triangles, Lines, or Points")
        sys.exit(1)
      elif s.split()[0].lower() == 'points':
        read_points(sif, num_points, num_elements, ug, three_d_file)
      elif s.split()[0].lower() == 'triangles' or \
           s.split()[0].lower() == 'lines':
        read_elements(sif, num_points, num_elements, ug, three_d_file)
  
  if not three_d_file:
    polygon = vtk.vtkPolygon()

def read_time_steps(result_file_list, time_steps_dict):
  for f in result_file_list:
    try:
      fh = open(f, "r")
    except IOError:
      print("Unable to open SPARTA result file: ", f)
      sys.exit(1)

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

def read_time_step_data(time_step_file_list, ug):
  for f in time_step_file_list:
    try:
      fh = open(f, "r")
    except IOError:
      print("Unable to open SPARTA result file: ", f)
      sys.exit(1)

    array_names = []
    for line in fh:
      s = clean_line(line)
      if s.lower().replace(" ", "")[:10] == "item:surfs":
        for name in s.split()[2:]: 
          array_names.append(name)
        break

    if not ug.GetCellData().GetNumberOfArrays():
      for name in array_names:
        array = vtk.vtkDoubleArray()
        array.SetName(name)
        array.SetNumberOfComponents(1)
        array.SetNumberOfTuples(ug.GetNumberOfCells())
        array.FillComponent(0, 0.0)
        ug.GetCellData().AddArray(array)

    if ug.GetCellData().GetNumberOfArrays() != len(array_names):
      print("Error reading SPARTA result file: ", f)
      print("Expected data columns:  ", ug.GetCellData().GetNumberOfArrays())
      print("Found data columns:  ", len(array_names))
      sys.exit(1)

    for line in fh:
      s = clean_line(line)
      if len(s.split()) == len(array_names):
        index = int(s.split()[0]) - 1
        if index < 0 or index >= ug.GetNumberOfCells():
          print("Error reading SPARTA result file: ", f)
          print("Surface index out of range: ", index)
          print("Number of expected surfaces: ", ug.GetNumberOfCells())
          sys.exit(1)
        
        for idx, val in enumerate(array_names):
          array = ug.GetCellData().GetArray(idx)
          array.SetValue(index, float(s.split()[idx]))
      else:
        print("Error reading SPARTA result file: ", f)
        print("Surf data line cannot be processed:  ", line)
        sys.exit(1)

    fh.close()


def write_pvd_file(time_steps, file_name):
  fh = open(file_name + ".pvd", "w")
  fh.write('<?xml version="1.0"?>\n')
  fh.write('<VTKFile type="Collection" version="0.1"\n')
  fh.write('               byte_order="LittleEndian"\n')
  fh.write('               compressor="vtkZLibDataCompressor">\n')
  fh.write('   <Collection>    \n')
  
  for time in time_steps:
    filepath = os.path.join(file_name, file_name + '_' + str(time) + '.vtu') 
    fh.write('    <DataSet timestep="' + str(time) + '" group="" part="0"   \n')
    fh.write('             file="' + filepath + '"/>\n')

  fh.write('   </Collection>    \n')
  fh.write('</VTKFile>    \n')
  fh.close()

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("sparta_surf_file", help="SPARTA surface geometry input file name")
  parser.add_argument("paraview_output_file", help="ParaView output file name")
  group = parser.add_mutually_exclusive_group()
  group.add_argument('-r', '--result', help="Optional list of SPARTA dump result files", nargs='+')
  group.add_argument('-f', '--resultfile', help="Optional filename containing path names of SPARTA dump result files")
  parser.add_argument('-e', '--exodus', default=False, action='store_true',
    help="Write output to Exodus II format file")
  args = parser.parse_args()

  try:
    sif = open(args.sparta_surf_file, "r")
  except IOError:
    print("Unable to open SPARTA surf input file: ", args.sparta_surf_file)
    sys.exit(1)

  if os.path.isfile(args.paraview_output_file + '.pvd') and not args.exodus:
    print("ParaView output file exists: ", args.paraview_output_file + '.pvd')
    sys.exit(1)

  if os.path.isdir(args.paraview_output_file) and not args.exodus:
    print("ParaView output directory exists: ", args.paraview_output_file)
    sys.exit(1)

  print("Processing SPARTA surface information.")

  ug = vtk.vtkUnstructuredGrid()
  read_surf_file(sif, ug)
  sif.close()

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
      print("Unable to open SPARTA result file input list file: ", args.result_file)
      sys.exit(1)

  if not time_steps_file_list:
    time_steps_dict[0] = []

  read_time_steps(time_steps_file_list, time_steps_dict)

  if args.exodus:
    writer = vtk.vtkExodusIIWriter()
    writer.WriteAllTimeStepsOn()
    writer.SetFileName(args.paraview_output_file + ".ex2")
    ugs = UGridSource(ug, time_steps_dict)
    writer.SetInputConnection(ugs.GetOutputPort())
    writer.Write()
  else:
    os.mkdir(args.paraview_output_file)
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(ug)
    for time in sorted(time_steps_dict.keys()):
      print("Processing dump result file: ", time_steps_dict[time])
      read_time_step_data(time_steps_dict[time], ug)
      filepath = os.path.join(args.paraview_output_file, args.paraview_output_file + '_' + str(time) + '.vtu')
      writer.SetFileName(filepath)
      writer.Write()
    write_pvd_file(sorted(time_steps_dict.keys()), args.paraview_output_file)
