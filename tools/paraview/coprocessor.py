
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

from paraview.catalyst import bridge

def initialize():
    bridge.initialize()

def finalize():
    bridge.finalize()

def addscript(filename):
    bridge.add_pipeline(filename)

def coprocess(time, timeStep, grid, inputName):
    from paraview import vtk
    from paraview.modules import vtkPVCatalyst as catalyst
    from paraview.vtk.util import numpy_support

    coProcessor = bridge.coprocessor
    dataDescription = catalyst.vtkCPDataDescription()
    dataDescription.SetTimeData(time, timeStep)
    dataDescription.AddInput(inputName)

    if coProcessor.RequestDataDescription(dataDescription):
        dataDescription.GetInputDescriptionByName(inputName).SetGrid(grid)
        coProcessor.CoProcess(dataDescription)
