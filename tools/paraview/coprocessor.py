
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

coProcessor = None

def initialize():
    global coProcessor
    import paraview
    paraview.options.batch = True
    paraview.options.symmetric = True
    from paraview.modules import vtkPVCatalyst as catalyst
    coProcessor = catalyst.vtkCPProcessor()

def finalize():
    global coProcessor
    coProcessor.Finalize()

def addscript(name):
    global coProcessor
    from paraview.modules import vtkPVPythonCatalyst as pythoncatalyst
    pipeline = pythoncatalyst.vtkCPPythonScriptPipeline()
    pipeline.Initialize(name)
    coProcessor.AddPipeline(pipeline)

def coprocess(time, timeStep, grid, inputName):
    global coProcessor
    from paraview.modules import vtkPVCatalyst as catalyst
    dataDescription = catalyst.vtkCPDataDescription()
    dataDescription.SetTimeData(time, timeStep)
    dataDescription.AddInput(inputName)

    if coProcessor.RequestDataDescription(dataDescription):
        dataDescription.GetInputDescriptionByName(inputName).SetGrid(grid)
        coProcessor.CoProcess(dataDescription)
