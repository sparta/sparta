coProcessor = None

def initialize():
    global coProcessor
    import paraview
    paraview.options.batch = True
    paraview.options.symmetric = True
    from paraview.vtk import vtkPVCatalyst as catalyst
    coProcessor = catalyst.vtkCPProcessor()

def finalize():
    global coProcessor
    coProcessor.Finalize()

def addscript(name):
    global coProcessor
    from paraview.vtk import vtkPVPythonCatalystPython as pythoncatalyst
    pipeline = pythoncatalyst.vtkCPPythonScriptPipeline()
    pipeline.Initialize(name)
    coProcessor.AddPipeline(pipeline)

def coprocess(time, timeStep, grid, inputName):
    global coProcessor
    from paraview.vtk import vtkPVCatalyst as catalyst
    dataDescription = catalyst.vtkCPDataDescription()
    dataDescription.SetTimeData(time, timeStep)
    dataDescription.AddInput(inputName)

    if coProcessor.RequestDataDescription(dataDescription):
        dataDescription.GetInputDescriptionByName(inputName).SetGrid(grid)
        coProcessor.CoProcess(dataDescription)
