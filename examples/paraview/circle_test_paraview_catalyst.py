# script-version: 2.0
# Catalyst state generated using paraview version 5.10.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1132, 679]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [5.0, 5.0, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [5.0, 5.0, 27.320508075688775]
renderView1.CameraFocalPoint = [5.0, 5.0, 0.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 7.0710678118654755
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1132, 679)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVD Reader'
circle_grid_pvbatchpvd = PVDReader(registrationName='circle_grid_pvbatch.pvd', FileName='/data/tjotaha/src/sparta_work/test_suite/sparta/build/examples/circle/circle_grid_pvbatch.pvd')
circle_grid_pvbatchpvd.CellArrays = ['id', 'f_1[1]', 'f_1[2]', 'f_1[3]', 'f_1[4]', 'f_1[5]', 'f_1[6]', 'f_1[7]', 'f_1[8]', 'f_1[9]', 'f_1[10]', 'f_1[11]', 'f_1[12]', 'f_1[13]', 'f_1[14]']

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from circle_grid_pvbatchpvd
circle_grid_pvbatchpvdDisplay = Show(circle_grid_pvbatchpvd, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'f_17'
f_17LUT = GetColorTransferFunction('f_17')
f_17LUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 170186.0, 0.865003, 0.865003, 0.865003, 340372.0, 0.705882, 0.0156863, 0.14902]
f_17LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'f_17'
f_17PWF = GetOpacityTransferFunction('f_17')
f_17PWF.Points = [0.0, 0.0, 0.5, 0.0, 340372.0, 1.0, 0.5, 0.0]
f_17PWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
circle_grid_pvbatchpvdDisplay.Representation = 'Surface With Edges'
circle_grid_pvbatchpvdDisplay.ColorArrayName = ['CELLS', 'f_1[7]']
circle_grid_pvbatchpvdDisplay.LookupTable = f_17LUT
circle_grid_pvbatchpvdDisplay.SelectTCoordArray = 'None'
circle_grid_pvbatchpvdDisplay.SelectNormalArray = 'None'
circle_grid_pvbatchpvdDisplay.SelectTangentArray = 'None'
circle_grid_pvbatchpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
circle_grid_pvbatchpvdDisplay.SelectOrientationVectors = 'None'
circle_grid_pvbatchpvdDisplay.SelectScaleArray = 'None'
circle_grid_pvbatchpvdDisplay.GlyphType = 'Arrow'
circle_grid_pvbatchpvdDisplay.GlyphTableIndexArray = 'None'
circle_grid_pvbatchpvdDisplay.GaussianRadius = 0.05
circle_grid_pvbatchpvdDisplay.SetScaleArray = [None, '']
circle_grid_pvbatchpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
circle_grid_pvbatchpvdDisplay.OpacityArray = [None, '']
circle_grid_pvbatchpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
circle_grid_pvbatchpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
circle_grid_pvbatchpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
circle_grid_pvbatchpvdDisplay.ScalarOpacityFunction = f_17PWF
circle_grid_pvbatchpvdDisplay.ScalarOpacityUnitDistance = 1.9193831036664846
circle_grid_pvbatchpvdDisplay.OpacityArrayName = ['CELLS', 'f_1[10]']

# setup the color legend parameters for each legend in this view

# get color legend/bar for f_17LUT in view renderView1
f_17LUTColorBar = GetScalarBar(f_17LUT, renderView1)
f_17LUTColorBar.Title = 'f_1[7]'
f_17LUTColorBar.ComponentTitle = ''

# set color bar visibility
f_17LUTColorBar.Visibility = 1

# show color legend
circle_grid_pvbatchpvdDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
pNG1.Trigger = 'TimeStep'

# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'RenderView1_{timestep:06d}{camera}.png'
pNG1.Writer.ImageResolution = [1132, 679]
pNG1.Writer.Format = 'PNG'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(pNG1)
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.GlobalTrigger = 'TimeStep'
options.CatalystLiveTrigger = 'TimeStep'

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
