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
renderView1.ViewSize = [844, 539]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [10.282934247554099, 5.721930840307251, 6.376457439874225]
renderView1.CameraViewUp = [-0.3704888297091945, 0.9039254487437967, -0.2136745426438009]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 3.4641016151377544
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(844, 539)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVD Reader'
sphere_grid_pvpythonpvd = PVDReader(registrationName='sphere_grid_pvpython.pvd', FileName='/ascldap/users/tjotaha/sparta_development/sparta/build/examples/sphere/sphere_grid_pvpython.pvd')
sphere_grid_pvpythonpvd.CellArrays = ['id', 'f_1[1]', 'f_1[2]', 'f_1[3]', 'f_1[4]', 'f_1[5]', 'f_1[6]', 'f_1[7]', 'f_1[8]', 'f_1[9]', 'f_1[10]', 'f_1[11]', 'f_1[12]', 'f_1[13]', 'f_1[14]']

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=sphere_grid_pvpythonpvd)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['CELLS', 'f_1[10]']
clip1.Value = 45.29399999999998

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from clip1
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'f_15'
f_15LUT = GetColorTransferFunction('f_15')
f_15LUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 213793.5, 0.865003, 0.865003, 0.865003, 427587.0, 0.705882, 0.0156863, 0.14902]
f_15LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'f_15'
f_15PWF = GetOpacityTransferFunction('f_15')
f_15PWF.Points = [0.0, 0.0, 0.5, 0.0, 427587.0, 1.0, 0.5, 0.0]
f_15PWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
clip1Display.Representation = 'Surface With Edges'
clip1Display.ColorArrayName = ['CELLS', 'f_1[5]']
clip1Display.LookupTable = f_15LUT
clip1Display.SelectTCoordArray = 'None'
clip1Display.SelectNormalArray = 'None'
clip1Display.SelectTangentArray = 'None'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'None'
clip1Display.ScaleFactor = 0.4
clip1Display.SelectScaleArray = 'None'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'None'
clip1Display.GaussianRadius = 0.02
clip1Display.SetScaleArray = [None, '']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = [None, '']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityFunction = f_15PWF
clip1Display.ScalarOpacityUnitDistance = 0.37797631496846196
clip1Display.OpacityArrayName = ['CELLS', 'f_1[10]']

# setup the color legend parameters for each legend in this view

# get color legend/bar for f_15LUT in view renderView1
f_15LUTColorBar = GetScalarBar(f_15LUT, renderView1)
f_15LUTColorBar.Title = 'f_1[5]'
f_15LUTColorBar.ComponentTitle = ''

# set color bar visibility
f_15LUTColorBar.Visibility = 1

# show color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

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
pNG1.Writer.ImageResolution = [844, 539]
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
