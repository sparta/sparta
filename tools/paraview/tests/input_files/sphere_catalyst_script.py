
#--------------------------------------------------------------

# Global timestep output options
timeStepToStartOutputAt=0
forceOutputAtFirstCall=False

# Global screenshot output options
imageFileNamePadding=0
rescale_lookuptable=False

# Whether or not to request specific arrays from the adaptor.
requestSpecificArrays=False

# a root directory under which all Catalyst output goes
rootDirectory=''

# makes a cinema D index table
make_cinema_table=False

#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# paraview version 5.6.0
#--------------------------------------------------------------

from paraview.simple import *
from paraview import coprocessing

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.6.0

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      # trace generated using paraview version 5.6.0
      #
      # To ensure correct image size when batch processing, please search 
      # for and uncomment the line `# renderView*.ViewSize = [*,*]`

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # get the material library
      materialLibrary1 = GetMaterialLibrary()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1512, 1568]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.StereoType = 0
      renderView1.CameraPosition = [11.693224453960685, 0.24893581001714482, 7.413155031681012]
      renderView1.CameraViewUp = [0.019443143110028937, 0.9977493262176635, -0.06417356323441484]
      renderView1.CameraParallelScale = 3.5924016749576713
      renderView1.Background = [0.32, 0.34, 0.43]
      renderView1.OSPRayMaterialLibrary = materialLibrary1

      # init the 'GridAxes3DActor' selected for 'AxesGrid'
      renderView1.AxesGrid.XTitleFontFile = ''
      renderView1.AxesGrid.YTitleFontFile = ''
      renderView1.AxesGrid.ZTitleFontFile = ''
      renderView1.AxesGrid.XLabelFontFile = ''
      renderView1.AxesGrid.YLabelFontFile = ''
      renderView1.AxesGrid.ZLabelFontFile = ''

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='RenderView1_%t.png', freq=1, fittoscreen=0, magnification=1, width=1512, height=1568, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # restore active view
      SetActiveView(renderView1)
      # ----------------------------------------------------------------

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'PVD Reader'
      # create a producer from a simulation input
      sphere_grid_pvbatchpvd = coprocessor.CreateProducer(datadescription, 'sphere_grid_pvbatch.pvd')

      # create a new 'Clip'
      clip1 = Clip(Input=sphere_grid_pvbatchpvd)
      clip1.ClipType = 'Plane'
      clip1.Scalars = ['CELLS', 'f_1[10]']

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from clip1
      clip1Display = Show(clip1, renderView1)

      # get color transfer function/color map for 'f_11'
      f_11LUT = GetColorTransferFunction('f_11')
      f_11LUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 4.555, 0.865003, 0.865003, 0.865003, 9.11, 0.705882, 0.0156863, 0.14902]
      f_11LUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'f_11'
      f_11PWF = GetOpacityTransferFunction('f_11')
      f_11PWF.Points = [0.0, 0.0, 0.5, 0.0, 9.11, 1.0, 0.5, 0.0]
      f_11PWF.ScalarRangeInitialized = 1

      # trace defaults for the display properties.
      clip1Display.Representation = 'Surface With Edges'
      clip1Display.ColorArrayName = ['CELLS', 'f_1[1]']
      clip1Display.LookupTable = f_11LUT
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
      clip1Display.SelectionCellLabelFontFile = ''
      clip1Display.SelectionPointLabelFontFile = ''
      clip1Display.PolarAxes = 'PolarAxesRepresentation'
      clip1Display.ScalarOpacityFunction = f_11PWF
      clip1Display.ScalarOpacityUnitDistance = 0.377976314968462

      # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
      clip1Display.DataAxesGrid.XTitleFontFile = ''
      clip1Display.DataAxesGrid.YTitleFontFile = ''
      clip1Display.DataAxesGrid.ZTitleFontFile = ''
      clip1Display.DataAxesGrid.XLabelFontFile = ''
      clip1Display.DataAxesGrid.YLabelFontFile = ''
      clip1Display.DataAxesGrid.ZLabelFontFile = ''

      # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
      clip1Display.PolarAxes.PolarAxisTitleFontFile = ''
      clip1Display.PolarAxes.PolarAxisLabelFontFile = ''
      clip1Display.PolarAxes.LastRadialAxisTextFontFile = ''
      clip1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for f_11LUT in view renderView1
      f_11LUTColorBar = GetScalarBar(f_11LUT, renderView1)
      f_11LUTColorBar.Title = 'f_1[1]'
      f_11LUTColorBar.ComponentTitle = ''
      f_11LUTColorBar.TitleFontFile = ''
      f_11LUTColorBar.LabelFontFile = ''

      # set color bar visibility
      f_11LUTColorBar.Visibility = 1

      # show color legend
      clip1Display.SetScalarBarVisibility(renderView1, True)

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(clip1)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'sphere_grid_pvbatch.pvd': [1, 1]}
  coprocessor.SetUpdateFrequencies(freqs)
  if requestSpecificArrays:
    arrays = [['f_1[10]', 1], ['f_1[11]', 1], ['f_1[12]', 1], ['f_1[13]', 1], ['f_1[14]', 1], ['f_1[1]', 1], ['f_1[2]', 1], ['f_1[3]', 1], ['f_1[4]', 1], ['f_1[5]', 1], ['f_1[6]', 1], ['f_1[7]', 1], ['f_1[8]', 1], ['f_1[9]', 1], ['id', 1]]
    coprocessor.SetRequestedArrays('sphere_grid_pvbatch.pvd', arrays)
  coprocessor.SetInitialOutputOptions(timeStepToStartOutputAt,forceOutputAtFirstCall)

  if rootDirectory:
      coprocessor.SetRootDirectory(rootDirectory)

  if make_cinema_table:
      coprocessor.EnableCinemaDTable()

  return coprocessor


#--------------------------------------------------------------
# Global variable that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView and the update frequency
coprocessor.EnableLiveVisualization(False, 1)

# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable,
        image_quality=0, padding_amount=imageFileNamePadding)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
