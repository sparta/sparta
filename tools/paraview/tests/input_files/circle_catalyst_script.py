
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
      renderView1.InteractionMode = '2D'
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [5.0, 5.0, 0.0]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [5.0, 5.0, 10000.0]
      renderView1.CameraFocalPoint = [5.0, 5.0, 0.0]
      renderView1.CameraParallelScale = 7.332959212304938
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
      circle_grid_pvbatchpvd = coprocessor.CreateProducer(datadescription, 'circle_grid_pvbatch.pvd')

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from circle_grid_pvbatchpvd
      circle_grid_pvbatchpvdDisplay = Show(circle_grid_pvbatchpvd, renderView1)

      # get color transfer function/color map for 'f_12'
      f_12LUT = GetColorTransferFunction('f_12')
      f_12LUT.RGBPoints = [-27.991, 0.231373, 0.298039, 0.752941, 263.7225, 0.865003, 0.865003, 0.865003, 555.436, 0.705882, 0.0156863, 0.14902]
      f_12LUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'f_12'
      f_12PWF = GetOpacityTransferFunction('f_12')
      f_12PWF.Points = [-27.991, 0.0, 0.5, 0.0, 555.436, 1.0, 0.5, 0.0]
      f_12PWF.ScalarRangeInitialized = 1

      # trace defaults for the display properties.
      circle_grid_pvbatchpvdDisplay.Representation = 'Surface With Edges'
      circle_grid_pvbatchpvdDisplay.ColorArrayName = ['CELLS', 'f_1[2]']
      circle_grid_pvbatchpvdDisplay.LookupTable = f_12LUT
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
      circle_grid_pvbatchpvdDisplay.SelectionCellLabelFontFile = ''
      circle_grid_pvbatchpvdDisplay.SelectionPointLabelFontFile = ''
      circle_grid_pvbatchpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
      circle_grid_pvbatchpvdDisplay.ScalarOpacityFunction = f_12PWF
      circle_grid_pvbatchpvdDisplay.ScalarOpacityUnitDistance = 1.9193831036664846

      # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
      circle_grid_pvbatchpvdDisplay.DataAxesGrid.XTitleFontFile = ''
      circle_grid_pvbatchpvdDisplay.DataAxesGrid.YTitleFontFile = ''
      circle_grid_pvbatchpvdDisplay.DataAxesGrid.ZTitleFontFile = ''
      circle_grid_pvbatchpvdDisplay.DataAxesGrid.XLabelFontFile = ''
      circle_grid_pvbatchpvdDisplay.DataAxesGrid.YLabelFontFile = ''
      circle_grid_pvbatchpvdDisplay.DataAxesGrid.ZLabelFontFile = ''

      # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
      circle_grid_pvbatchpvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
      circle_grid_pvbatchpvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
      circle_grid_pvbatchpvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
      circle_grid_pvbatchpvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for f_12LUT in view renderView1
      f_12LUTColorBar = GetScalarBar(f_12LUT, renderView1)
      f_12LUTColorBar.Title = 'f_1[2]'
      f_12LUTColorBar.ComponentTitle = ''
      f_12LUTColorBar.TitleFontFile = ''
      f_12LUTColorBar.LabelFontFile = ''

      # set color bar visibility
      f_12LUTColorBar.Visibility = 1

      # show color legend
      circle_grid_pvbatchpvdDisplay.SetScalarBarVisibility(renderView1, True)

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(circle_grid_pvbatchpvd)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'circle_grid_pvbatch.pvd': [1]}
  coprocessor.SetUpdateFrequencies(freqs)
  if requestSpecificArrays:
    arrays = [['f_1[10]', 1], ['f_1[11]', 1], ['f_1[12]', 1], ['f_1[13]', 1], ['f_1[14]', 1], ['f_1[1]', 1], ['f_1[2]', 1], ['f_1[3]', 1], ['f_1[4]', 1], ['f_1[5]', 1], ['f_1[6]', 1], ['f_1[7]', 1], ['f_1[8]', 1], ['f_1[9]', 1], ['id', 1]]
    coprocessor.SetRequestedArrays('circle_grid_pvbatch.pvd', arrays)
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
