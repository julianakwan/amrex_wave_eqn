# trace generated using paraview version 5.11.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
import glob
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

#get filenames
filenames = glob.glob("/home/dc-kwan1/rds/rds-dirac-dp002/dc-kwan1/AMReX/wave/sine_gordon_cell_quartic_interp/plt_cell_quartic_log_err*")
# create a new 'AMReX/BoxLib Grid Reader'


print(filenames)

plt00000 = AMReXBoxLibGridReader(registrationName='plt00000*', FileNames=filenames)
plt00000.CellArrayStatus = []

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on plt00000
plt00000.CellArrayStatus = ['dphi0', 'frac_error', 'phi0']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
plt00000Display = Show(plt00000, renderView1, 'AMRRepresentation')

# trace defaults for the display properties.
plt00000Display.Representation = 'Outline'
plt00000Display.ColorArrayName = [None, '']
plt00000Display.SelectTCoordArray = 'None'
plt00000Display.SelectNormalArray = 'None'
plt00000Display.SelectTangentArray = 'None'
plt00000Display.OSPRayScaleFunction = 'PiecewiseFunction'
plt00000Display.SelectOrientationVectors = 'None'
plt00000Display.ScaleFactor = 0.1
plt00000Display.SelectScaleArray = 'None'
plt00000Display.GlyphType = 'Arrow'
plt00000Display.GlyphTableIndexArray = 'None'
plt00000Display.GaussianRadius = 0.005
plt00000Display.SetScaleArray = [None, '']
plt00000Display.ScaleTransferFunction = 'PiecewiseFunction'
plt00000Display.OpacityArray = [None, '']
plt00000Display.OpacityTransferFunction = 'PiecewiseFunction'
plt00000Display.DataAxesGrid = 'GridAxesRepresentation'
plt00000Display.PolarAxes = 'PolarAxesRepresentation'
plt00000Display.ScalarOpacityUnitDistance = 0.013531646934131853

# reset view to fit data
renderView1.ResetCamera(False)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(plt00000Display, ('CELLS', 'frac_err'))

# rescale color and/or opacity maps used to include current data range
plt00000Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
plt00000Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'phi0'
LUT = GetColorTransferFunction('frac_err')

# get opacity transfer function/opacity map for 'phi0'
PWF = GetOpacityTransferFunction('frac_err')

# get 2D transfer function for 'phi0'
TF2D = GetTransferFunction2D('frac_err')

#Rescale transfer functions
LUT.RescaleTransferFunction(-2.0,1.0)
PWF.RescaleTransferFunction(-2.0,1.0)

# change representation type
plt00000Display.SetRepresentationType('Surface')

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=plt00000)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.5, 0.5, 0.5]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [0.5, 0.5, 0.5]

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 0.0, -1.0]

# show data in view
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['CELLS', 'frac_err']
slice1Display.LookupTable = LUT
slice1Display.SelectTCoordArray = 'None'
slice1Display.SelectNormalArray = 'None'
slice1Display.SelectTangentArray = 'None'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 0.1
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = 0.005
slice1Display.SetScaleArray = [None, '']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = [None, '']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'
slice1Display.SelectInputVectors = [None, '']
slice1Display.WriteLog = ''

# hide data in view
Hide(plt00000, renderView1)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on animationScene1
animationScene1.AnimationTime = 0.0015625

animationScene1.Play()

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(819, 539)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [50.0, 50.0, 384.60652149512318]
renderView1.CameraFocalPoint = [50., 50., 50.]
renderView1.CameraParallelScale = 86.60254037844386

#--------------------------------------------
# uncomment the following to render all views
#RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
#SaveScreenshot("plot.png")
SaveAnimation('/home/dc-kwan1/rds/rds-dirac-dp002/dc-kwan1/AMReX/wave/sine_gordon_cell_quartic_interp/frac_err.jpeg', GetActiveView(), FrameWindow=[1,100], FrameRate=1)
#AnimateReader(plt00000, filename="test_movie.jpeg")
