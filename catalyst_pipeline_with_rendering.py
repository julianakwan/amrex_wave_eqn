# script-version: 2.0
from paraview.simple import *
from paraview import catalyst
import time

#paraview.simple._DisableFirstRenderCameraReset()

# registrationName must match the channel name used in the
# 'CatalystAdaptor'.
producer = PVTrivialProducer(registrationName="input")


# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [800,800]
renderView1.CameraPosition = [50., 50., 200.0]
renderView1.CameraFocalPoint = [50.0, 50.0, 50.0]
renderView1.CameraParallelScale = 86.60254037844386

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=producer)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [50., 50., 50.]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [50., 50., 50.]

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 0.0, -1.0]


# show data from grid
#gridDisplay = Show(slice1, renderView1, 'UniformGridRepresentation')
gridDisplay = Show(slice1, renderView1, 'AMRRepresentation')

# trace defaults for the display properties.
gridDisplay.Representation = 'Surface'
gridDisplay.ColorArrayName = [None, '']
gridDisplay.SelectTCoordArray = 'None'
gridDisplay.SelectNormalArray = 'None'
gridDisplay.SelectTangentArray = 'None'
gridDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
gridDisplay.SelectOrientationVectors = 'None'
gridDisplay.ScaleFactor = 0.1
gridDisplay.SelectScaleArray = 'None'
gridDisplay.GlyphType = 'Arrow'
gridDisplay.GlyphTableIndexArray = 'None'
gridDisplay.GaussianRadius = 0.005
gridDisplay.SetScaleArray = [None, '']
gridDisplay.ScaleTransferFunction = 'PiecewiseFunction'
gridDisplay.OpacityArray = [None, '']
gridDisplay.OpacityTransferFunction = 'PiecewiseFunction'
gridDisplay.DataAxesGrid = 'GridAxesRepresentation'
gridDisplay.PolarAxes = 'PolarAxesRepresentation'
gridDisplay.ScalarOpacityUnitDistance = 0.0065053373254344085


# get color transfer function/color map for 'phi0'
phi0_LUT = GetColorTransferFunction('dphi0')
gridDisplay.LookupTable = phi0_LUT


# rescale color and/or opacity maps used to include current data range
gridDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
gridDisplay.SetScalarBarVisibility(renderView1, True)


# get opacity transfer function/opacity map for 'phi0'
phi0_PWF = GetOpacityTransferFunction('phi0')

# get 2D transfer function for 'phi0'
phi0_TF2D = GetTransferFunction2D('phi0')


# get color legend/bar for velocityLUT in view renderView1
phi0_LUTColorBar = GetScalarBar(phi0_LUT, renderView1)
phi0_LUTColorBar.Title = 'phi0'

# set color bar visibility
phi0_LUTColorBar.Visibility = 1

# set scalar coloring
ColorBy(gridDisplay, ('CELLS', 'phi0'))

# Rescale transfer function
phi0_LUT.RescaleTransferFunction(-3.5, 3.5)

# Rescale transfer function
phi0_PWF.RescaleTransferFunction(-3.5, 3.5)

# Rescale 2D transfer function
phi0_TF2D.RescaleTransferFunction(-3.5, 3.5, 0.0, 1.0)
# show color legend
gridDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

SetActiveView(renderView1)
# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
pNG1.Trigger = 'TimeStep'

# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'test_{timestep:06d}.png'
pNG1.Writer.ImageResolution = [1600,800]
pNG1.Writer.Format = 'PNG'

SetActiveSource(pNG1)

# ------------------------------------------------------------------------------
# Catalyst options
options = catalyst.Options()
options.GlobalTrigger = 'TimeStep'
options.ExtractsOutputDirectory = "/home/dc-kwan1/rds/rds-dirac-dp002/dc-kwan1/AMReX/wave/insitu-viz/paraview_5p12"

if "--enable-live" in catalyst.get_args():
  options.EnableCatalystLive = 1
  options.CatalystLiveTrigger = 'TimeStep'
  options.CatalystLiveURL = "localhost:22222"


# Greeting to ensure that ctest knows this script is being imported
print("executing catalyst_pipeline")

def catalyst_execute(info):
    global producer, grid
    producer.UpdatePipeline()
    print(producer.GetDataInformation().GetDataAssembly())


    print("-----------------------------------")
    print("executing (cycle={}, time={})".format(info.cycle, info.time))

    SaveExtractsUsingCatalystOptions(options)

    # In a real simulation sleep is not needed. We use it here to slow down the
    # "simulation" and make sure ParaView client can catch up with the produced
    # results instead of having all of them flashing at once.
    if options.EnableCatalystLive:
        time.sleep(1)
