# script-version: 2.0
from paraview.simple import *
from paraview import catalyst
import time

# registrationName must match the channel name used in the
# 'CatalystAdaptor'.
producer = PVTrivialProducer(registrationName="input")


# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [800,800]
renderView1.CameraPosition = [0.5, 0.5, 3.8460652149512318]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.5]
renderView1.CameraParallelScale = 0.8660254037844386


# get color transfer function/color map for 'velocity'
#velocityLUT = GetColorTransferFunction('velocity')
#velocityLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 29.205000000000002, 0.865003, 0.865003, 0.865003, 58.410000000000004, 0.705882, 0.0156863, 0.14902]
#velocityLUT.ScalarRangeInitialized = 1.0

clip1 = Clip(registrationName="Clip1", Input=producer)
clip1.ClipType='Plane'
clip1.HyperTreeGridClipper='Plane'
clip1.Scalars = ['CELLS', 'phi0']


clip1.ClipType.Origin = [0.0,0.0,0.0]
clip1.ClipType.Normal = [0.0,1.0,0.0]

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=producer)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.5, 0.5, 0.5]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [0.5, 0.5, 0.5]

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 0.0, -1.0]


# show data from grid
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
phi0LUT = GetColorTransferFunction('phi0')
gridDisplay.LookupTable = phi0LUT


# rescale color and/or opacity maps used to include current data range
gridDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
gridDisplay.SetScalarBarVisibility(renderView1, True)


# get opacity transfer function/opacity map for 'phi0'
phi0_PWF = GetOpacityTransferFunction('phi0')

# get 2D transfer function for 'phi0'
phi0_TF2D = GetTransferFunction2D('phi0')


# get color legend/bar for velocityLUT in view renderView1
phi0_LUTColorBar = GetScalarBar(phi0LUT, renderView1)
phi0_LUTColorBar.Title = 'phi0'

# set color bar visibility
phi0_LUTColorBar.Visibility = 1

# set scalar coloring
ColorBy(gridDisplay, ('CELLS', 'phi0'))


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
pNG1.Writer.FileName = 'screenshot_{timestep:06d}.png'
pNG1.Writer.ImageResolution = [1600,800]
pNG1.Writer.Format = 'PNG'

SetActiveSource(pNG1)

# ------------------------------------------------------------------------------
# Catalyst options
options = catalyst.Options()
options.GlobalTrigger = 'TimeStep'
options.ExtractsOutputDirectory = "/home/dc-kwan1/rds/rds-dirac-dp002/dc-kwan1/AMReX/wave/insitu-viz/"

if "--enable-live" in catalyst.get_args():
  options.EnableCatalystLive = 1
  options.CatalystLiveTrigger = 'Timestep'


# Greeting to ensure that ctest knows this script is being imported
print("executing catalyst_pipeline")

def catalyst_execute(info):
    global producer, grid
    producer.UpdatePipeline()
    print(producer.GetDataInformation().GetDataAssembly())
#    gridInfo = producer.GetSubsetDataInformation(0, "//grid", "Hierarchy");

    print("-----------------------------------")
    print("executing (cycle={}, time={})".format(info.cycle, info.time))
#    print("bounds:", gridInfo.GetBounds())

    SaveExtractsUsingCatalystOptions(options)
#    print("phi0-range:", producer.PointData["phi0"].GetRange(-1))
#    print("fractional error-range:", producer.CellData["phi0"].GetRange(-1))
    # In a real simulation sleep is not needed. We use it here to slow down the
    # "simulation" and make sure ParaView client can catch up with the produced
    # results instead of having all of them flashing at once.
    if options.EnableCatalystLive:
        time.sleep(1)
