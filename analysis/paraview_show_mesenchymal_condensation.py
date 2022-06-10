#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
resultspvd = PVDReader(FileName='/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/mesenchymal_condensation/orientation/20220223-152730/0/results_from_time_0/results.pvd')

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1229, 612]

# show data in view
resultspvdDisplay = Show(resultspvd, renderView1)
# trace defaults for the display properties.
resultspvdDisplay.Representation = 'Surface'
resultspvdDisplay.ColorArrayName = [None, '']
resultspvdDisplay.OSPRayScaleArray = 'Ancestors'
resultspvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
resultspvdDisplay.SelectOrientationVectors = 'Ancestors'
resultspvdDisplay.ScaleFactor = 0.0
resultspvdDisplay.SelectScaleArray = 'Ancestors'
resultspvdDisplay.GlyphType = 'Arrow'
resultspvdDisplay.GlyphTableIndexArray = 'Ancestors'
resultspvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
resultspvdDisplay.PolarAxes = 'PolarAxesRepresentation'
resultspvdDisplay.ScalarOpacityUnitDistance = 13.590116619704348

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
resultspvdDisplay.OSPRayScaleFunction.Points = [-1.0, 0.0, 0.5, 0.0, 147.0, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Glyph'
glyph1 = Glyph(Input=resultspvd,
    GlyphType='Arrow')
glyph1.Scalars = ['POINTS', 'Ancestors']
glyph1.Vectors = ['POINTS', 'None']
glyph1.ScaleFactor = 0.0
glyph1.GlyphTransform = 'Transform2'

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.Visibility = 1

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XTitle = 'x [d]'
renderView1.AxesGrid.YTitle = 'y [d]'
renderView1.AxesGrid.ZTitle = 'z [d]'
renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XTitleFontSize = 48
renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YTitleFontSize = 48
renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZTitleFontSize = 48
renderView1.AxesGrid.GridColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelFontSize = 42
renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YLabelFontSize = 42
renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZLabelFontSize = 42

# Properties modified on renderView1
renderView1.Background = [1.0, 1.0, 1.0]

# Properties modified on glyph1
glyph1.GlyphType = 'Sphere'
glyph1.Scalars = ['POINTS', 'None']
glyph1.ScaleMode = 'scalar'
glyph1.ScaleFactor = 1.0
glyph1.GlyphMode = 'All Points'

# show data in view
glyph1Display = Show(glyph1, renderView1)
# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = [None, '']
glyph1Display.OSPRayScaleArray = 'Ancestors'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'Ancestors'
glyph1Display.ScaleFactor = 0.0
glyph1Display.SelectScaleArray = 'Ancestors'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'Ancestors'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
glyph1Display.OSPRayScaleFunction.Points = [-1.0, 0.0, 0.5, 0.0, 147.0, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(glyph1Display, ('POINTS', 'Ancestors'))

# rescale color and/or opacity maps used to include current data range
glyph1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Ancestors'
ancestorsLUT = GetColorTransferFunction('Ancestors')

# Properties modified on ancestorsLUT
ancestorsLUT.NumberOfTableValues = 258

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# reset view to fit data
renderView1.ResetCamera()

# create a new 'Plane'
plane1 = Plane()

# Properties modified on plane1
plane1.Point1 = [8.5, -0.5, 0.0]
plane1.Point2 = [-0.5, 11.5, 0.0]
plane1.XResolution = 8
plane1.YResolution = 11

# show data in view
plane1Display = Show(plane1, renderView1)
# trace defaults for the display properties.
plane1Display.Representation = 'Surface'
plane1Display.ColorArrayName = [None, '']
plane1Display.OSPRayScaleArray = 'Normals'
plane1Display.OSPRayScaleFunction = 'PiecewiseFunction'
plane1Display.SelectOrientationVectors = 'None'
plane1Display.ScaleFactor = 0.0
plane1Display.SelectScaleArray = 'None'
plane1Display.GlyphType = 'Arrow'
plane1Display.GlyphTableIndexArray = 'None'
plane1Display.DataAxesGrid = 'GridAxesRepresentation'
plane1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
plane1Display.OSPRayScaleFunction.Points = [-1.0, 0.0, 0.5, 0.0, 147.0, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# change representation type
plane1Display.SetRepresentationType('Surface With Edges')

# Properties modified on plane1Display
plane1Display.Opacity = 0.25

# Properties modified on plane1Display
plane1Display.EdgeColor = [0.0, 0.0, 0.0]

# create a new 'Plane'
plane2 = Plane()

# Properties modified on plane2
plane2.Origin = [-0.5, -0.5, 3.5]
plane2.Point1 = [8.5, -0.5, 3.5]
plane2.Point2 = [-0.5, 11.5, 3.5]
plane2.XResolution = 8
plane2.YResolution = 11

# show data in view
plane2Display = Show(plane2, renderView1)
# trace defaults for the display properties.
plane2Display.Representation = 'Surface'
plane2Display.ColorArrayName = [None, '']
plane2Display.OSPRayScaleArray = 'Normals'
plane2Display.OSPRayScaleFunction = 'PiecewiseFunction'
plane2Display.SelectOrientationVectors = 'None'
plane2Display.ScaleFactor = 0.0
plane2Display.SelectScaleArray = 'None'
plane2Display.GlyphType = 'Arrow'
plane2Display.GlyphTableIndexArray = 'None'
plane2Display.DataAxesGrid = 'GridAxesRepresentation'
plane2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
plane2Display.OSPRayScaleFunction.Points = [-1.0, 0.0, 0.5, 0.0, 147.0, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# change representation type
plane2Display.SetRepresentationType('Surface With Edges')

# Properties modified on plane2Display
plane2Display.Opacity = 0.25

# Properties modified on plane2Display
plane2Display.EdgeColor = [0.0, 0.0, 0.0]

# set active source
SetActiveSource(glyph1)

# create a new 'Threshold'
threshold1 = Threshold(Input=glyph1)
threshold1.Scalars = ['POINTS', 'Ancestors']
threshold1.ThresholdRange = [0.0, 95.0]

# Properties modified on threshold1
threshold1.ThresholdRange = [0.0, 47.0]

# get opacity transfer function/opacity map for 'Ancestors'
ancestorsPWF = GetOpacityTransferFunction('Ancestors')

# show data in view
threshold1Display = Show(threshold1, renderView1)
# trace defaults for the display properties.
threshold1Display.Representation = 'Surface'
threshold1Display.ColorArrayName = ['POINTS', 'Ancestors']
threshold1Display.LookupTable = ancestorsLUT
threshold1Display.OSPRayScaleArray = 'Ancestors'
threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
threshold1Display.SelectOrientationVectors = 'Ancestors'
threshold1Display.ScaleFactor = 0.0
threshold1Display.SelectScaleArray = 'Ancestors'
threshold1Display.GlyphType = 'Arrow'
threshold1Display.GlyphTableIndexArray = 'Ancestors'
threshold1Display.DataAxesGrid = 'GridAxesRepresentation'
threshold1Display.PolarAxes = 'PolarAxesRepresentation'
threshold1Display.ScalarOpacityFunction = ancestorsPWF
threshold1Display.ScalarOpacityUnitDistance = 0.6967045325829462

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
threshold1Display.OSPRayScaleFunction.Points = [-1.0, 0.0, 0.5, 0.0, 147.0, 1.0, 0.5, 0.0]

# hide data in view
Hide(glyph1, renderView1)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# hide color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, False)

# reset view to fit data
renderView1.ResetCamera()

animationScene1.GoToLast()

animationScene1.GoToFirst()

# hide data in view
Hide(threshold1, renderView1)

# set active source
SetActiveSource(glyph1)

# show data in view
glyph1Display = Show(glyph1, renderView1)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# set active source
SetActiveSource(resultspvd)

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.ZAxisUseCustomLabels = 1
renderView1.AxesGrid.ZAxisLabels = [0.0, 4.0]

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.YAxisUseCustomLabels = 1
renderView1.AxesGrid.YAxisLabels = [0.0, 10.0]

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XAxisUseCustomLabels = 1
renderView1.AxesGrid.XAxisLabels = [0.0, 2.0, 6.0, 8.0]

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XAxisLabels = [0.0, 8.0]

# current camera placement for renderView1
renderView1.CameraPosition = [-12.851045891633762, 28.865497695574422, 11.258452192395092]
renderView1.CameraFocalPoint = [3.9931107759475735, 5.48482006788253, 1.8251452594995499]
renderView1.CameraViewUp = [0.20844682580113302, -0.23306404468433012, 0.949858448343127]
renderView1.CameraParallelScale = 7.847677613477703

# save screenshot
SaveScreenshot('/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/plots/initial_mesenchymal_condensation.png', renderView1, ImageResolution=[2458, 1224],
    TransparentBackground=1)

# reset view to fit data
renderView1.ResetCamera()

# current camera placement for renderView1
renderView1.CameraPosition = [4.011303067207336, 5.484820067882538, -34.91278801761328]
renderView1.CameraFocalPoint = [4.011303067207336, 5.484820067882538, 1.8251452594995499]
renderView1.CameraParallelScale = 7.858245297373936

# save screenshot
SaveScreenshot('/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/plots/initial_mesenchymal_condensation_xy.png', renderView1, ImageResolution=[2458, 1224],
    TransparentBackground=1)

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# current camera placement for renderView1
renderView1.CameraPosition = [4.011303067207336, 35.84674839607496, 1.8251452594995499]
renderView1.CameraFocalPoint = [4.011303067207336, 5.484820067882538, 1.8251452594995499]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 7.858245297373936

# save screenshot
SaveScreenshot('/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/plots/initial_mesenchymal_condensation_xz.png', renderView1, ImageResolution=[2458, 1224],
    TransparentBackground=1)

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# set active source
SetActiveSource(threshold1)

# show data in view
threshold1Display = Show(threshold1, renderView1)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(glyph1, renderView1)

# set active source
SetActiveSource(glyph1)

# show data in view
glyph1Display = Show(glyph1, renderView1)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(threshold1, renderView1)

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# hide data in view
Hide(glyph1, renderView1)

# set active source
SetActiveSource(threshold1)

# show data in view
threshold1Display = Show(threshold1, renderView1)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(threshold1, renderView1)

# set active source
SetActiveSource(threshold1)

# set active source
SetActiveSource(glyph1)

# show data in view
glyph1Display = Show(glyph1, renderView1)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# reset view to fit data
renderView1.ResetCamera()

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# set active source
SetActiveSource(threshold1)

# show data in view
threshold1Display = Show(threshold1, renderView1)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(glyph1, renderView1)

# set active source
SetActiveSource(glyph1)

# show data in view
glyph1Display = Show(glyph1, renderView1)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(threshold1, renderView1)

# reset view to fit data
renderView1.ResetCamera()

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# current camera placement for renderView1
renderView1.CameraPosition = [4.331651791321505, 29.712948833992805, 20.12109686071556]
renderView1.CameraFocalPoint = [4.011303067207336, 5.484820067882539, 1.8251452594995503]
renderView1.CameraViewUp = [-0.0015757472035050888, -0.6026147055443349, 0.7980307222672979]
renderView1.CameraParallelScale = 7.858245297373936

# save screenshot
SaveScreenshot('/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/plots/initial_mesenchymal_condensation_view2.png', renderView1, ImageResolution=[2458, 1224],
    TransparentBackground=1)

# hide data in view
Hide(glyph1, renderView1)

# set active source
SetActiveSource(threshold1)

# show data in view
threshold1Display = Show(threshold1, renderView1)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# hide color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, False)

# current camera placement for renderView1
renderView1.CameraPosition = [4.331651791321505, 29.712948833992805, 20.12109686071556]
renderView1.CameraFocalPoint = [4.011303067207336, 5.484820067882539, 1.8251452594995503]
renderView1.CameraViewUp = [-0.0015757472035050888, -0.6026147055443349, 0.7980307222672979]
renderView1.CameraParallelScale = 7.858245297373936

# save screenshot
SaveScreenshot('/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/plots/initial_mesenchymal_condensation_view2_threshold.png', renderView1, ImageResolution=[2458, 1224],
    TransparentBackground=1)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [4.331651791321505, 29.712948833992805, 20.12109686071556]
renderView1.CameraFocalPoint = [4.011303067207336, 5.484820067882539, 1.8251452594995503]
renderView1.CameraViewUp = [-0.0015757472035050888, -0.6026147055443349, 0.7980307222672979]
renderView1.CameraParallelScale = 7.858245297373936

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).