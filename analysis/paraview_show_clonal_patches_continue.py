#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
resultspvd = PVDReader(FileName='/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/mesenchymal_condensation/continue/20220301-143843/0/results_from_time_50/results.pvd')

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
resultspvdDisplay.ScalarOpacityUnitDistance = 14.318492701569367

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

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# Properties modified on ancestorsLUT
ancestorsLUT.NumberOfTableValues = 258

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
threshold1Display.ScalarOpacityUnitDistance = 0.46745159917040213

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
threshold1Display.OSPRayScaleFunction.Points = [-1.0, 0.0, 0.5, 0.0, 147.0, 1.0, 0.5, 0.0]

# hide data in view
Hide(glyph1, renderView1)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# reset view to fit data
renderView1.ResetCamera()

animationScene1.GoToLast()

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

# set active source
SetActiveSource(resultspvd)

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XAxisUseCustomLabels = 1
renderView1.AxesGrid.XAxisLabels = [0.0, 2.0, 6.0, 8.0]

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.ZTitle = '  z [d]  '

# set active source
SetActiveSource(threshold1)

# hide color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, False)

# current camera placement for renderView1
renderView1.CameraPosition = [4.223685383796692, 29.434890919771497, 1.75]
renderView1.CameraFocalPoint = [4.223685383796692, 2.146424949169159, 1.75]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 7.062774704823929

# save screenshot
SaveScreenshot('/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/mesenchymal_condensation/continue/20220301-143843/0/threshold_t800.png', renderView1, ImageResolution=[2458, 1224],
    TransparentBackground=1)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [4.223685383796692, 29.434890919771497, 1.75]
renderView1.CameraFocalPoint = [4.223685383796692, 2.146424949169159, 1.75]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 7.062774704823929

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).