#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
resultspvd = PVDReader(FileName='/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/0/results_from_time_0/results.pvd')

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
resultspvdDisplay.ScalarOpacityUnitDistance = 7.527427615808969

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
resultspvdDisplay.OSPRayScaleFunction.Points = [-1.0, 0.0, 0.5, 0.0, 147.0, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on renderView1
renderView1.Background = [1.0, 1.0, 1.0]

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.Visibility = 1

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XTitle = 'x [d]'
renderView1.AxesGrid.YTitle = 'y [d]'
renderView1.AxesGrid.ZTitle = 'z [d]'
renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XTitleFontSize = 46
renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YTitleFontSize = 46
renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZTitleFontSize = 46
renderView1.AxesGrid.GridColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelFontSize = 42
renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YLabelFontSize = 42
renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZLabelFontSize = 42

# create a new 'Glyph'
glyph1 = Glyph(Input=resultspvd,
    GlyphType='Arrow')
glyph1.Scalars = ['POINTS', 'Ancestors']
glyph1.Vectors = ['POINTS', 'None']
glyph1.ScaleFactor = 0.0
glyph1.GlyphTransform = 'Transform2'

# Properties modified on glyph1
glyph1.GlyphType = 'Sphere'
glyph1.Scalars = ['POINTS', 'None']
glyph1.ScaleMode = 'scalar'
glyph1.ScaleFactor = 1.0

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

# Properties modified on glyph1
glyph1.GlyphMode = 'All Points'

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

# reset view to fit data
renderView1.ResetCamera()

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

animationScene1.GoToLast()

animationScene1.GoToFirst()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

animationScene1.Play()

animationScene1.GoToLast()

animationScene1.GoToFirst()

# set active source
SetActiveSource(resultspvd)

animationScene1.GoToLast()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/1/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/2/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/3/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(glyph1)

# create a new 'Threshold'
threshold1 = Threshold(Input=glyph1)
threshold1.Scalars = ['POINTS', 'Ancestors']
threshold1.ThresholdRange = [-1.0, 6.0]

# set active source
SetActiveSource(glyph1)

# create a new 'Threshold'
threshold2 = Threshold(Input=glyph1)
threshold2.Scalars = ['POINTS', 'Ancestors']
threshold2.ThresholdRange = [-1.0, 6.0]

# set active source
SetActiveSource(threshold1)

# Properties modified on threshold1
threshold1.ThresholdRange = [0.0, 6.0]

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
threshold1Display.ScalarOpacityUnitDistance = 0.5572078110901317

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
threshold1Display.OSPRayScaleFunction.Points = [-1.0, 0.0, 0.5, 0.0, 147.0, 1.0, 0.5, 0.0]

# hide data in view
Hide(glyph1, renderView1)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# show data in view
threshold2Display = Show(threshold2, renderView1)
# trace defaults for the display properties.
threshold2Display.Representation = 'Surface'
threshold2Display.ColorArrayName = ['POINTS', 'Ancestors']
threshold2Display.LookupTable = ancestorsLUT
threshold2Display.OSPRayScaleArray = 'Ancestors'
threshold2Display.OSPRayScaleFunction = 'PiecewiseFunction'
threshold2Display.SelectOrientationVectors = 'Ancestors'
threshold2Display.ScaleFactor = 0.0
threshold2Display.SelectScaleArray = 'Ancestors'
threshold2Display.GlyphType = 'Arrow'
threshold2Display.GlyphTableIndexArray = 'Ancestors'
threshold2Display.DataAxesGrid = 'GridAxesRepresentation'
threshold2Display.PolarAxes = 'PolarAxesRepresentation'
threshold2Display.ScalarOpacityFunction = ancestorsPWF
threshold2Display.ScalarOpacityUnitDistance = 0.4018409767529291

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
threshold2Display.OSPRayScaleFunction.Points = [-1.0, 0.0, 0.5, 0.0, 147.0, 1.0, 0.5, 0.0]

# hide data in view
Hide(glyph1, renderView1)

# show color bar/color legend
threshold2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# hide data in view
Hide(threshold1, renderView1)

# set active source
SetActiveSource(threshold1)

# show data in view
threshold1Display = Show(threshold1, renderView1)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(threshold2)

# hide data in view
Hide(threshold2, renderView1)

# reset view to fit data
renderView1.ResetCamera()

animationScene1.GoToLast()

# set active source
SetActiveSource(resultspvd)

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/2/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/4/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/5/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/6/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/7/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/8/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/9/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/10/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/11/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/12/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/13/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/14/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/15/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/16/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/17/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/18/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/19/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on resultspvd
resultspvd.FileName = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/2/results_from_time_0/results.pvd'

# update the view to ensure updated data information
renderView1.Update()

animationScene1.GoToFirst()

animationScene1.GoToLast()

# reset view to fit data
renderView1.ResetCamera()

# set active source
SetActiveSource(threshold2)

# Properties modified on threshold2Display
threshold2Display.Opacity = 0.1

# Properties modified on threshold2
threshold2.ThresholdRange = [-1.0, 0.0]

# set active source
SetActiveSource(threshold2)

# show data in view
threshold2Display = Show(threshold2, renderView1)

# show color bar/color legend
threshold2Display.SetScalarBarVisibility(renderView1, True)

# reset view to fit data
renderView1.ResetCamera()

# hide color bar/color legend
threshold2Display.SetScalarBarVisibility(renderView1, False)

# hide data in view
Hide(threshold1, renderView1)

# hide data in view
Hide(threshold2, renderView1)

# set active source
SetActiveSource(glyph1)

# show data in view
glyph1Display = Show(glyph1, renderView1)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

animationScene1.GoToFirst()

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XAxisUseCustomLabels = 1
renderView1.AxesGrid.XAxisLabels = [0.0, 1.0, 4.0, 5.0]

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.ZTitle = '  z [d]  '

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.ZAxisUseCustomLabels = 1
renderView1.AxesGrid.ZAxisLabels = [0.0, 1.0, 2.0, 3.0, 4.0]

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# current camera placement for renderView1
renderView1.CameraPosition = [2.426445722579956, 19.45732457406211, 2.0]
renderView1.CameraFocalPoint = [2.426445722579956, 1.8432855606079102, 2.0]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 4.558848757860763

# save screenshot
SaveScreenshot('/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/2/whole_sheet_-y_t0.png', renderView1, ImageResolution=[2458, 1224],
    TransparentBackground=1)

# set active source
SetActiveSource(resultspvd)

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.YTitle = '  y [d]  '
renderView1.AxesGrid.YAxisUseCustomLabels = 1

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.YAxisUseCustomLabels = 0
renderView1.AxesGrid.YAxisLabels = [0.0]

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.YAxisUseCustomLabels = 1
renderView1.AxesGrid.YAxisLabels = [0.0, 1.0, 2.0, 3.0, 4.0]

# current camera placement for renderView1
renderView1.CameraPosition = [2.4273074567317963, 1.8387418687343597, -15.690503561541114]
renderView1.CameraFocalPoint = [2.4273074567317963, 1.8387418687343597, 1.9992891550064087]
renderView1.CameraParallelScale = 4.578455258958356

# save screenshot
SaveScreenshot('/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/2/whole_sheet_+z_t0.png', renderView1, ImageResolution=[2458, 1224],
    TransparentBackground=1)

# reset view to fit data
renderView1.ResetCamera()

animationScene1.GoToLast()

# current camera placement for renderView1
renderView1.CameraPosition = [2.4273074567317963, 19.528534585281882, 1.9992891550064087]
renderView1.CameraFocalPoint = [2.4273074567317963, 1.8387418687343597, 1.9992891550064087]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 4.578455258958356

# save screenshot
SaveScreenshot('/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/2/whole_sheet_-y_t80.png', renderView1, ImageResolution=[2458, 1224],
    TransparentBackground=1)

# hide data in view
Hide(glyph1, renderView1)

# set active source
SetActiveSource(threshold1)

# show data in view
threshold1Display = Show(threshold1, renderView1)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(threshold2)

# show data in view
threshold2Display = Show(threshold2, renderView1)

# show color bar/color legend
threshold2Display.SetScalarBarVisibility(renderView1, True)

# hide color bar/color legend
threshold2Display.SetScalarBarVisibility(renderView1, False)

# current camera placement for renderView1
renderView1.CameraPosition = [2.4273074567317963, 19.528534585281882, 1.9992891550064087]
renderView1.CameraFocalPoint = [2.4273074567317963, 1.8387418687343597, 1.9992891550064087]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 4.578455258958356

# save screenshot
SaveScreenshot('/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/missing_column/20220612-151424/2/column_only_-y_t80.png', renderView1, ImageResolution=[2458, 1224],
    TransparentBackground=1)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [2.4273074567317963, 19.528534585281882, 1.9992891550064087]
renderView1.CameraFocalPoint = [2.4273074567317963, 1.8387418687343597, 1.9992891550064087]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 4.578455258958356

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).