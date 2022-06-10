#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
resultspvd = PVDReader(FileName='/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/baseline/20220609-141147/1/results_from_time_0/results.pvd')

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1229, 612]

# Properties modified on renderView1
renderView1.Background = [1.0, 1.0, 1.0]

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.Visibility = 1

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XTitle = 'x [d]'
renderView1.AxesGrid.YTitle = 'y [d]'
renderView1.AxesGrid.ZTitle = 'y [d]'
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
resultspvdDisplay.ScalarOpacityUnitDistance = 14.226415837691832

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

# reset view to fit data
renderView1.ResetCamera()

# set scalar coloring
ColorBy(glyph1Display, ('POINTS', 'Ancestors'))

# rescale color and/or opacity maps used to include current data range
glyph1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Ancestors'
ancestorsLUT = GetColorTransferFunction('Ancestors')

# Properties modified on ancestorsLUT
ancestorsLUT.NumberOfTableValues = 238

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# reset view to fit data
renderView1.ResetCamera()

# set scalar coloring
ColorBy(glyph1Display, ('POINTS', 'Cell tissue types'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(ancestorsLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
glyph1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Celltissuetypes'
celltissuetypesLUT = GetColorTransferFunction('Celltissuetypes')

# create a new 'Plane'
plane1 = Plane()

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

# Properties modified on plane1
plane1.Point1 = [8.5, -0.5, 0.0]
plane1.Point2 = [-0.5, 10.5, 0.0]
plane1.XResolution = 11
plane1.YResolution = 9

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on plane1
plane1.Origin = [-1.0, -1.0, 0.0]
plane1.Point1 = [9.0, -1.0, 0.0]
plane1.Point2 = [-1.0, 11.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# change representation type
plane1Display.SetRepresentationType('Surface With Edges')

# Properties modified on plane1Display
plane1Display.Opacity = 0.5

# change solid color
plane1Display.DiffuseColor = [0.0, 0.0, 0.0]

# Properties modified on plane1Display
plane1Display.Opacity = 0.25

# Properties modified on plane1Display
plane1Display.EdgeColor = [0.0, 0.0, 0.0]

# change solid color
plane1Display.DiffuseColor = [1.0, 1.0, 1.0]

# create a new 'Plane'
plane2 = Plane()

# Properties modified on plane2
plane2.Origin = [-1.0, -1.0, 5.5]
plane2.Point1 = [9.0, -1.0, 5.5]
plane2.Point2 = [-1.0, 11.0, 5.5]

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

# Properties modified on plane2
plane2.Origin = [-1.0, -1.0, 4.5]
plane2.Point1 = [9.0, -1.0, 4.5]
plane2.Point2 = [-1.0, 11.0, 4.5]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on plane2
plane2.Origin = [-1.0, -1.0, 5.0]
plane2.Point1 = [9.0, -1.0, 5.0]
plane2.Point2 = [-1.0, 11.0, 5.0]

# update the view to ensure updated data information
renderView1.Update()

# change representation type
plane2Display.SetRepresentationType('Surface With Edges')

# Properties modified on plane2Display
plane2Display.Opacity = 0.25

# Properties modified on plane2Display
plane2Display.EdgeColor = [0.0, 0.0, 0.0]

# Properties modified on plane2
plane2.Origin = [-1.0, -1.0, 5.5]
plane2.Point1 = [9.0, -1.0, 5.5]
plane2.Point2 = [-1.0, 11.0, 5.5]

# update the view to ensure updated data information
renderView1.Update()

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

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.YAxisUseCustomLabels = 1
renderView1.AxesGrid.YAxisLabels = [0.0, 2.0, 8.0, 10.0]

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.ZTitle = 'z [d]'

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.ZTitle = ' z [d] '

# set active source
SetActiveSource(plane1)

# set active source
SetActiveSource(glyph1)

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# current camera placement for renderView1
renderView1.CameraPosition = [-28.44413149940346, 5.0, 2.507368355989456]
renderView1.CameraFocalPoint = [4.0, 5.0, 2.507368355989456]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 8.397159133856219

# save screenshot
SaveScreenshot('/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/baseline/20220609-141147/1/cartilage_sheet_configuration_+x_cell_tissue_types.png', renderView1, ImageResolution=[2458, 1224],
    TransparentBackground=1)

# set scalar coloring
ColorBy(glyph1Display, ('POINTS', 'Ancestors'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(celltissuetypesLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
glyph1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

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

# set scalar coloring
ColorBy(glyph1Display, ('POINTS', 'Cell division directions'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(ancestorsLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
glyph1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Celldivisiondirections'
celldivisiondirectionsLUT = GetColorTransferFunction('Celldivisiondirections')

# set scalar coloring
ColorBy(glyph1Display, ('POINTS', 'Cell tissue types'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(celldivisiondirectionsLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
glyph1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XAxisUseCustomLabels = 1
renderView1.AxesGrid.XAxisLabels = [0.0, 2.0, 6.0, 8.0]

# reset view to fit data
renderView1.ResetCamera()

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# current camera placement for renderView1
renderView1.CameraPosition = [4.0, 37.44413149940346, 2.507368355989456]
renderView1.CameraFocalPoint = [4.0, 5.0, 2.507368355989456]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 8.397159133856219

# save screenshot
SaveScreenshot('/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/baseline/20220609-141147/1/cartilage_sheet_configuration_-y_cell_tissue_types.png', renderView1, ImageResolution=[2458, 1224],
    TransparentBackground=1)

# set scalar coloring
ColorBy(glyph1Display, ('POINTS', 'Ancestors'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(celltissuetypesLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
glyph1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# create a new 'Threshold'
threshold1 = Threshold(Input=glyph1)
threshold1.Scalars = ['POINTS', 'Ancestors']
threshold1.ThresholdRange = [-1.0, 558.0]

# Properties modified on threshold1
threshold1.ThresholdRange = [-1.0, 0.0]

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
threshold1Display.ScalarOpacityUnitDistance = 0.4188118875381537

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
threshold1Display.OSPRayScaleFunction.Points = [-1.0, 0.0, 0.5, 0.0, 147.0, 1.0, 0.5, 0.0]

# hide data in view
Hide(glyph1, renderView1)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on threshold1Display
threshold1Display.Opacity = 0.2

# Properties modified on threshold1Display
threshold1Display.Opacity = 0.3

# create a new 'Threshold'
threshold2 = Threshold(Input=threshold1)
threshold2.Scalars = ['POINTS', 'Ancestors']
threshold2.ThresholdRange = [-1.0, -1.0]

# set active source
SetActiveSource(threshold1)

# destroy threshold2
Delete(threshold2)
del threshold2

# set active source
SetActiveSource(glyph1)

# create a new 'Threshold'
threshold2 = Threshold(Input=glyph1)
threshold2.Scalars = ['POINTS', 'Ancestors']
threshold2.ThresholdRange = [-1.0, 558.0]

# Properties modified on threshold2
threshold2.ThresholdRange = [0.0, 558.0]

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
threshold2Display.ScalarOpacityUnitDistance = 1.2257701524375688

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
threshold2Display.OSPRayScaleFunction.Points = [-1.0, 0.0, 0.5, 0.0, 147.0, 1.0, 0.5, 0.0]

# hide data in view
Hide(glyph1, renderView1)

# show color bar/color legend
threshold2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

animationScene1.GoToLast()

# set active source
SetActiveSource(threshold1)

# Properties modified on threshold1Display
threshold1Display.Opacity = 0.1

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

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# Properties modified on threshold1Display
threshold1Display.Opacity = 0.05

# hide color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, False)

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XAxisLabels = [0.0, 6.0, 8.0]

# current camera placement for renderView1
renderView1.CameraPosition = [29.299012578491467, -14.634416126818712, 10.920274591442785]
renderView1.CameraFocalPoint = [4.078896045684811, 4.962783098220823, 2.7500000000000004]
renderView1.CameraViewUp = [-0.20795676204059807, 0.13619854567421769, 0.9686092820522721]
renderView1.CameraParallelScale = 8.532624037305347

# save screenshot
SaveScreenshot('/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/baseline/20220609-141147/1/cartilage_sheet_columns_t80.png', renderView1, ImageResolution=[2458, 1224],
    TransparentBackground=1)

animationScene1.GoToFirst()

# current camera placement for renderView1
renderView1.CameraPosition = [29.299012578491467, -14.634416126818712, 10.920274591442785]
renderView1.CameraFocalPoint = [4.078896045684811, 4.962783098220823, 2.7500000000000004]
renderView1.CameraViewUp = [-0.20795676204059807, 0.13619854567421769, 0.9686092820522721]
renderView1.CameraParallelScale = 8.532624037305347

# save screenshot
SaveScreenshot('/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/baseline/20220609-141147/1/cartilage_sheet_columns_t0.png', renderView1, ImageResolution=[2458, 1224],
    TransparentBackground=1)

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

# current camera placement for renderView1
renderView1.CameraPosition = [29.299012578491467, -14.634416126818712, 10.920274591442785]
renderView1.CameraFocalPoint = [4.078896045684811, 4.962783098220823, 2.7500000000000004]
renderView1.CameraViewUp = [-0.20795676204059807, 0.13619854567421769, 0.9686092820522721]
renderView1.CameraParallelScale = 8.532624037305347

# save screenshot
SaveScreenshot('/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/baseline/20220609-141147/1/cartilage_sheet_t..png', renderView1, ImageResolution=[2458, 1224],
    TransparentBackground=1)

animationScene1.GoToLast()

# current camera placement for renderView1
renderView1.CameraPosition = [29.299012578491467, -14.634416126818712, 10.920274591442785]
renderView1.CameraFocalPoint = [4.078896045684811, 4.962783098220823, 2.7500000000000004]
renderView1.CameraViewUp = [-0.20795676204059807, 0.13619854567421769, 0.9686092820522721]
renderView1.CameraParallelScale = 8.532624037305347

# save screenshot
SaveScreenshot('/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/baseline/20220609-141147/1/cartilage_sheet_t80.png', renderView1, ImageResolution=[2458, 1224],
    TransparentBackground=1)

# reset view to fit data
renderView1.ResetCamera()

animationScene1.GoToFirst()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# hide data in view
Hide(plane2, renderView1)

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# current camera placement for renderView1
renderView1.CameraPosition = [4.0, 5.0, 34.95149985539292]
renderView1.CameraFocalPoint = [4.0, 5.0, 2.507368355989456]
renderView1.CameraParallelScale = 8.397159133856219

# save screenshot
SaveScreenshot('/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/cartilage_sheet/baseline/20220609-141147/1/cartilage_sheet_-z_t0.png', renderView1, ImageResolution=[2458, 1224],
    TransparentBackground=1)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [4.0, 5.0, 34.95149985539292]
renderView1.CameraFocalPoint = [4.0, 5.0, 2.507368355989456]
renderView1.CameraParallelScale = 8.397159133856219

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).