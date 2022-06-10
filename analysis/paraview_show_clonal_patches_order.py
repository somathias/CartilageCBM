import sys, getopt

#### import the simple module from the paraview
from paraview.simple import *


def screenshot(argv):
    
    path = '/home/kubuntu1804/Documents/sf_simulation_results/exp-draft/mesenchymal_condensation/order/20220301-125119/7/'
    time_step = 450
    start='0'
    upper_boundary=3.5
    
    try:
        opts, args = getopt.getopt(argv,"hi:t:s:u:",["ifile=", "time=", "start=", "upper_boundary="])
    except getopt.GetoptError:
        print('No path to input provided via -i flag. Using default.')
    for opt, arg in opts:
        if opt == '-h':
            print('paraview_show_clonal_patches.py -i <path_to_inputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            path = arg
        elif opt in ("-t", "--time"):
            time_step = int(arg)*10
        elif opt in ("-s", "--start"):
            start = arg
        elif opt in ("-u", "--upper_boundary"):
            upper_boundary = float(arg)
            
    print('Input file path is '+ path)
    print('Time step is '+ str(time_step))
    print('Start time is '+ start)
    print('Upper boundary will be placed at '+str(upper_boundary))


    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'PVD Reader'
    resultspvd = PVDReader(FileName=path+'results_from_time_'+start+'/results.pvd')

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
    resultspvdDisplay.ScalarOpacityUnitDistance = 15.098001312075798

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

    # Properties modified on ancestorsLUT
    ancestorsLUT.NumberOfTableValues = 258

    # hide color bar/color legend
    glyph1Display.SetScalarBarVisibility(renderView1, False)

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
    threshold1Display.ScalarOpacityUnitDistance = 0.7885756289502218

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

    # hide color bar/color legend
    threshold1Display.SetScalarBarVisibility(renderView1, False)

    # Properties modified on renderView1.AxesGrid
    renderView1.AxesGrid.ZTitle = ' z [d] '

    # Properties modified on renderView1.AxesGrid
    renderView1.AxesGrid.XAxisUseCustomLabels = 1
    renderView1.AxesGrid.XAxisLabels = [-1.0, 0.0, 1.0, 2.0, 3.0, 6.0, 7.0, 8.0, 9.0]

    # Properties modified on renderView1.AxesGrid
    renderView1.AxesGrid.XAxisLabels = [-1.0, 0.0, 1.0, 2.0, 6.0, 7.0, 8.0, 9.0]

    # create a new 'Plane'
    plane1 = Plane()

    # Properties modified on plane1
    plane1.Origin = [-1.5, -0.5, 0.0]
    plane1.Point1 = [9.5, -0.5, 0.0]
    plane1.Point2 = [-1.5, 6.5, 0.0]
    plane1.XResolution = 9
    plane1.YResolution = 6

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

    # reset view to fit data
    renderView1.ResetCamera()

    # reset view to fit data
    renderView1.ResetCamera()

    # hide data in view
    Hide(threshold1, renderView1)

    # set active source
    SetActiveSource(glyph1)

    # show data in view
    glyph1Display = Show(glyph1, renderView1)

    # show color bar/color legend
    glyph1Display.SetScalarBarVisibility(renderView1, True)

    # reset view to fit data
    renderView1.ResetCamera()

    animationScene1.GoToFirst()

    # reset view to fit data
    renderView1.ResetCamera()

    # reset view to fit data
    renderView1.ResetCamera()

    # reset view to fit data
    renderView1.ResetCamera()

    # set active source
    SetActiveSource(plane1)

    # change representation type
    plane1Display.SetRepresentationType('Surface With Edges')

    # Properties modified on plane1Display
    plane1Display.Opacity = 0.25

    # Properties modified on plane1Display
    plane1Display.EdgeColor = [0.0, 0.0, 0.0]

    # create a new 'Plane'
    plane2 = Plane()

    # Properties modified on plane2
    plane2.Origin = [-1.5, -0.5, upper_boundary]
    plane2.Point1 = [9.5, -0.5, upper_boundary]
    plane2.Point2 = [-1.5, 6.5, upper_boundary]
    plane2.XResolution = 9
    plane2.YResolution = 6

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

    #animationScene1.GoToLast()

    # reset view to fit data
    renderView1.ResetCamera()

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(threshold1)

    # hide color bar/color legend
    threshold1Display.SetScalarBarVisibility(renderView1, False)

    # current camera placement for renderView1
    renderView1.CameraPosition = [4.073004245758057, 34.037772034163694, 3.75]
    renderView1.CameraFocalPoint = [4.073004245758057, 2.628475606441498, 3.75]
    renderView1.CameraViewUp = [0.0, 0.0, 1.0]
    renderView1.CameraParallelScale = 8.129324108765074
    
    tk = GetTimeKeeper()
    timesteps = tk.TimestepValues
    animationScene1.AnimationTime = timesteps[time_step]

    # save screenshot
    SaveScreenshot(path+'threshold_48_minus_y_t'+str(time_step)+'.png', renderView1, ImageResolution=[2458, 1464],
        TransparentBackground=1)
    
    animationScene1.AnimationTime = timesteps[0]
    # save screenshot
    SaveScreenshot(path+'threshold_48_minus_y_t'+str(0)+'.png', renderView1, ImageResolution=[2458, 1464], TransparentBackground=1)
  

    # set active source
    SetActiveSource(resultspvd)

    #### saving camera placements for all active views

    # current camera placement for renderView1
    renderView1.CameraPosition = [4.073004245758057, 34.037772034163694, 3.75]
    renderView1.CameraFocalPoint = [4.073004245758057, 2.628475606441498, 3.75]
    renderView1.CameraViewUp = [0.0, 0.0, 1.0]
    renderView1.CameraParallelScale = 8.129324108765074

    #### uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).
    
    
if __name__ == "__main__":
   screenshot(sys.argv[1:])    
