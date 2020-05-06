import sys, getopt
#### import the simple module from the paraview
from paraview.simple import *


def screenshot(argv): 
       
    path = '/home/kubuntu1804/Documents/sf_simulation_results/exp-optimal_adhesion/20190306-151409/0.5/0.5/0/'
    time_step = 450
    
    try:
        opts, args = getopt.getopt(argv,"hi:t:",["ifile=", "time="])
    except getopt.GetoptError:
        print('No path to input provided via -i flag. Using default.')
    for opt, arg in opts:
        if opt == '-h':
            print('paraview_show_clonal_patches.py -i <path_to_inputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            path = arg
        elif opt in ("-t", "--time"):
            time_step = int(arg)
            
    print('Input file path is '+ path)
    print('Time step is '+ str(time_step))


    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()
    
    # create a new 'PVD Reader'
    resultspvd = PVDReader(FileName=path+'results_from_time_0/results.pvd')
    #resultspvd = PVDReader(FileName='/home/kubuntu1804/Documents/sf_simulation_results/exp-random_direction/20190219-155036/0.0/0.5/0/results_from_time_0/results.pvd')
    
    # get animation scene
    animationScene1 = GetAnimationScene()
    
    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()
    
    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    renderView1.ViewSize = [1229, 732]
    
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
    resultspvdDisplay.ScalarOpacityUnitDistance = 12.72729312787153
    
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
    renderView1.UseGradientBackground = 1
    
    # Properties modified on renderView1
    renderView1.Background = [0.0, 0.0, 0.0]
    
    # Properties modified on renderView1
    renderView1.Background = [1.0, 1.0, 1.0]
    
    # Properties modified on renderView1
    renderView1.Background2 = [0.3058823529411765, 0.3058823529411765, 0.3058823529411765]
    
    # Properties modified on renderView1.AxesGrid
    renderView1.AxesGrid.Visibility = 1
    
    # Properties modified on renderView1.AxesGrid
    renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.GridColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
    
    # Properties modified on glyph1
    glyph1.GlyphType = 'Sphere'
    glyph1.ScaleFactor = 1.0
    glyph1.GlyphMode = 'All Points'
    
    # get color transfer function/color map for 'GlyphScale'
    glyphScaleLUT = GetColorTransferFunction('GlyphScale')
    glyphScaleLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 47.5, 0.865003, 0.865003, 0.865003, 95.0, 0.705882, 0.0156863, 0.14902]
    glyphScaleLUT.ScalarRangeInitialized = 1.0
    
    # show data in view
    glyph1Display = Show(glyph1, renderView1)
    # trace defaults for the display properties.
    glyph1Display.Representation = 'Surface'
    glyph1Display.ColorArrayName = ['POINTS', 'GlyphScale']
    glyph1Display.LookupTable = glyphScaleLUT
    glyph1Display.OSPRayScaleArray = 'GlyphScale'
    glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    glyph1Display.SelectOrientationVectors = 'Ancestors'
    glyph1Display.ScaleFactor = 0.0
    glyph1Display.SelectScaleArray = 'GlyphScale'
    glyph1Display.GlyphType = 'Arrow'
    glyph1Display.GlyphTableIndexArray = 'GlyphScale'
    glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
    glyph1Display.PolarAxes = 'PolarAxesRepresentation'
    
    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    glyph1Display.OSPRayScaleFunction.Points = [-1.0, 0.0, 0.5, 0.0, 147.0, 1.0, 0.5, 0.0]
    
    # show color bar/color legend
    glyph1Display.SetScalarBarVisibility(renderView1, True)
    
    # update the view to ensure updated data information
    renderView1.Update()
    
    # reset view to fit data
    renderView1.ResetCamera()
    
    # hide color bar/color legend
    glyph1Display.SetScalarBarVisibility(renderView1, False)
    
    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    glyphScaleLUT.ApplyPreset('periodic', True)
    
    animationScene1.GoToLast()
    
    # reset view to fit data
    renderView1.ResetCamera()
    
    # Properties modified on glyphScaleLUT
    glyphScaleLUT.NumberOfTableValues = 283
    
    # get opacity transfer function/opacity map for 'GlyphScale'
    glyphScalePWF = GetOpacityTransferFunction('GlyphScale')
    glyphScalePWF.Points = [0.0, 0.0, 0.5, 0.0, 95.0, 1.0, 0.5, 0.0]
    glyphScalePWF.ScalarRangeInitialized = 1
    
    # reset view to fit data
    renderView1.ResetCamera()
    
    #animationScene1.GoToLast()
    
    # current camera placement for renderView1
    renderView1.CameraPosition = [-19.125742985611605, 4.73473185300827, 1.7598326057195663]
    renderView1.CameraFocalPoint = [3.751069515943527, 4.73473185300827, 1.7598326057195663]
    renderView1.CameraViewUp = [0.0, 0.0, 1.0]
    renderView1.CameraParallelScale = 7.16435526763671
    
    tk = GetTimeKeeper()
    timesteps = tk.TimestepValues
    animationScene1.AnimationTime = timesteps[time_step]
    
    # save screenshot
    SaveScreenshot(path+'plus_x_t'+str(time_step)+'.png', renderView1, ImageResolution=[1229, 732])
    
    # create a new 'Threshold'
    threshold1 = Threshold(Input=glyph1)
    threshold1.Scalars = ['POINTS', 'GlyphScale']
    threshold1.ThresholdRange = [0.0, 95.0]
    
    # set active source
    #SetActiveSource(glyph1)
    
    
    # reset view to fit data
    renderView1.ResetCamera()
    
    
    # get opacity transfer function/opacity map for 'GlyphScale'
    #glyphScalePWF = GetOpacityTransferFunction('GlyphScale')
    #glyphScalePWF.Points = [0.0, 0.0, 0.5, 0.0, 95.0, 1.0, 0.5, 0.0]
    #glyphScalePWF.ScalarRangeInitialized = 1
    
    # show data in view
    threshold1Display = Show(threshold1, renderView1)
    # trace defaults for the display properties.
    threshold1Display.Representation = 'Surface'
    threshold1Display.ColorArrayName = ['POINTS', 'GlyphScale']
    threshold1Display.LookupTable = glyphScaleLUT
    threshold1Display.OSPRayScaleArray = 'GlyphScale'
    threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    threshold1Display.SelectOrientationVectors = 'Ancestors'
    threshold1Display.ScaleFactor = 0.0
    threshold1Display.SelectScaleArray = 'GlyphScale'
    threshold1Display.GlyphType = 'Arrow'
    threshold1Display.GlyphTableIndexArray = 'GlyphScale'
    threshold1Display.DataAxesGrid = 'GridAxesRepresentation'
    threshold1Display.PolarAxes = 'PolarAxesRepresentation'
    threshold1Display.ScalarOpacityFunction = glyphScalePWF
    threshold1Display.ScalarOpacityUnitDistance = 0.48851975457835795
    
    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    threshold1Display.OSPRayScaleFunction.Points = [-1.0, 0.0, 0.5, 0.0, 147.0, 1.0, 0.5, 0.0]
    
    # hide data in view
    Hide(glyph1, renderView1)
    
    
    # update the view to ensure updated data information
    renderView1.Update()
    
    # set active source
    SetActiveSource(threshold1)
    
    # reset view to fit data
    renderView1.ResetCamera()
    
    # Properties modified on threshold1
    threshold1.ThresholdRange = [0.0, 48.0]
    
    # update the view to ensure updated data information
    renderView1.Update()
    
    # reset view to fit data
    renderView1.ResetCamera()
    
    # show data in view
    threshold1Display = Show(threshold1, renderView1)
    
    # hide data in view
    Hide(glyph1, renderView1)
    
    # update the view to ensure updated data information
    renderView1.Update()
    
    # reset view to fit data
    renderView1.ResetCamera()
    
    #animationScene1.GoToLast()
    animationScene1.AnimationTime = timesteps[time_step]
    
    # current camera placement for renderView1
    renderView1.CameraPosition = [3.854806661605835, 19.1787146788673, 1.75]
    renderView1.CameraFocalPoint = [3.854806661605835, 2.17948997020721, 1.75]
    renderView1.CameraViewUp = [0.0, 0.0, 1.0]
    renderView1.CameraParallelScale = 6.441634600341691
    
    # save screenshot
    SaveScreenshot(path+'threshold_48_minus_y_t'+str(time_step)+'.png', renderView1, ImageResolution=[1229, 732])
     
    #### uncomment the following to render all views
    #RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).
    
if __name__ == "__main__":
   screenshot(sys.argv[1:])    
