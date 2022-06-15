import sys, getopt
#### import the simple module from the paraview
from paraview.simple import *


def screenshot(argv): 
       
    path = '/home/kubuntu1804/Documents/sf_simulation_results/exp-optimal_adhesion/20190306-151409/0.5/0.5/0/'
    time_step = 450
    start='0'
    
    try:
        opts, args = getopt.getopt(argv,"hi:t:s:",["ifile=", "time=", "start="])
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
        elif opt in ("-s", "--start"):
            start = arg
            
    print('Input file path is '+ path)
    print('Time step is '+ str(time_step))
    print('Start time is '+ start)


    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()
    
    # create a new 'PVD Reader'
    resultspvd = PVDReader(FileName=path+'results_from_time_'+start+'/results.pvd')
    #resultspvd = PVDReader(FileName='/home/kubuntu1804/Documents/sf_simulation_results/exp-random_direction/20190219-155036/0.0/0.5/0/results_from_time_0/results.pvd')
    
    # get animation scene
    animationScene1 = GetAnimationScene()
    
    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()
    
    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    renderView1.ViewSize = [1625, 618]
    
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
    resultspvdDisplay.ScalarOpacityUnitDistance = 20.89849674222659
    
    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    resultspvdDisplay.OSPRayScaleFunction.Points = [-1.0, 0.0, 0.5, 0.0, 147.0, 1.0, 0.5, 0.0]
    
    # reset view to fit data
    renderView1.ResetCamera()
    
    # update the view to ensure updated data information
    renderView1.Update()
    
    # create a new 'Glyph'
    glyph1 = Glyph(Input=resultspvd,
        GlyphType='Sphere')
    glyph1.Scalars = ['POINTS', 'None']
    glyph1.Vectors = ['POINTS', 'None']
    glyph1.ScaleFactor = 1.0
    glyph1.GlyphTransform = 'Transform2'
    glyph1.GlyphMode = 'All Points'

    
#    # Properties modified on renderView1
#    renderView1.UseGradientBackground = 1
#    
#    # Properties modified on renderView1
#    renderView1.Background = [0.0, 0.0, 0.0]
#    
#    # Properties modified on renderView1
#    renderView1.Background = [1.0, 1.0, 1.0]
#    
#    # Properties modified on renderView1
#    renderView1.Background2 = [0.3058823529411765, 0.3058823529411765, 0.3058823529411765]
    
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
    
    renderView1.AxesGrid.XTitleFontSize = 80
    renderView1.AxesGrid.YTitleFontSize = 80
    renderView1.AxesGrid.ZTitleFontSize = 80
    renderView1.AxesGrid.XLabelFontSize = 60
    renderView1.AxesGrid.YLabelFontSize = 60
    renderView1.AxesGrid.ZLabelFontSize = 60
    
    renderView1.AxesGrid.XTitle = 'x [d]'
    renderView1.AxesGrid.YTitle = 'y [d]'
    renderView1.AxesGrid.ZTitle = '  z [d]  '
    
    renderView1.AxesGrid.XAxisUseCustomLabels = 1
    renderView1.AxesGrid.XAxisLabels= [0, 2, 6, 8]
    renderView1.AxesGrid.ZAxisUseCustomLabels = 1
    renderView1.AxesGrid.ZAxisLabels= [0,1, 2, 3]
    
        
    # get color transfer function/color map for 'GlyphScale'
#    glyphScaleLUT = GetColorTransferFunction('GlyphScale')
#    glyphScaleLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 47.5, 0.865003, 0.865003, 0.865003, 95.0, 0.705882, 0.0156863, 0.14902]
#    glyphScaleLUT.ScalarRangeInitialized = 1.0
    
    # show data in view
    glyph1Display = Show(glyph1, renderView1)
    # trace defaults for the display properties.
    glyph1Display.Representation = 'Surface'
    glyph1Display.ColorArrayName = ['POINTS', 'Ancestors']
#    glyph1Display.LookupTable = glyphScaleLUT
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
    
    # get color transfer function/color map for 'Ancestors'
    ancestorsLUT = GetColorTransferFunction('Ancestors')
    ancestorsLUT.RGBPoints = [0.0, 0.4961, 0.5078, 0.7539, 0.39879600000000004, 0.929, 0.694, 0.125, 0.7975920000000001, 0.494, 0.184, 0.556, 1.196388, 0.466, 0.674, 0.188, 1.5951840000000002, 0.301, 0.745, 0.933, 1.9939799999999999, 0.635, 0.078, 0.184, 2.392776, 0.0, 0.447, 0.741, 2.791572, 0.85, 0.325, 0.098, 3.1903680000000003, 0.929, 0.694, 0.125, 3.589164, 0.494, 0.184, 0.556, 3.9879599999999997, 0.466, 0.674, 0.188, 4.386756, 0.301, 0.745, 0.933, 4.785552, 0.635, 0.078, 0.184, 5.184347999999999, 0.0, 0.447, 0.741, 5.583144, 0.85, 0.325, 0.098, 5.98194, 0.929, 0.694, 0.125, 6.380736000000001, 0.494, 0.184, 0.556, 6.779532000000001, 0.466, 0.674, 0.188, 7.178328, 0.301, 0.745, 0.933, 7.5771239999999995, 0.635, 0.078, 0.184, 7.9759199999999995, 0.0, 0.447, 0.741, 8.374716000000001, 0.85, 0.325, 0.098, 8.773512, 0.929, 0.694, 0.125, 9.172308000000001, 0.494, 0.184, 0.556, 9.571104, 0.466, 0.674, 0.188, 9.969899999999999, 0.301, 0.745, 0.933, 10.368695999999998, 0.635, 0.078, 0.184, 10.767492, 0.0, 0.447, 0.741, 11.166288, 0.85, 0.325, 0.098, 11.565084, 0.929, 0.694, 0.125, 11.96388, 0.494, 0.184, 0.556, 12.362676, 0.466, 0.674, 0.188, 12.761472000000001, 0.301, 0.745, 0.933, 13.160267999999999, 0.635, 0.078, 0.184, 13.559064000000001, 0.0, 0.447, 0.741, 13.957859999999998, 0.85, 0.325, 0.098, 14.356656, 0.929, 0.694, 0.125, 14.755452000000002, 0.494, 0.184, 0.556, 15.154247999999999, 0.466, 0.674, 0.188, 15.553044, 0.301, 0.745, 0.933, 15.951839999999999, 0.635, 0.078, 0.184, 16.350636, 0.0, 0.447, 0.741, 16.749432000000002, 0.85, 0.325, 0.098, 17.148228, 0.929, 0.694, 0.125, 17.547024, 0.494, 0.184, 0.556, 17.945819999999998, 0.466, 0.674, 0.188, 18.344616000000002, 0.301, 0.745, 0.933, 18.743411999999996, 0.635, 0.078, 0.184, 19.142208, 0.0, 0.447, 0.741, 19.541004, 0.85, 0.325, 0.098, 19.939799999999998, 0.929, 0.694, 0.125, 20.338596, 0.494, 0.184, 0.556, 20.737391999999996, 0.466, 0.674, 0.188, 21.136188, 0.301, 0.745, 0.933, 21.534984, 0.635, 0.078, 0.184, 21.93378, 0.0, 0.447, 0.741, 22.332576, 0.85, 0.325, 0.098, 22.731371999999997, 0.929, 0.694, 0.125, 23.130168, 0.494, 0.184, 0.556, 23.528964, 0.466, 0.674, 0.188, 23.92776, 0.301, 0.745, 0.933, 24.326556, 0.635, 0.078, 0.184, 24.725352, 0.0, 0.447, 0.741, 25.124347, 0.85, 0.325, 0.098, 25.523143, 0.929, 0.694, 0.125, 25.921939, 0.494, 0.184, 0.556, 26.320735, 0.466, 0.674, 0.188, 26.719531, 0.301, 0.745, 0.933, 27.118327, 0.635, 0.078, 0.184, 27.517123, 0.0, 0.447, 0.741, 27.915919, 0.85, 0.325, 0.098, 28.314715, 0.929, 0.694, 0.125, 28.713511, 0.494, 0.184, 0.556, 29.112307, 0.466, 0.674, 0.188, 29.511103000000002, 0.301, 0.745, 0.933, 29.909899, 0.635, 0.078, 0.184, 30.308695, 0.0, 0.447, 0.741, 30.707491, 0.85, 0.325, 0.098, 31.106287000000002, 0.929, 0.694, 0.125, 31.505083000000003, 0.494, 0.184, 0.556, 31.903879, 0.466, 0.674, 0.188, 32.302675, 0.301, 0.745, 0.933, 32.701471, 0.635, 0.078, 0.184, 33.100267, 0.0, 0.447, 0.741, 33.499063, 0.85, 0.325, 0.098, 33.897859, 0.929, 0.694, 0.125, 34.296655, 0.494, 0.184, 0.556, 34.695451000000006, 0.466, 0.674, 0.188, 35.094246999999996, 0.301, 0.745, 0.933, 35.493043, 0.635, 0.078, 0.184, 35.89183899999999, 0.0, 0.447, 0.741, 36.290635, 0.85, 0.325, 0.098, 36.689431, 0.929, 0.694, 0.125, 37.088227, 0.494, 0.184, 0.556, 37.487023, 0.466, 0.674, 0.188, 37.885819, 0.301, 0.745, 0.933, 38.284615, 0.635, 0.078, 0.184, 38.68341100000001, 0.0, 0.447, 0.741, 39.082207, 0.85, 0.325, 0.098, 39.481002999999994, 0.929, 0.694, 0.125, 39.879799, 0.494, 0.184, 0.556, 40.278595, 0.466, 0.674, 0.188, 40.677391, 0.301, 0.745, 0.933, 41.076187000000004, 0.635, 0.078, 0.184, 41.474983, 0.0, 0.447, 0.741, 41.873779, 0.85, 0.325, 0.098, 42.272575, 0.929, 0.694, 0.125, 42.67137099999999, 0.494, 0.184, 0.556, 43.070167000000005, 0.466, 0.674, 0.188, 43.468962999999995, 0.301, 0.745, 0.933, 43.867759, 0.635, 0.078, 0.184, 44.266555000000004, 0.0, 0.447, 0.741, 44.665351, 0.85, 0.325, 0.098, 45.064147, 0.929, 0.694, 0.125, 45.462942999999996, 0.494, 0.184, 0.556, 45.861739, 0.466, 0.674, 0.188, 46.260535000000004, 0.301, 0.745, 0.933, 46.659330999999995, 0.635, 0.078, 0.184, 47.058127000000006, 0.0, 0.447, 0.741, 47.456922999999996, 0.85, 0.325, 0.098, 47.855719, 0.929, 0.694, 0.125, 48.254515, 0.494, 0.184, 0.556, 48.653311, 0.466, 0.674, 0.188, 49.052107, 0.301, 0.745, 0.933, 49.450903, 0.635, 0.078, 0.184, 49.849698999999994, 0.0, 0.447, 0.741, 50.248495, 0.85, 0.325, 0.098, 50.647290999999996, 0.929, 0.694, 0.125, 51.046087, 0.494, 0.184, 0.556, 51.444883, 0.466, 0.674, 0.188, 51.843679, 0.301, 0.745, 0.933, 52.242475, 0.635, 0.078, 0.184, 52.641271, 0.0, 0.447, 0.741, 53.040067, 0.85, 0.325, 0.098, 53.438863000000005, 0.929, 0.694, 0.125, 53.837658999999995, 0.494, 0.184, 0.556, 54.236455, 0.466, 0.674, 0.188, 54.635251, 0.301, 0.745, 0.933, 55.034047, 0.635, 0.078, 0.184, 55.432843, 0.0, 0.447, 0.741, 55.831639, 0.85, 0.325, 0.098, 56.230435, 0.929, 0.694, 0.125, 56.629231000000004, 0.494, 0.184, 0.556, 57.028027, 0.466, 0.674, 0.188, 57.42682299999999, 0.301, 0.745, 0.933, 57.825618999999996, 0.635, 0.078, 0.184, 58.22441499999999, 0.0, 0.447, 0.741, 58.623211, 0.85, 0.325, 0.098, 59.022007, 0.929, 0.694, 0.125, 59.420803, 0.494, 0.184, 0.556, 59.819599000000004, 0.466, 0.674, 0.188, 60.218395, 0.301, 0.745, 0.933, 60.617191000000005, 0.635, 0.078, 0.184, 61.015987, 0.0, 0.447, 0.741, 61.41478299999999, 0.85, 0.325, 0.098, 61.813579, 0.929, 0.694, 0.125, 62.212374999999994, 0.494, 0.184, 0.556, 62.611171, 0.466, 0.674, 0.188, 63.009966999999996, 0.301, 0.745, 0.933, 63.408763, 0.635, 0.078, 0.184, 63.807559000000005, 0.0, 0.447, 0.741, 64.206355, 0.85, 0.325, 0.098, 64.605151, 0.929, 0.694, 0.125, 65.00394700000001, 0.494, 0.184, 0.556, 65.402743, 0.466, 0.674, 0.188, 65.80153899999999, 0.301, 0.745, 0.933, 66.200335, 0.635, 0.078, 0.184, 66.599131, 0.0, 0.447, 0.741, 66.997927, 0.85, 0.325, 0.098, 67.396723, 0.929, 0.694, 0.125, 67.79551899999998, 0.494, 0.184, 0.556, 68.194315, 0.466, 0.674, 0.188, 68.59311100000001, 0.301, 0.745, 0.933, 68.991907, 0.635, 0.078, 0.184, 69.39070299999999, 0.0, 0.447, 0.741, 69.78949899999999, 0.85, 0.325, 0.098, 70.188295, 0.929, 0.694, 0.125, 70.587091, 0.494, 0.184, 0.556, 70.985887, 0.466, 0.674, 0.188, 71.38468300000001, 0.301, 0.745, 0.933, 71.783479, 0.635, 0.078, 0.184, 72.18227499999999, 0.0, 0.447, 0.741, 72.58107100000001, 0.85, 0.325, 0.098, 72.97986700000001, 0.929, 0.694, 0.125, 73.378663, 0.494, 0.184, 0.556, 73.777459, 0.466, 0.674, 0.188, 74.176255, 0.301, 0.745, 0.933, 74.575051, 0.635, 0.078, 0.184, 74.974046, 0.0, 0.447, 0.741, 75.37284199999999, 0.85, 0.325, 0.098, 75.771638, 0.929, 0.694, 0.125, 76.170434, 0.494, 0.184, 0.556, 76.56923, 0.466, 0.674, 0.188, 76.968026, 0.301, 0.745, 0.933, 77.36682200000001, 0.635, 0.078, 0.184, 77.765618, 0.0, 0.447, 0.741, 78.164414, 0.85, 0.325, 0.098, 78.56321, 0.929, 0.694, 0.125, 78.96200599999999, 0.494, 0.184, 0.556, 79.36080199999999, 0.466, 0.674, 0.188, 79.759598, 0.301, 0.745, 0.933, 80.158394, 0.635, 0.078, 0.184, 80.55719, 0.0, 0.447, 0.741, 80.95598600000001, 0.85, 0.325, 0.098, 81.354782, 0.929, 0.694, 0.125, 81.75357799999999, 0.494, 0.184, 0.556, 82.15237400000001, 0.466, 0.674, 0.188, 82.55117000000001, 0.301, 0.745, 0.933, 82.949966, 0.635, 0.078, 0.184, 83.348762, 0.0, 0.447, 0.741, 83.747558, 0.85, 0.325, 0.098, 84.146354, 0.929, 0.694, 0.125, 84.54515, 0.494, 0.184, 0.556, 84.943946, 0.466, 0.674, 0.188, 85.34274199999999, 0.301, 0.745, 0.933, 85.741538, 0.635, 0.078, 0.184, 86.14033400000001, 0.0, 0.447, 0.741, 86.53913, 0.85, 0.325, 0.098, 86.93792599999999, 0.929, 0.694, 0.125, 87.336722, 0.494, 0.184, 0.556, 87.735518, 0.466, 0.674, 0.188, 88.134314, 0.301, 0.745, 0.933, 88.53311000000001, 0.635, 0.078, 0.184, 88.93190600000001, 0.0, 0.447, 0.741, 89.330702, 0.85, 0.325, 0.098, 89.72949799999999, 0.929, 0.694, 0.125, 90.128294, 0.494, 0.184, 0.556, 90.52708999999999, 0.466, 0.674, 0.188, 90.92588599999999, 0.301, 0.745, 0.933, 91.324682, 0.635, 0.078, 0.184, 91.723478, 0.0, 0.447, 0.741, 92.122274, 0.85, 0.325, 0.098, 92.52107000000001, 0.929, 0.694, 0.125, 92.919866, 0.494, 0.184, 0.556, 93.31866199999999, 0.466, 0.674, 0.188, 93.71745800000001, 0.301, 0.745, 0.933, 94.11625400000001, 0.635, 0.078, 0.184, 94.51505, 0.0, 0.447, 0.741, 94.91384599999999, 0.85, 0.325, 0.098, 95.312642, 0.929, 0.694, 0.125, 95.711438, 0.494, 0.184, 0.556, 96.110234, 0.466, 0.674, 0.188, 96.50903, 0.301, 0.745, 0.933, 96.90782599999999, 0.635, 0.078, 0.184, 97.306622, 0.0, 0.447, 0.741, 97.70541800000001, 0.85, 0.325, 0.098, 98.104214, 0.929, 0.694, 0.125, 98.50300999999999, 0.494, 0.184, 0.556, 98.901806, 0.466, 0.674, 0.188, 99.300602, 0.301, 0.745, 0.933, 99.69939799999999, 0.635, 0.078, 0.184, 100.09819399999999, 0.0, 0.447, 0.741, 100.49699, 0.85, 0.325, 0.098, 100.89578599999999, 0.929, 0.694, 0.125, 101.29458199999999, 0.494, 0.184, 0.556, 101.693378, 0.466, 0.674, 0.188, 102.092174, 0.301, 0.745, 0.933, 102.49097, 0.635, 0.078, 0.184, 102.889766, 0.0, 0.447, 0.741, 103.288562, 0.85, 0.325, 0.098, 103.687358, 0.929, 0.694, 0.125, 104.08615400000001, 0.494, 0.184, 0.556, 104.48495, 0.466, 0.674, 0.188, 104.883746, 0.301, 0.745, 0.933, 105.282542, 0.635, 0.078, 0.184, 105.68133800000001, 0.0, 0.447, 0.741, 106.080134, 0.85, 0.325, 0.098, 106.47893, 0.929, 0.694, 0.125, 106.87772600000001, 0.494, 0.184, 0.556, 107.27652199999999, 0.466, 0.674, 0.188, 107.67531799999999, 0.301, 0.745, 0.933, 108.074114, 0.635, 0.078, 0.184, 108.47291, 0.0, 0.447, 0.741, 108.87170599999999, 0.85, 0.325, 0.098, 109.270502, 0.929, 0.694, 0.125, 109.669298, 0.494, 0.184, 0.556, 110.068094, 0.466, 0.674, 0.188, 110.46688999999999, 0.301, 0.745, 0.933, 110.865686, 0.635, 0.078, 0.184, 111.264482, 0.0, 0.447, 0.741, 111.663278, 0.85, 0.325, 0.098, 112.06207400000001, 0.929, 0.694, 0.125, 112.46087, 0.494, 0.184, 0.556, 112.859666, 0.466, 0.674, 0.188, 113.25846200000001, 0.301, 0.745, 0.933, 113.65725800000001, 0.635, 0.078, 0.184, 114.056054, 0.0, 0.447, 0.741, 114.45485000000001, 0.85, 0.325, 0.098, 114.85364599999998, 0.929, 0.694, 0.125, 115.25244199999999, 0.494, 0.184, 0.556, 115.65123799999999, 0.466, 0.674, 0.188, 116.050034, 0.301, 0.745, 0.933, 116.44882999999999, 0.635, 0.078, 0.184, 116.84762599999999, 0.0, 0.447, 0.741, 117.246422, 0.85, 0.325, 0.098, 117.645218, 0.929, 0.694, 0.125, 118.044014, 0.494, 0.184, 0.556, 118.44281, 0.466, 0.674, 0.188, 118.841606, 0.301, 0.745, 0.933, 119.240402, 0.635, 0.078, 0.184, 119.63919800000001, 0.0, 0.447, 0.741, 120.037994, 0.85, 0.325, 0.098, 120.43679, 0.929, 0.694, 0.125, 120.835586, 0.494, 0.184, 0.556, 121.23438200000001, 0.466, 0.674, 0.188, 121.63317800000002, 0.301, 0.745, 0.933, 122.031974, 0.635, 0.078, 0.184, 122.43077000000001, 0.0, 0.447, 0.741, 122.82956599999999, 0.85, 0.325, 0.098, 123.22836199999999, 0.929, 0.694, 0.125, 123.627158, 0.494, 0.184, 0.556, 124.025954, 0.466, 0.674, 0.188, 124.424949, 0.301, 0.745, 0.933, 124.823745, 0.635, 0.078, 0.184, 125.222541, 0.0, 0.447, 0.741, 125.621337, 0.85, 0.325, 0.098, 126.020133, 0.929, 0.694, 0.125, 126.418929, 0.494, 0.184, 0.556, 126.81772500000001, 0.466, 0.674, 0.188, 127.21652100000001, 0.301, 0.745, 0.933, 127.615317, 0.635, 0.078, 0.184, 128.014113, 0.0, 0.447, 0.741, 128.41290899999998, 0.85, 0.325, 0.098, 128.811705, 0.929, 0.694, 0.125, 129.210501, 0.494, 0.184, 0.556, 129.609297, 0.466, 0.674, 0.188, 130.008093, 0.301, 0.745, 0.933, 130.406889, 0.635, 0.078, 0.184, 130.80568499999998, 0.0, 0.447, 0.741, 131.204481, 0.85, 0.325, 0.098, 131.603277, 0.929, 0.694, 0.125, 132.002073, 0.494, 0.184, 0.556, 132.400869, 0.466, 0.674, 0.188, 132.799665, 0.301, 0.745, 0.933, 133.198461, 0.635, 0.078, 0.184, 133.597257, 0.0, 0.447, 0.741, 133.99605300000002, 0.85, 0.325, 0.098, 134.394849, 0.929, 0.694, 0.125, 134.793645, 0.494, 0.184, 0.556, 135.192441, 0.466, 0.674, 0.188, 135.59123700000004, 0.301, 0.745, 0.933, 135.990033, 0.635, 0.078, 0.184, 136.388829, 0.0, 0.447, 0.741, 136.787625, 0.85, 0.325, 0.098, 137.186421, 0.929, 0.694, 0.125, 137.58521699999997, 0.494, 0.184, 0.556, 137.984013, 0.466, 0.674, 0.188, 138.382809, 0.301, 0.745, 0.933, 138.781605, 0.635, 0.078, 0.184, 139.180401, 0.0, 0.447, 0.741, 139.579197, 0.85, 0.325, 0.098, 139.977993, 0.929, 0.694, 0.125, 140.376789, 0.494, 0.184, 0.556, 140.775585, 0.466, 0.674, 0.188, 141.174381, 0.301, 0.745, 0.933, 141.573177, 0.635, 0.078, 0.184, 141.97197300000002, 0.0, 0.447, 0.741, 142.370769, 0.85, 0.325, 0.098, 142.76956500000003, 0.929, 0.694, 0.125, 143.168361, 0.494, 0.184, 0.556, 143.56715699999998, 0.466, 0.674, 0.188, 143.96595299999998, 0.301, 0.745, 0.933, 144.364749, 0.635, 0.078, 0.184, 144.76354499999997, 0.0, 0.447, 0.741, 145.162341, 0.85, 0.325, 0.098, 145.561137, 0.929, 0.694, 0.125, 145.959933, 0.494, 0.184, 0.556, 146.358729, 0.466, 0.674, 0.188, 146.757525, 0.301, 0.745, 0.933, 147.156321, 0.635, 0.078, 0.184, 147.55511699999997, 0.0, 0.447, 0.741, 147.953913, 0.85, 0.325, 0.098, 148.352709, 0.929, 0.694, 0.125, 148.75150500000004, 0.494, 0.184, 0.556, 149.150301, 0.466, 0.674, 0.188, 149.54909700000002, 0.301, 0.745, 0.933, 149.947893, 0.635, 0.078, 0.184, 150.346689, 0.0, 0.447, 0.741, 150.745485, 0.85, 0.325, 0.098, 151.144281, 0.929, 0.694, 0.125, 151.54307699999998, 0.494, 0.184, 0.556, 151.94187299999996, 0.466, 0.674, 0.188, 152.340669, 0.301, 0.745, 0.933, 152.739465, 0.635, 0.078, 0.184, 153.13826100000003, 0.0, 0.447, 0.741, 153.537057, 0.85, 0.325, 0.098, 153.935853, 0.929, 0.694, 0.125, 154.33464899999998, 0.494, 0.184, 0.556, 154.733445, 0.466, 0.674, 0.188, 155.132241, 0.301, 0.745, 0.933, 155.531037, 0.635, 0.078, 0.184, 155.92983300000003, 0.0, 0.447, 0.741, 156.328629, 0.85, 0.325, 0.098, 156.727425, 0.929, 0.694, 0.125, 157.12622100000002, 0.494, 0.184, 0.556, 157.525017, 0.466, 0.674, 0.188, 157.923813, 0.301, 0.745, 0.933, 158.322609, 0.635, 0.078, 0.184, 158.72140500000003, 0.0, 0.447, 0.741, 159.120201, 0.85, 0.325, 0.098, 159.51899699999998, 0.929, 0.694, 0.125, 159.917793, 0.494, 0.184, 0.556, 160.316589, 0.466, 0.674, 0.188, 160.71538499999997, 0.301, 0.745, 0.933, 161.114181, 0.635, 0.078, 0.184, 161.512977, 0.0, 0.447, 0.741, 161.911773, 0.85, 0.325, 0.098, 162.310569, 0.929, 0.694, 0.125, 162.709365, 0.494, 0.184, 0.556, 163.108161, 0.466, 0.674, 0.188, 163.506957, 0.301, 0.745, 0.933, 163.905753, 0.635, 0.078, 0.184, 164.304549, 0.0, 0.447, 0.741, 164.70334500000004, 0.85, 0.325, 0.098, 165.10214100000002, 0.929, 0.694, 0.125, 165.50093700000002, 0.494, 0.184, 0.556, 165.899733, 0.466, 0.674, 0.188, 166.298529, 0.301, 0.745, 0.933, 166.69732499999998, 0.635, 0.078, 0.184, 167.09612099999998, 0.0, 0.447, 0.741, 167.494917, 0.85, 0.325, 0.098, 167.89371299999996, 0.929, 0.694, 0.125, 168.292509, 0.494, 0.184, 0.556, 168.691305, 0.466, 0.674, 0.188, 169.09010100000003, 0.301, 0.745, 0.933, 169.488897, 0.635, 0.078, 0.184, 169.887693, 0.0, 0.447, 0.741, 170.286489, 0.85, 0.325, 0.098, 170.68528499999996, 0.929, 0.694, 0.125, 171.084081, 0.494, 0.184, 0.556, 171.482877, 0.466, 0.674, 0.188, 171.88167300000003, 0.301, 0.745, 0.933, 172.280469, 0.635, 0.078, 0.184, 172.67926500000002, 0.0, 0.447, 0.741, 173.07806100000002, 0.85, 0.325, 0.098, 173.476857, 0.929, 0.694, 0.125, 173.875653, 0.494, 0.184, 0.556, 174.27464799999998, 0.466, 0.674, 0.188, 174.673444, 0.301, 0.745, 0.933, 175.07224000000002, 0.635, 0.078, 0.184, 175.471036, 0.0, 0.447, 0.741, 175.869832, 0.85, 0.325, 0.098, 176.268628, 0.929, 0.694, 0.125, 176.66742399999998, 0.494, 0.184, 0.556, 177.06622000000002, 0.466, 0.674, 0.188, 177.465016, 0.301, 0.745, 0.933, 177.86381200000002, 0.635, 0.078, 0.184, 178.262608, 0.0, 0.447, 0.741, 178.661404, 0.85, 0.325, 0.098, 179.0602, 0.929, 0.694, 0.125, 179.45899599999998, 0.494, 0.184, 0.556, 179.857792, 0.466, 0.674, 0.188, 180.256588, 0.301, 0.745, 0.933, 180.655384, 0.635, 0.078, 0.184, 181.05417999999997, 0.0, 0.447, 0.741, 181.452976, 0.85, 0.325, 0.098, 181.85177199999998, 0.929, 0.694, 0.125, 182.25056800000002, 0.494, 0.184, 0.556, 182.649364, 0.466, 0.674, 0.188, 183.04816, 0.301, 0.745, 0.933, 183.446956, 0.635, 0.078, 0.184, 183.84575199999998, 0.0, 0.447, 0.741, 184.244548, 0.85, 0.325, 0.098, 184.643344, 0.929, 0.694, 0.125, 185.04214000000002, 0.494, 0.184, 0.556, 185.440936, 0.466, 0.674, 0.188, 185.839732, 0.301, 0.745, 0.933, 186.238528, 0.635, 0.078, 0.184, 186.63732399999998, 0.0, 0.447, 0.741, 187.03612, 0.85, 0.325, 0.098, 187.43491600000002, 0.929, 0.694, 0.125, 187.833712, 0.494, 0.184, 0.556, 188.23250800000002, 0.466, 0.674, 0.188, 188.631304, 0.301, 0.745, 0.933, 189.0301, 0.635, 0.078, 0.184, 189.42889599999998, 0.0, 0.447, 0.741, 189.82769199999998, 0.85, 0.325, 0.098, 190.226488, 0.929, 0.694, 0.125, 190.625284, 0.494, 0.184, 0.556, 191.02408000000003, 0.466, 0.674, 0.188, 191.422876, 0.301, 0.745, 0.933, 191.821672, 0.635, 0.078, 0.184, 192.220468, 0.0, 0.447, 0.741, 192.61926400000002, 0.85, 0.325, 0.098, 193.01806, 0.929, 0.694, 0.125, 193.416856, 0.494, 0.184, 0.556, 193.81565199999997, 0.466, 0.674, 0.188, 194.214448, 0.301, 0.745, 0.933, 194.613244, 0.635, 0.078, 0.184, 195.01204, 0.0, 0.447, 0.741, 195.41083600000002, 0.85, 0.325, 0.098, 195.809632, 0.929, 0.694, 0.125, 196.208428, 0.494, 0.184, 0.556, 196.607224, 0.466, 0.674, 0.188, 197.00601999999998, 0.301, 0.745, 0.933, 197.40481599999998, 0.635, 0.078, 0.184, 197.803612, 0.0, 0.447, 0.741, 198.20240800000002, 0.85, 0.325, 0.098, 198.601204, 0.929, 0.694, 0.125, 199.0, 0.494, 0.184, 0.556]
    ancestorsLUT.ColorSpace = 'RGB'
    ancestorsLUT.NumberOfTableValues = 258 #296
    ancestorsLUT.ScalarRangeInitialized = 1.0

    
    # hide color bar/color legend
    glyph1Display.SetScalarBarVisibility(renderView1, False)
    
    # update the view to ensure updated data information
    renderView1.Update()
    
    # reset view to fit data
    renderView1.ResetCamera()
        
    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    ancestorsLUT.ApplyPreset('periodic', True)
    
    animationScene1.GoToLast()
    
    # reset view to fit data
    renderView1.ResetCamera()
    
#    # get opacity transfer function/opacity map for 'GlyphScale'
#    glyphScalePWF = GetOpacityTransferFunction('GlyphScale')
#    glyphScalePWF.Points = [0.0, 0.0, 0.5, 0.0, 95.0, 1.0, 0.5, 0.0]
#    glyphScalePWF.ScalarRangeInitialized = 1
    # get opacity transfer function/opacity map for 'Ancestors'
    ancestorsPWF = GetOpacityTransferFunction('Ancestors')
    ancestorsPWF.Points = [0.0, 0.0, 0.5, 0.0, 199.0, 1.0, 0.5, 0.0]
    ancestorsPWF.ScalarRangeInitialized = 1
    
    
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
    
    # reset view to fit data
    renderView1.ResetCamera()
    
    #animationScene1.GoToLast()
    
    # current camera placement for renderView1
    renderView1.CameraPosition = [29.51439755018888, 8.852873533964157, 1.7701021134853363]
    renderView1.CameraFocalPoint = [5.0945106744766235, 8.852873533964157, 1.7701021134853363]
    renderView1.CameraViewUp = [0.0, 0.0, 1.0]
    renderView1.CameraParallelScale = 11.196853328693654
    
    tk = GetTimeKeeper()
    timesteps = tk.TimestepValues
    animationScene1.AnimationTime = timesteps[time_step]
    
    # save screenshot
    #SaveScreenshot(path+'plus_x_t'+str(time_step)+'.png', renderView1, ImageResolution=[1625, 618], TransparentBackground=1)
    
    # create a new 'Threshold'
    threshold1 = Threshold(Input=glyph1)
    threshold1.Scalars = ['POINTS', 'Ancestors']
    threshold1.ThresholdRange = [0.0, 100.0]
    
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
    threshold1Display.ScalarOpacityUnitDistance =0.5074161884921996
    
    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    threshold1Display.OSPRayScaleFunction.Points = [-1.0, 0.0, 0.5, 0.0, 147.0, 1.0, 0.5, 0.0]
    
    # hide data in view
    Hide(glyph1, renderView1)
    
    # show color bar/color legend
    threshold1Display.SetScalarBarVisibility(renderView1, False)
    
    
    # update the view to ensure updated data information
    renderView1.Update()
    
    # set active source
    SetActiveSource(threshold1)
    
    # reset view to fit data
    renderView1.ResetCamera()
       
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
    renderView1.CameraPosition = [4.094428837299347, 25.090791181239215, 1.75]
    renderView1.CameraFocalPoint = [4.094428837299347, 5.0814430713653564, 1.75]
    renderView1.CameraViewUp = [0.0, 0.0, 1.0]
    renderView1.CameraParallelScale = 9.174560763910034
    
    # save screenshot
    SaveScreenshot(path+'threshold_48_minus_y_t'+str(time_step)+'.png', renderView1, ImageResolution=[3250, 1236], TransparentBackground=1)
    
    animationScene1.AnimationTime = timesteps[0]
    # save screenshot
    SaveScreenshot(path+'threshold_48_minus_y_t'+str(0)+'.png', renderView1, ImageResolution=[3250, 1236], TransparentBackground=1)
    
     
    #### uncomment the following to render all views
    #RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).
    
if __name__ == "__main__":
   screenshot(sys.argv[1:])    
