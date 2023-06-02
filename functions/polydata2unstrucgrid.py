# trace generated using paraview version 5.10.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *

import argparse
import sys


def get_arguments(input_args):
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("input_file_name", help='Enter output name with extension')
    parser.add_argument("path_in", help='Enter full input path')
    parser.add_argument("path_out", help='Enter full output path')
    return parser.parse_args()


if __name__ == "__main__":
    
    args = get_arguments(sys.argv)
    filename = args.input_file_name
    path_in = args.path_in
    path_out = args.path_out
    
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'Legacy VTK Reader'
    struct1_septum_boundaryvtk = LegacyVTKReader(registrationName=f'{filename}', FileNames=[f'{path_in}'])

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    struct1_septum_boundaryvtkDisplay = Show(struct1_septum_boundaryvtk, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    struct1_septum_boundaryvtkDisplay.Representation = 'Surface'
    struct1_septum_boundaryvtkDisplay.ColorArrayName = [None, '']
    struct1_septum_boundaryvtkDisplay.SelectTCoordArray = 'None'
    struct1_septum_boundaryvtkDisplay.SelectNormalArray = 'None'
    struct1_septum_boundaryvtkDisplay.SelectTangentArray = 'None'
    struct1_septum_boundaryvtkDisplay.OSPRayScaleArray = 'guess'
    struct1_septum_boundaryvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    struct1_septum_boundaryvtkDisplay.SelectOrientationVectors = 'None'
    struct1_septum_boundaryvtkDisplay.ScaleFactor = 5.9568474903
    struct1_septum_boundaryvtkDisplay.SelectScaleArray = 'None'
    struct1_septum_boundaryvtkDisplay.GlyphType = 'Arrow'
    struct1_septum_boundaryvtkDisplay.GlyphTableIndexArray = 'None'
    struct1_septum_boundaryvtkDisplay.GaussianRadius = 0.297842374515
    struct1_septum_boundaryvtkDisplay.SetScaleArray = ['POINTS', 'guess']
    struct1_septum_boundaryvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    struct1_septum_boundaryvtkDisplay.OpacityArray = ['POINTS', 'guess']
    struct1_septum_boundaryvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    struct1_septum_boundaryvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
    struct1_septum_boundaryvtkDisplay.PolarAxes = 'PolarAxesRepresentation'

    # reset view to fit data
    renderView1.ResetCamera(False)

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Clean to Grid'
    cleantoGrid1 = CleantoGrid(registrationName='CleantoGrid1', Input=struct1_septum_boundaryvtk)

    # show data in view
    cleantoGrid1Display = Show(cleantoGrid1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    cleantoGrid1Display.Representation = 'Surface'
    cleantoGrid1Display.ColorArrayName = [None, '']
    cleantoGrid1Display.SelectTCoordArray = 'None'
    cleantoGrid1Display.SelectNormalArray = 'None'
    cleantoGrid1Display.SelectTangentArray = 'None'
    cleantoGrid1Display.OSPRayScaleArray = 'guess'
    cleantoGrid1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    cleantoGrid1Display.SelectOrientationVectors = 'None'
    cleantoGrid1Display.ScaleFactor = 5.9568474903
    cleantoGrid1Display.SelectScaleArray = 'None'
    cleantoGrid1Display.GlyphType = 'Arrow'
    cleantoGrid1Display.GlyphTableIndexArray = 'None'
    cleantoGrid1Display.GaussianRadius = 0.297842374515
    cleantoGrid1Display.SetScaleArray = ['POINTS', 'guess']
    cleantoGrid1Display.ScaleTransferFunction = 'PiecewiseFunction'
    cleantoGrid1Display.OpacityArray = ['POINTS', 'guess']
    cleantoGrid1Display.OpacityTransferFunction = 'PiecewiseFunction'
    cleantoGrid1Display.DataAxesGrid = 'GridAxesRepresentation'
    cleantoGrid1Display.PolarAxes = 'PolarAxesRepresentation'
    cleantoGrid1Display.ScalarOpacityUnitDistance = 7.5310525534515165
    cleantoGrid1Display.OpacityArrayName = ['POINTS', 'guess']

    # hide data in view
    Hide(struct1_septum_boundaryvtk, renderView1)

    # update the view to ensure updated data information
    renderView1.Update()

    # save data
    SaveData(f'{path_out}', proxy=cleantoGrid1, PointDataArrays=['guess', 'initial'],
        FileType='Ascii')

    #================================================================
    # addendum: following script captures some of the application
    # state to faithfully reproduce the visualization during playback
    #================================================================

    # get layout
    layout1 = GetLayout()

    #--------------------------------
    # saving layout sizes for layouts

    # layout/tab size in pixels
    layout1.SetSize(2638, 1794)

    #-----------------------------------
    # saving camera placements for views

    # current camera placement for renderView1
    renderView1.CameraPosition = [-24.22396884705, -1.7270419979999998, 154.4061303626508]
    renderView1.CameraFocalPoint = [-24.22396884705, -1.7270419979999998, 16.158228806500002]
    renderView1.CameraParallelScale = 35.78118986819024

    #--------------------------------------------
    # uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).