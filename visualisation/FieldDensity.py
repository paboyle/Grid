#!/usr/bin/env python

# noinspection PyUnresolvedReferences
import math
import vtk
import vtkmodules.vtkInteractionStyle
# noinspection PyUnresolvedReferences
import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkCommonCore import (
    VTK_VERSION_NUMBER,
    vtkVersion
)
from vtkmodules.vtkCommonCore import VTK_DOUBLE
from vtkmodules.vtkCommonDataModel import vtkImageData
from vtkmodules.vtkFiltersCore import (
    vtkMarchingCubes,
    vtkStripper
)
from vtkmodules.vtkFiltersModeling import vtkOutlineFilter
from vtkmodules.vtkIOImage import (
    vtkMetaImageReader,
    vtkJPEGWriter,
    vtkPNGWriter
)
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkCamera,
    vtkPolyDataMapper,
    vtkProperty,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkRenderer,
    vtkWindowToImageFilter
)


class vtkTimerCallback():
    def __init__(self, steps, imageData, iren):
        self.timer_count = 0
        self.steps = steps
        self.imageData = imageData
        self.iren = iren
        self.timerId = None
        self.step = 0

    def execute(self, obj, event):

        print(self.timer_count)

        dims = self.imageData.GetDimensions()

        t=self.step/10.0

        z0 = 2
        y0 = 4
        x0 = 4
        z1 = 14
        y1 = 12
        x1 = 12
        for z in range(dims[2]):
            for y in range(dims[1]):
                for x in range(dims[0]):
                    self.imageData.SetScalarComponentFromDouble(x, y, z, 0,
                                                                math.sin(t)*math.exp(-0.25*((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0)))
                                                                - math.cos(t)*math.exp(-0.25*((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1))))

        self.imageData.Modified()

        iren = obj
        iren.GetRenderWindow().Render()
        self.timer_count += 1
        self.step += 1

        if self.step >= self.steps :
            iren.DestroyTimer(self.timerId)

def WriteImage(fileName, renWin):
    '''
    '''

    import os

    if fileName:
        # Select the writer to use.
        path, ext = os.path.splitext(fileName)
        ext = ext.lower()
        if not ext:
            ext = '.png'
            fileName = fileName + ext
        elif ext == '.jpg':
            writer = vtkJPEGWriter()
        else:
            writer = vtkPNGWriter()

        windowto_image_filter = vtkWindowToImageFilter()
        windowto_image_filter.SetInput(renWin)
        windowto_image_filter.SetScale(1)  # image quality
        windowto_image_filter.SetInputBufferTypeToRGBA()

        writer.SetFileName(fileName)
        writer.SetInputConnection(windowto_image_filter.GetOutputPort())
        writer.Write()
    else:
        raise RuntimeError('Need a filename.')


def main():

    colors = vtkNamedColors()

    file_name = get_program_parameters()

    colors.SetColor('InstantonColor', [240, 184, 160, 255])
    colors.SetColor('BackfaceColor', [255, 229, 200, 255])
    colors.SetColor('BkgColor', [51, 77, 102, 255])

    # Create the renderer, the render window, and the interactor. The renderer
    # draws into the render window, the interactor enables mouse- and
    # keyboard-based interaction with the data within the render window.
    #
    a_renderer = vtkRenderer()
    ren_win = vtkRenderWindow()
    ren_win.AddRenderer(a_renderer)

    iren = vtkRenderWindowInteractor()
    iren.SetRenderWindow(ren_win)

    # The following reader is used to read a series of 2D slices (images)
    # that compose the volume. The slice dimensions are set, and the
    # pixel spacing. The data Endianness must also be specified. The reader
    # uses the FilePrefix in combination with the slice number to construct
    # filenames using the format FilePrefix.%d. (In this case the FilePrefix
    # is the root name of the file: quarter.)
    imageData = vtkImageData()
    imageData.SetDimensions(16, 16, 16)
    imageData.AllocateScalars(VTK_DOUBLE, 1)

    dims = imageData.GetDimensions()

    # Fill every entry of the image data with '2.0'
    # Set the input data
    for z in range(dims[2]):
        z0 = dims[2]/2
        for y in range(dims[1]):
            y0 = dims[1]/2
            for x in range(dims[0]):
                x0 = dims[0]/2
                imageData.SetScalarComponentFromDouble(x, y, z, 0, math.exp(-0.25*((x-x0)*(x-x0)+(y-y0)*(y-y0)+z*z)) -  math.exp(-0.25*((x-x0)*(x-x0)+y*y+(z-z0)*(z-z0))))

    instanton_extractor = vtkMarchingCubes()
    instanton_extractor.SetInputData(imageData)
    instanton_extractor.SetValue(0, 0.1)
    
    instanton_stripper = vtkStripper()
    instanton_stripper.SetInputConnection(instanton_extractor.GetOutputPort())
    
    instanton_mapper = vtkPolyDataMapper()
    instanton_mapper.SetInputConnection(instanton_stripper.GetOutputPort())
    instanton_mapper.ScalarVisibilityOff()
    
    instanton = vtkActor()
    instanton.SetMapper(instanton_mapper)
    instanton.GetProperty().SetDiffuseColor(colors.GetColor3d('InstantonColor'))
    instanton.GetProperty().SetSpecular(0.3)
    instanton.GetProperty().SetSpecularPower(20)
    instanton.GetProperty().SetOpacity(0.5)

    # The triangle stripper is used to create triangle strips from the
    # isosurface these render much faster on may systems.
    antiinstanton_extractor = vtkMarchingCubes()
    antiinstanton_extractor.SetInputData(imageData)
    antiinstanton_extractor.SetValue(0, -0.1)
    
    antiinstanton_stripper = vtkStripper()
    antiinstanton_stripper.SetInputConnection(antiinstanton_extractor.GetOutputPort())
    
    antiinstanton_mapper = vtkPolyDataMapper()
    antiinstanton_mapper.SetInputConnection(antiinstanton_stripper.GetOutputPort())
    antiinstanton_mapper.ScalarVisibilityOff()
    
    antiinstanton = vtkActor()
    antiinstanton.SetMapper(antiinstanton_mapper)
    antiinstanton.GetProperty().SetDiffuseColor(colors.GetColor3d('Ivory'))

    # An outline provides box around the data.
    outline_data = vtkOutlineFilter()
    outline_data.SetInputData(imageData)
    
    map_outline = vtkPolyDataMapper()
    map_outline.SetInputConnection(outline_data.GetOutputPort())
    
    outline = vtkActor()
    outline.SetMapper(map_outline)
    outline.GetProperty().SetColor(colors.GetColor3d('Black'))
    
    # It is convenient to create an initial view of the data. The FocalPoint
    # and Position form a vector direction. Later on (ResetCamera() method)
    # this vector is used to position the camera to look at the data in
    # this direction.
    a_camera = vtkCamera()
    a_camera.SetViewUp(0, 0, -1)
    a_camera.SetPosition(0, -100, 0)
    a_camera.SetFocalPoint(0, 0, 0)
    a_camera.ComputeViewPlaneNormal()
    a_camera.Azimuth(30.0)
    a_camera.Elevation(30.0)

    # Actors are added to the renderer. An initial camera view is created.
    # The Dolly() method moves the camera towards the FocalPoint,
    # thereby enlarging the image.
    a_renderer.AddActor(outline)
    a_renderer.AddActor(instanton)
    a_renderer.AddActor(antiinstanton)
    a_renderer.SetActiveCamera(a_camera)
    a_renderer.ResetCamera()
    a_camera.Dolly(1.0)

    # Set a background color for the renderer and set the size of the
    # render window (expressed in pixels).
    a_renderer.SetBackground(colors.GetColor3d('BkgColor'))
    ren_win.SetSize(1024, 1024)
    ren_win.SetWindowName('ExpoDemo')

    # Note that when camera movement occurs (as it does in the Dolly()
    # method), the clipping planes often need adjusting. Clipping planes
    # consist of two planes: near and far along the view direction. The
    # near plane clips out objects in front of the plane the far plane
    # clips out objects behind the plane. This way only what is drawn
    # between the planes is actually rendered.
    a_renderer.ResetCameraClippingRange()

    # write image
    #    WriteImage('exp.jpg',ren_win)

    # Sign up to receive TimerEvent
    cb = vtkTimerCallback(200, imageData, iren)
    iren.AddObserver('TimerEvent', cb.execute)
    cb.timerId = iren.CreateRepeatingTimer(50)

    # start the interaction and timer
    ren_win.Render()
    
    # Initialize the event loop and then start it.
    iren.Initialize()
    iren.Start()


def get_program_parameters():
    import argparse
    description = 'Simple lattice volumetric demo'
    epilogue = '''
    Derived from VTK/Examples/Cxx/Medical2.cxx
    '''
    parser = argparse.ArgumentParser(description=description, epilog=epilogue,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename', help='FieldDensity.py')
    args = parser.parse_args()
    return args.filename


def vtk_version_ok(major, minor, build):
    """
    Check the VTK version.

    :param major: Major version.
    :param minor: Minor version.
    :param build: Build version.
    :return: True if the requested VTK version is greater or equal to the actual VTK version.
    """
    needed_version = 10000000000 * int(major) + 100000000 * int(minor) + int(build)
    try:
        vtk_version_number = VTK_VERSION_NUMBER
    except AttributeError:  # as error:
        ver = vtkVersion()
        vtk_version_number = 10000000000 * ver.GetVTKMajorVersion() + 100000000 * ver.GetVTKMinorVersion() \
                             + ver.GetVTKBuildVersion()
    if vtk_version_number >= needed_version:
        return True
    else:
        return False


if __name__ == '__main__':
    main()
