// Derived from VTK/Examples/Cxx/Medical2.cxx
// The example reads a volume dataset, extracts two isosurfaces that
// represent the skin and bone, and then displays them.
//
// Modified heavily by Peter Boyle to display lattice field theory data as movies and compare multiple files

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkMetaImageReader.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkOutlineFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkStripper.h>
#include <vtkImageData.h>
#include <vtkVersion.h>
#include <vtkCallbackCommand.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>

#define MPEG
#ifdef MPEG
#include <vtkFFMPEGWriter.h>
#endif

#include <vtkProperty2D.h>
#include <vtkSliderWidget.h>
#include <vtkSliderRepresentation2D.h>
#include <vtkWindowToImageFilter.h>

#include <array>
#include <string>

#include <Grid/Grid.h>

#define USE_FLYING_EDGES
#ifdef USE_FLYING_EDGES
#include <vtkFlyingEdges3D.h>
typedef vtkFlyingEdges3D isosurface;
#else
#include <vtkMarchingCubes.h>
typedef vtkMarchingCubes isosurface;
#endif

int mpeg = 0 ;
int xlate = 0 ;
int framerate = 10;

template <class T> void readFile(T& out, std::string const fname){
  Grid::emptyUserRecord record;
  Grid::ScidacReader RD;
  RD.open(fname);
  RD.readScidacFieldRecord(out,record);
  RD.close();
}
using namespace Grid;

class FrameUpdater : public vtkCallbackCommand
{
public:

  FrameUpdater() {
    TimerCount = 0;
    xoff       = 0;
    t          = 0;
    imageData = nullptr;
    grid_data = nullptr;
    timerId = 0;
    maxCount = -1;
  }
  
  static FrameUpdater* New()
  {
    FrameUpdater* cb = new FrameUpdater;
    cb->TimerCount = 0;
    return cb;
  }

  virtual void Execute(vtkObject* caller, unsigned long eventId,void* vtkNotUsed(callData))
  {
    const int max=256;
    char text_string[max];

    if (this->TimerCount < this->maxCount) {

      if (vtkCommand::TimerEvent == eventId)
	{
	  ++this->TimerCount;
	  
	  // Make a new frame
	  auto latt_size = grid_data->Grid()->GlobalDimensions();
	  for(int xx=0;xx<latt_size[0];xx++){
	    for(int yy=0;yy<latt_size[1];yy++){
	      for(int zz=0;zz<latt_size[2];zz++){
		int x = (xx+xoff)%latt_size[0];
		Coordinate site({x,yy,zz,t});
		RealD value = real(peekSite(*grid_data,site));
		imageData->SetScalarComponentFromDouble(xx,yy,zz,0,value);
	  }}}

	  if ( xlate ) { 
	    xoff = (xoff + 1)%latt_size[0];
	    if ( xoff== 0 ) t = (t+1)%latt_size[3];
	  } else {
	    t = (t+1)%latt_size[3];
	    if ( t== 0 ) 	xoff = (xoff + 1)%latt_size[0];
	  }

	  snprintf(text_string,max,"T=%d",t);
	  text->SetInput(text_string);
      
	  std::cout << this->TimerCount<<"/"<<maxCount<< " xoff "<<xoff<<" t "<<t  <<std::endl;
	  imageData->Modified();

	  vtkRenderWindowInteractor* iren = dynamic_cast<vtkRenderWindowInteractor*>(caller);
	  iren->GetRenderWindow()->Render();
	  
	}
    }
    
    if (this->TimerCount >= this->maxCount) {
      vtkRenderWindowInteractor* iren = dynamic_cast<vtkRenderWindowInteractor*>(caller);
      if (this->timerId > -1)
      {
        iren->DestroyTimer(this->timerId);
      }
    }
  }
  

private:
  int TimerCount;
  int xoff;
  int t;
public:
  Grid::LatticeComplexD * grid_data;
  vtkImageData* imageData = nullptr;
  vtkTextActor* text = nullptr;
  vtkFFMPEGWriter *writer = nullptr;
  int timerId ;
  int maxCount ;
  double rms;
  isosurface * posExtractor;
  isosurface * negExtractor;
};
class SliderCallback : public vtkCommand
{
public:
    static SliderCallback* New()
    {
        return new SliderCallback;
    }
    virtual void Execute(vtkObject* caller, unsigned long eventId, void* callData)
    {
        vtkSliderWidget *sliderWidget = vtkSliderWidget::SafeDownCast(caller);
        if (sliderWidget)
        {
	  contour = ((vtkSliderRepresentation *)sliderWidget->GetRepresentation())->GetValue();
        }
	for(int i=0;i<fu_list.size();i++){
	  fu_list[i]->posExtractor->SetValue(0,  SliderCallback::contour*fu_list[i]->rms);
	  fu_list[i]->negExtractor->SetValue(0, -SliderCallback::contour*fu_list[i]->rms);
	  fu_list[i]->posExtractor->Modified();
	  fu_list[i]->negExtractor->Modified();
	}
    }
public:
  static double contour;
  std::vector<FrameUpdater *> fu_list;
};


double SliderCallback::contour;

int main(int argc, char* argv[])
{
  using namespace Grid;

  Grid_init(&argc, &argv);
  GridLogLayout();

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  GridCartesian    Grid(latt_size, simd_layout, mpi_layout);

 
  std::cout << argc << " command Line arguments "<<std::endl;
  for(int c=0;c<argc;c++) {
    std::cout << " - "<<argv[c]<<std::endl;
  }

  std::vector<std::string> file_list;
  double default_contour = 1.0;
  std::string arg;
#ifdef MPEG
  if( GridCmdOptionExists(argv,argv+argc,"--mpeg") ){
    mpeg = 1;
  }
#endif

  if( GridCmdOptionExists(argv,argv+argc,"--xlate") ){
    xlate = 1;
  }

  if( GridCmdOptionExists(argv,argv+argc,"--fps") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--fps");
    GridCmdOptionInt(arg,framerate);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--isosurface") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--isosurface");
    GridCmdOptionFloat(arg,default_contour);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--file1") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--file1");
    file_list.push_back(arg);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--file2") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--file2");
    file_list.push_back(arg);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--file3") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--file3");
    file_list.push_back(arg);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--file4") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--file4");
    file_list.push_back(arg);
  }
  for(int c=0;c<file_list.size();c++) {
    std::cout << " file: "<<file_list[c]<<std::endl;
  }

  // Common things:
  vtkNew<vtkNamedColors> colors;
  std::array<unsigned char, 4> posColor{{240, 184, 160, 255}};  colors->SetColor("posColor", posColor.data());
  std::array<unsigned char, 4> bkg{{51, 77, 102, 255}};         colors->SetColor("BkgColor", bkg.data());

  // Create the renderer, the render window, and the interactor. The renderer
  // draws into the render window, the interactor enables mouse- and
  // keyboard-based interaction with the data within the render window.
  //
  vtkNew<vtkRenderWindow> renWin;
  vtkNew<vtkRenderWindowInteractor> iren;
  iren->SetRenderWindow(renWin);

  
  std::vector<LatticeComplexD> data(file_list.size(),&Grid);
  FieldMetaData header;

  
  int frameCount;
  if ( mpeg ) frameCount = latt_size[3];
  else        frameCount = latt_size[0] * latt_size[3];

  std::vector<FrameUpdater *> fu_list;
  for (int f=0;f<file_list.size();f++){

    // It is convenient to create an initial view of the data. The FocalPoint
    // and Position form a vector direction. Later on (ResetCamera() method)
    // this vector is used to position the camera to look at the data in
    // this direction.
    vtkNew<vtkCamera> aCamera;
    aCamera->SetViewUp(0, 0, -1);
    aCamera->SetPosition(0, -1000, 0);
    aCamera->SetFocalPoint(0, 0, 0);
    aCamera->ComputeViewPlaneNormal();
    aCamera->Azimuth(30.0);
    aCamera->Elevation(30.0);

    
    vtkNew<vtkRenderer> aRenderer;
    renWin->AddRenderer(aRenderer);
    
    double vol = data[f].Grid()->gSites();

    std::cout << "Reading "<<file_list[f]<<std::endl;
    readFile(data[f],file_list[f]);

    auto nrm    = norm2(data[f]);
    auto nrmbar = nrm/vol;
    auto rms    = sqrt(nrmbar);

    double contour = default_contour * rms; // default to 1 x RMS

    // The following reader is used to read a series of 2D slices (images)
    // that compose the volume. The slice dimensions are set, and the
    // pixel spacing. The data Endianness must also be specified. The reader
    // uses the FilePrefix in combination with the slice number to construct
    // filenames using the format FilePrefix.%d. (In this case the FilePrefix
    // is the root name of the file: quarter.)
    vtkNew<vtkImageData> imageData;
    imageData->SetDimensions(latt_size[0],latt_size[1],latt_size[2]);
    imageData->AllocateScalars(VTK_DOUBLE, 1);
    for(int xx=0;xx<latt_size[0];xx++){
      for(int yy=0;yy<latt_size[1];yy++){
	for(int zz=0;zz<latt_size[2];zz++){
	  Coordinate site({xx,yy,zz,0});
	  RealD value = real(peekSite(data[f],site));
	  imageData->SetScalarComponentFromDouble(xx,yy,zz,0,value);
    }}}

    vtkNew<isosurface> posExtractor;
    posExtractor->SetInputData(imageData);
    posExtractor->SetValue(0, contour);
  
    vtkNew<vtkStripper> posStripper;
    posStripper->SetInputConnection(posExtractor->GetOutputPort());

    vtkNew<vtkPolyDataMapper> posMapper;
    posMapper->SetInputConnection(posStripper->GetOutputPort());
    posMapper->ScalarVisibilityOff();

    vtkNew<vtkActor> pos;
    pos->SetMapper(posMapper);
    pos->GetProperty()->SetDiffuseColor(colors->GetColor3d("posColor").GetData());
    pos->GetProperty()->SetSpecular(0.3);
    pos->GetProperty()->SetSpecularPower(20);
    pos->GetProperty()->SetOpacity(0.5);

    // An isosurface, or contour value is set
    // The triangle stripper is used to create triangle strips from the
    // isosurface; these render much faster on may systems.
    vtkNew<isosurface> negExtractor;
    negExtractor->SetInputData(imageData);
    negExtractor->SetValue(0, -contour);

    vtkNew<vtkStripper> negStripper;
    negStripper->SetInputConnection(negExtractor->GetOutputPort());

    vtkNew<vtkPolyDataMapper> negMapper;
    negMapper->SetInputConnection(negStripper->GetOutputPort());
    negMapper->ScalarVisibilityOff();

    vtkNew<vtkActor> neg;
    neg->SetMapper(negMapper);
    neg->GetProperty()->SetDiffuseColor(colors->GetColor3d("Ivory").GetData());

    // An outline provides context around the data.
    vtkNew<vtkOutlineFilter> outlineData;
    outlineData->SetInputData(imageData);

    vtkNew<vtkPolyDataMapper> mapOutline;
    mapOutline->SetInputConnection(outlineData->GetOutputPort());

    vtkNew<vtkActor> outline;
    outline->SetMapper(mapOutline);
    outline->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());

    vtkNew<vtkTextActor> Text;
    Text->SetInput(file_list[f].c_str());
    Text->SetPosition2(0,0);
    Text->GetTextProperty()->SetFontSize(48);
    Text->GetTextProperty()->SetColor(colors->GetColor3d("Gold").GetData());

    vtkNew<vtkTextActor> TextT;
    TextT->SetInput("T=0");
    TextT->SetPosition(0,.9*1025);
    TextT->GetTextProperty()->SetFontSize(48);
    TextT->GetTextProperty()->SetColor(colors->GetColor3d("Gold").GetData());
    
  
    // Actors are added to the renderer. An initial camera view is created.
    // The Dolly() method moves the camera towards the FocalPoint,
    // thereby enlarging the image.
    aRenderer->AddActor(Text);
    aRenderer->AddActor(TextT);
    aRenderer->AddActor(outline);
    aRenderer->AddActor(pos);
    aRenderer->AddActor(neg);

    // Sign up to receive TimerEvent
    vtkNew<FrameUpdater> fu;
    fu->imageData = imageData;
    fu->grid_data = &data[f];
    fu->text      = TextT;
    fu->maxCount = frameCount;
    fu->posExtractor = posExtractor;
    fu->negExtractor = negExtractor;
    fu->rms = rms;
      
    iren->AddObserver(vtkCommand::TimerEvent, fu);

    aRenderer->SetActiveCamera(aCamera);
    aRenderer->ResetCamera();
    aRenderer->SetBackground(colors->GetColor3d("BkgColor").GetData());
    aCamera->Dolly(1.0);

    double nf = file_list.size();
    std::cout << " Adding renderer " <<f<<" of "<<nf<<std::endl;
    aRenderer->SetViewport((1.0/nf)*f, 0.0,(1.0/nf)*(f+1) , 1.0);

    // Note that when camera movement occurs (as it does in the Dolly()
    // method), the clipping planes often need adjusting. Clipping planes
    // consist of two planes: near and far along the view direction. The
    // near plane clips out objects in front of the plane; the far plane
    // clips out objects behind the plane. This way only what is drawn
    // between the planes is actually rendered.
    aRenderer->ResetCameraClippingRange();
    
    fu_list.push_back(fu);
  }


  // Set a background color for the renderer and set the size of the
  // render window (expressed in pixels).
  // Initialize the event loop and then start it.
  renWin->SetSize(1024*file_list.size(), 1024);
  renWin->SetWindowName("FieldDensity");
  renWin->Render();

  iren->Initialize();

  if ( mpeg ) {
#ifdef MPEG
    vtkWindowToImageFilter *imageFilter = vtkWindowToImageFilter::New();
    imageFilter->SetInput( renWin );
    imageFilter->SetInputBufferTypeToRGB();
    
    vtkFFMPEGWriter *writer = vtkFFMPEGWriter::New();
    writer->SetFileName("movie.avi");
    writer->SetRate(framerate);
    writer->SetInputConnection(imageFilter->GetOutputPort());
    writer->Start();

    for(int i=0;i<fu_list[0]->maxCount;i++){
      for(int f=0;f<fu_list.size();f++){
	fu_list[f]->Execute(iren,vtkCommand::TimerEvent,nullptr);
      }
      imageFilter->Modified();
      writer->Write();
    }
    writer->End();
    writer->Delete();
#else
    assert(-1 && "MPEG support not compiled");
#endif
  } else { 
  
    // Add control of contour threshold
    // Create a slider widget
    vtkSmartPointer<vtkSliderRepresentation2D> sliderRep = vtkSmartPointer<vtkSliderRepresentation2D>::New();
    sliderRep->SetMinimumValue(0.1);
    sliderRep->SetMaximumValue(5.0);
    sliderRep->SetValue(1.0);
    sliderRep->SetTitleText("Fraction RMS");
    // Set color properties:

    // Change the color of the knob that slides
    //  sliderRep->GetSliderProperty()->SetColor(colors->GetColor3d("Green").GetData());
    sliderRep->GetTitleProperty()->SetColor(colors->GetColor3d("AliceBlue").GetData());
    sliderRep->GetLabelProperty()->SetColor(colors->GetColor3d("AliceBlue").GetData());
    sliderRep->GetSelectedProperty()->SetColor(colors->GetColor3d("DeepPink").GetData());

    // Change the color of the bar
    sliderRep->GetTubeProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());
    sliderRep->GetCapProperty()->SetColor(colors->GetColor3d("Yellow").GetData());
    sliderRep->SetSliderLength(0.05);
    sliderRep->SetSliderWidth(0.025);
    sliderRep->SetEndCapLength(0.02);

    double nf = file_list.size();
    sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
    sliderRep->GetPoint1Coordinate()->SetValue(0.1, 0.1);
    sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
    sliderRep->GetPoint2Coordinate()->SetValue(0.9/nf, 0.1);
  
    vtkSmartPointer<vtkSliderWidget> sliderWidget = vtkSmartPointer<vtkSliderWidget>::New();
    sliderWidget->SetInteractor(iren);
    sliderWidget->SetRepresentation(sliderRep);
    sliderWidget->SetAnimationModeToAnimate();
    sliderWidget->EnabledOn();
  
    // Create the slider callback
    vtkSmartPointer<SliderCallback> slidercallback = vtkSmartPointer<SliderCallback>::New();
    slidercallback->fu_list = fu_list;
    sliderWidget->AddObserver(vtkCommand::InteractionEvent, slidercallback);

    int timerId = iren->CreateRepeatingTimer(1000/framerate);
    std::cout << "timerId: " << timerId << std::endl;

    // Start the interaction and timer
    iren->Start();
  }

  Grid_finalize();

  return EXIT_SUCCESS;
}
