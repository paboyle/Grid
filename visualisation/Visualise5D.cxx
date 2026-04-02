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
    ffile=0;
    TimerCount = 0;
    xoff       = 0;
    t          = 0;
    imageData = nullptr;
    timerId = 0;
    maxCount = -1;
    old_file=-1;
  }
  
  static FrameUpdater* New()
  {
    FrameUpdater* cb = new FrameUpdater;
    cb->TimerCount = 0;
    return cb;
  }

  //
  // Must map a x,y,z + frame index into
  // i)  a d-dimensional site Coordinate
  // ii) a file name
  // Need a:
  //     loop_ranges
  //     sum_ranges
  //     loop_vol  -- map loop_idx -> loop_coor
  //     sum_vol   -- map sum_idx -> sum_coor with Lexicographic
  //
  /*
   * Just set this up
   */
  int old_file ; // Cache, avoid reread
  
  Coordinate          latt;
  Coordinate          xyz_dims    ; // List lattice dimensions corresponding to xyz_dims displayed 
  Coordinate          xyz_ranges  ; // 3-vector
  Coordinate          g_xyz_ranges; // Nd-vector
  uint64_t            xyz_vol     ;
  
  Coordinate          loop_dims;    // List lattice dimensions put into movie time
  Coordinate          loop_ranges;  // movie time ranges
  uint64_t            loop_vol;

  Coordinate          sum_dims;      // List lattice dimensions summed
  Coordinate          sum_ranges;    // summation ranges
  uint64_t            sum_vol;

  Coordinate          slice_dims;      // List slice dimensions 

  Coordinate          Slice;
  
  std::vector<std::string> files;        // file list that is looped over

  int Nd;
  GridBase *grid;
  Grid::LatticeComplexD *grid_data;

  void SetGrid(GridBase *_grid)
  {
    grid = _grid;
    Nd=grid->Nd();
    latt = grid->GlobalDimensions();
    grid_data = new Grid::LatticeComplexD(grid);
  }

  void SetFiles(std::vector<std::string> list)      { files = list; old_file = -1; }
  void SetSlice(Coordinate _Slice)             { Slice = _Slice;} // Offset / skew for lattice coords
  void SetSumDimensions (Coordinate _SumDims ) {
    sum_ranges=Coordinate(Nd);
    sum_dims = _SumDims; // 1 hot for dimensions summed
    sum_vol = 1;
    for(int d=0;d<sum_dims.size();d++){
      if ( sum_dims[d] == 1 ) sum_ranges[d] = latt[d];
      else                    sum_ranges[d] = 1;
      sum_vol*=sum_ranges[d];
    }
  }
  void SetLoopDimensions(Coordinate _LoopDims) {
    loop_ranges=Coordinate(Nd);
    loop_dims= _LoopDims;
    loop_vol = 1;
    for(int d=0;d<loop_dims.size();d++){
      if ( loop_dims[d] == 1 ) loop_ranges[d] = latt[d];
      else                     loop_ranges[d] = 1;
      loop_vol*=loop_ranges[d];
    }
  } // 
  void SetDisplayDimensions(Coordinate _xyz_dims   ) {
    g_xyz_ranges=Coordinate(Nd);
    xyz_ranges=Coordinate(3);
    xyz_dims = _xyz_dims;
    xyz_vol  = 1;
    for(int d=0;d<3;d++){
      xyz_ranges[d] = latt[xyz_dims[d]];
      xyz_vol *= xyz_ranges[d];
    }
    // Find dim extents for grid
    int dd=0;
    for(int d=0;d<Nd;d++){
      g_xyz_ranges[d] = 1;
      for(int dd=0;dd<3;dd++) {
	if ( xyz_dims[dd]==d ) {
	  g_xyz_ranges[d] = latt[d];
	}
      }
    }
  }
  void SetSliceDimensions(void) {
    Coordinate _slice_dims;
    for ( int d=0;d<Nd;d++){
      int is_slice = 1;
      if(g_xyz_ranges[d]>1) is_slice = 0;
      if(loop_dims[d]) is_slice = 0;
      if(sum_dims[d] ) is_slice = 0;
      if(is_slice) _slice_dims.push_back(d);
    }
    slice_dims = _slice_dims;
    std::cout << " Setting Slice Dimensions to "<<slice_dims<<std::endl;
  }
  virtual void Execute(vtkObject* caller, unsigned long eventId,void* vtkNotUsed(callData))
  {
    const int max=256;
    char text_string[max];

    auto latt_size = grid->GlobalDimensions();
    
    if ( vtkCommand::KeyPressEvent == eventId ) {

      vtkRenderWindowInteractor* iren = static_cast<vtkRenderWindowInteractor*>(caller);
      std::string key = iren->GetKeySym();
      std::cout << "Pressed: " << key << std::endl;
      if (slice_dims.size()>0) {
	int vert = slice_dims[slice_dims.size()-1];
	int horz = slice_dims[0];
	if ( key == "Up" ) {
	  Slice[vert] = (Slice[vert]+1)%latt[vert];
	}
	if ( key == "Down" ) {
	  Slice[vert] = (Slice[vert]+latt[vert]-1)%latt[vert];
	}
	if ( key == "Right" ) {
	  Slice[horz] = (Slice[horz]+1)%latt[horz];
	}
	if ( key == "Left" ) {
	  Slice[horz] = (Slice[horz]+latt[horz]-1)%latt[horz];
	}
      }
      if ( key == "greater" ) {
	ffile = (ffile + 1) % files.size();
      }
      if ( key == "less" ) {
	ffile = (ffile - 1 + files.size()) % files.size();
      }
      std::cout <<"Slice " <<Slice <<std::endl;
      std::cout <<"File  " <<ffile <<std::endl;
    }
      
    // Make a new frame for frame index TimerCount
    if ( vtkCommand::TimerEvent == eventId || vtkCommand::KeyPressEvent == eventId)
      {
	int file     = ((this->TimerCount / loop_vol) + ffile )%files.size();

	if ( file != old_file ) {
	  readFile(*grid_data,files[file]);
	  old_file = file;
	}

	RealD max, min, max_abs,min_abs;
	Coordinate max_site;
	Coordinate min_site;
	Coordinate max_abs_site;
	Coordinate min_abs_site;
	for(int idx=0;idx<grid->gSites();idx++){
	  Coordinate site;
	  Lexicographic::CoorFromIndex (site,idx,latt);
	  RealD val=real(peekSite(*grid_data,site));
	  if (idx==0){
	    max = min = val;
	    max_abs = min_abs = fabs(val);
	    max_site = site;
	    min_site = site;
	    min_abs_site = site;
	    max_abs_site = site;
	  } else {
	    if ( val > max ) {
	      max=val;
	      max_site = site;
	    }
	    if ( fabs(val) > max_abs ) {
	      max_abs=fabs(val);
	      max_abs_site = site;
	    }
	    if ( val < min ) {
	      min=val;
	      min_site = site;
	    }	    
	    if ( fabs(val) < min_abs ) {
	      min_abs=fabs(val);
	      min_abs_site = site;
	    }	    
	  }
	}
	std::cout << " abs_max "<<max_abs<<" at " << max_abs_site<<std::endl;
	std::cout << " abs_min "<<min_abs<<" at " << min_abs_site<<std::endl;
	std::cout << " max "<<max<<" at " << max_site<<std::endl;
	std::cout << " min "<<min<<" at " << min_site<<std::endl;
	// Looped dimensions, map index to coordinate
	int loop_idx = this->TimerCount % loop_vol;
	Coordinate loop_coor;
	Lexicographic::CoorFromIndex (loop_coor,loop_idx,loop_ranges);
	
	// Loop over xyz sites
	Coordinate xyz_coor(3);
	Coordinate g_xyz_coor(Nd);
	Coordinate sum_coor(Nd);
	for(uint64_t xyz = 0 ; xyz< xyz_vol; xyz++){
	  Lexicographic::CoorFromIndex (xyz_coor,xyz,xyz_ranges);
	  Lexicographic::CoorFromIndex (g_xyz_coor,xyz,g_xyz_ranges);
	  
	  RealD sum_value = 0.0;
	  for(uint64_t sum_idx = 0 ; sum_idx< sum_vol; sum_idx++){
	    Lexicographic::CoorFromIndex (sum_coor,sum_idx,sum_ranges);
	    Coordinate site(Nd);
	    for(int d=0;d<Nd;d++){
	      site[d] = (sum_coor[d] + loop_coor[d] + g_xyz_coor[d] + Slice[d])%latt[d];
	    }
	    sum_value+= real(peekSite(*grid_data,site));
	    if(xyz==0) std::cout << "sum "<<sum_idx<<" "<<sum_value<<std::endl;
	  }
	  imageData->SetScalarComponentFromDouble(xyz_coor[0],xyz_coor[1],xyz_coor[2],0,sum_value);
	}
	imageData->Modified();

	std::stringstream ss;
	ss<< files[file] <<"\nSlice "<<Slice << "\nLoop  " <<loop_coor<<"\nSummed "<<sum_dims;

	text->SetInput(ss.str().c_str());
	

	vtkRenderWindowInteractor* iren = dynamic_cast<vtkRenderWindowInteractor*>(caller);
	iren->GetRenderWindow()->Render();
	  
      }

    if ( vtkCommand::TimerEvent == eventId ) {
      ++this->TimerCount;
      std::cout << " This was a timer event count "<<this->TimerCount << std::endl;
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
  int ffile;
  int xoff;
  int t;
public:
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
	fu->posExtractor->SetValue(0,  SliderCallback::contour*fu->rms);
	fu->negExtractor->SetValue(0, -SliderCallback::contour*fu->rms);
	fu->posExtractor->Modified();
	fu->negExtractor->Modified();
    }    
public:
  static double contour;
  FrameUpdater * fu;
};

FrameUpdater * KBfu;
void KeypressCallbackFunction(vtkObject* caller, long unsigned int eventId,
                              void* clientData, void* callData)
{
  std::cout << "Keypress callback" << std::endl;
  vtkRenderWindowInteractor* iren = static_cast<vtkRenderWindowInteractor*>(caller);
  std::cout << "Pressed: " << iren->GetKeySym() << std::endl;
  //  imageData->Modified();
}

double SliderCallback::contour;

int main(int argc, char* argv[])
{
  using namespace Grid;

  Grid_init(&argc, &argv);
  GridLogLayout();

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(latt_size.size(), vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  GridCartesian    Grid(latt_size, simd_layout, mpi_layout);

  double default_contour = 1.0;
  std::string arg;
 
  std::cout << argc << " command Line arguments "<<std::endl;
  for(int c=0;c<argc;c++) {
    std::cout << " - "<<argv[c]<<std::endl;
  }

  std::vector<std::string> file_list({
      "file1",
      "file2",
      "file3",
      "file4",
      "file5",
      "file6",
      "file7",
      "file8"
    });

  if( GridCmdOptionExists(argv,argv+argc,"--files") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--files");
    GridCmdOptionCSL(arg, file_list);
  }

#ifdef MPEG
  if( GridCmdOptionExists(argv,argv+argc,"--mpeg") ){
    mpeg = 1;
  }
#endif

  if( GridCmdOptionExists(argv,argv+argc,"--fps") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--fps");
    GridCmdOptionInt(arg,framerate);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--isosurface") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--isosurface");
    GridCmdOptionFloat(arg,default_contour);
  }
  for(int c=0;c<file_list.size();c++) {
    std::cout << " file: "<<file_list[c]<<std::endl;
  }

  int NoTime = 0;
  int Nd; Nd = Grid.Nd();
  Coordinate    Slice(Nd,0);
  Coordinate  SumDims(Nd,0);
  Coordinate LoopDims(Nd,0);
  Coordinate  XYZDims({0,1,2});
  if( GridCmdOptionExists(argv,argv+argc,"--slice") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--slice");
    GridCmdOptionIntVector(arg,Slice);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--sum") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--sum");
    GridCmdOptionIntVector(arg,SumDims);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--loop") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--loop");
    GridCmdOptionIntVector(arg,LoopDims);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--xyz") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--xyz");
    GridCmdOptionIntVector(arg,XYZDims);
    std::cout << "xyz : "<<XYZDims<<std::endl;
  }
  if( GridCmdOptionExists(argv,argv+argc,"--notime") ){
    NoTime = 1;
    std::cout << "Suppressing time loop"<<std::endl;
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

  
  //  std::vector<LatticeComplexD> data(file_list.size(),&Grid);
  //  FieldMetaData header;

  

  int frameCount = file_list.size();
  for(int d=0;d<Grid.Nd();d++) {
    if ( LoopDims[d] ) frameCount*= latt_size[d];

  }
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
    
  double vol = Grid.gSites();
  std::cout << "Reading "<<file_list[0]<<std::endl;
  double nrm, nrmbar,rms, contour;
  {
    LatticeComplexD data(&Grid);
    readFile(data,file_list[0]);
    nrm    = norm2(data);
  }
  nrmbar = nrm/vol;
  rms    = sqrt(nrmbar);
  
  contour = default_contour * rms; // default to 1 x RMS

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
	RealD value = 0;
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
  //  Text->SetInput(file_list[f].c_str());
  Text->SetPosition2(0,0);
  Text->GetTextProperty()->SetFontSize(24);
  Text->GetTextProperty()->SetColor(colors->GetColor3d("Gold").GetData());

  vtkNew<vtkTextActor> TextT;
  TextT->SetInput("T=0");
  TextT->SetPosition(0,.7*1025);
  TextT->GetTextProperty()->SetFontSize(24);
  TextT->GetTextProperty()->SetColor(colors->GetColor3d("Gold").GetData());
    
  
  // Actors are added to the renderer. An initial camera view is created.
  // The Dolly() method moves the camera towards the FocalPoint,
  // thereby enlarging the image.
  //    aRenderer->AddActor(Text);
  aRenderer->AddActor(TextT);
  aRenderer->AddActor(outline);
  aRenderer->AddActor(pos);
  aRenderer->AddActor(neg);

  // Sign up to receive TimerEvent
  vtkNew<FrameUpdater> fu;

  fu->SetGrid(&Grid);
  fu->SetFiles(file_list);
  fu->SetSlice(Slice);
  fu->SetSumDimensions (SumDims);
  fu->SetLoopDimensions(LoopDims);
  fu->SetDisplayDimensions(XYZDims);
  fu->SetSliceDimensions();
  

  fu->imageData = imageData;
  //    fu->grid_data = &data[f];
  fu->text      = TextT;
  fu->maxCount = frameCount;
  fu->posExtractor = posExtractor;
  fu->negExtractor = negExtractor;
  fu->rms = rms;
      
  iren->AddObserver(vtkCommand::TimerEvent, fu);
  iren->AddObserver(vtkCommand::KeyPressEvent, fu);

  aRenderer->SetActiveCamera(aCamera);
  aRenderer->ResetCamera();
  aRenderer->SetBackground(colors->GetColor3d("BkgColor").GetData());
  aCamera->Dolly(1.0);

  //    double nf = file_list.size();
  //    std::cout << " Adding renderer " <<f<<" of "<<nf<<std::endl;
  aRenderer->SetViewport(0.0, 0.0,1.0 , 1.0);

  // Note that when camera movement occurs (as it does in the Dolly()
  // method), the clipping planes often need adjusting. Clipping planes
  // consist of two planes: near and far along the view direction. The
  // near plane clips out objects in front of the plane; the far plane
  // clips out objects behind the plane. This way only what is drawn
  // between the planes is actually rendered.
  aRenderer->ResetCameraClippingRange();
    

  // Set a background color for the renderer and set the size of the
  // render window (expressed in pixels).
  // Initialize the event loop and then start it.
  renWin->SetSize(1024, 1024);
  renWin->SetWindowName("FieldDensity");
  renWin->Render();

  // Take a pointer to the FrameUpdater for keypress mgt.
  //  KBfu = fu;
  //  vtkNew<vtkCallbackCommand> keypressCallback;
  //  keypressCallback->SetCallback(KeypressCallbackFunction);
  //  iren->AddObserver(vtkCommand::KeyPressEvent,keypressCallback);

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

    for(int i=0;i<fu->maxCount;i++){
	fu->Execute(iren,vtkCommand::TimerEvent,nullptr);
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
    sliderRep->SetMinimumValue(0.0);
    sliderRep->SetMaximumValue(10.0);
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

    sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
    sliderRep->GetPoint1Coordinate()->SetValue(0.1, 0.1);
    sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
    sliderRep->GetPoint2Coordinate()->SetValue(0.9, 0.1);
  
    vtkSmartPointer<vtkSliderWidget> sliderWidget = vtkSmartPointer<vtkSliderWidget>::New();
    sliderWidget->SetInteractor(iren);
    sliderWidget->SetRepresentation(sliderRep);
    sliderWidget->SetAnimationModeToAnimate();
    sliderWidget->EnabledOn();
    
    // Create the slider callback
    vtkSmartPointer<SliderCallback> slidercallback = vtkSmartPointer<SliderCallback>::New();
    slidercallback->fu = fu;
    sliderWidget->AddObserver(vtkCommand::InteractionEvent, slidercallback);

    if ( NoTime==0 ) {
      int timerId = iren->CreateRepeatingTimer(10000/framerate);
      std::cout << "timerId "<<timerId<<std::endl;
    }

    // Start the interaction and timer
    iren->Start();
  }

  Grid_finalize();

  return EXIT_SUCCESS;
}
