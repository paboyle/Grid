// ForceAnalysis.cxx
//
// Reads a sequence of force snapshot files (LatticeComplexD, real part = force magnitude)
// and produces two outputs:
//
//   1. Tab-separated timeseries to stdout:
//        idx  Gauge_lat_rms  Gauge_lat_hot  Gauge_smr_rms  ...
//      where _rms is the lattice RMS and _hot is the value at --hotsite.
//
//   2. PNG images (one per force component per snapshot) rendered via VTK
//      as isosurfaces of the force density, using the same pipeline as
//      Visualise5D.  Images are written to --pngdir/<label>_<idx>.png.
//      These can be read back by Claude to interpret spatial structure.
//
// Usage:
//   ForceAnalysis --grid 32.32.32.32 --mpi 1.1.1.1
//                 --snapdir /path/to/snapshots
//                 --first 0 --last 1920 --step 10
//                 --hotsite x.y.z.t
//                 --pngdir /path/to/output/pngs
//                 --isosurface 1.0        (contour in units of field RMS)
//                 --fixediso 0.05         (fixed absolute contour, overrides --isosurface)
//                 --slice t               (which dimension to fix for 3D display, default: t)
//                 --sliceval 2            (value of that dimension, default: 0)
//
// Dimension order on the 32^4 lattice: x=0 y=1 z=2 t=3

#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkCamera.h>
#include <vtkImageActor.h>
#include <vtkImageMapper3D.h>
#include <vtkImageData.h>
#include <vtkImageMapToColors.h>
#include <vtkLookupTable.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkOutlineFilter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkProperty.h>
#include <vtkProperty2D.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkStripper.h>
#include <vtkCallbackCommand.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>

#define USE_FLYING_EDGES
#ifdef USE_FLYING_EDGES
#include <vtkFlyingEdges3D.h>
typedef vtkFlyingEdges3D isosurface;
#else
#include <vtkMarchingCubes.h>
typedef vtkMarchingCubes isosurface;
#endif

#include <Grid/Grid.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <memory>
#include <sys/stat.h>

using namespace Grid;

// ─── I/O ─────────────────────────────────────────────────────────────────────
template <class T>
bool tryReadFile(T& out, const std::string& fname)
{
    std::ifstream test(fname);
    if (!test.good()) return false;
    test.close();
    emptyUserRecord record;
    ScidacReader RD;
    RD.open(fname);
    RD.readScidacFieldRecord(out, record);
    RD.close();
    return true;
}

// ─── Fill a 3D vtkImageData slice from a 4D lattice field ────────────────────
// Sums over the sliced dimension at sliceval, displays the remaining 3 dims.
void fillImageData(vtkImageData* img,
                   LatticeComplexD& field,
                   const Coordinate& latt_size,
                   int slice_dim, int sliceval)
{
    // Display dims = all dims except slice_dim, in order
    std::vector<int> disp;
    for (int d = 0; d < 4; d++) if (d != slice_dim) disp.push_back(d);

    int Nx = latt_size[disp[0]];
    int Ny = latt_size[disp[1]];
    int Nz = latt_size[disp[2]];

    for (int ix = 0; ix < Nx; ix++)
    for (int iy = 0; iy < Ny; iy++)
    for (int iz = 0; iz < Nz; iz++) {
        Coordinate site(4);
        site[disp[0]]  = ix;
        site[disp[1]]  = iy;
        site[disp[2]]  = iz;
        site[slice_dim] = sliceval;
        RealD val = real(peekSite(field, site));
        img->SetScalarComponentFromDouble(ix, iy, iz, 0, val);
    }
    img->Modified();
}

// ─── 2D heatmap: persistent context ───────────────────────────────────────────
// Renders a fixed (dim1=v1, dim2=v2) slice of the 4D force field as a
// diverging blue→white→red colour map, with a fixed symmetric colour scale
// so brightness directly encodes force magnitude across all frames.
// A white cross-hair marks the hotsite projection onto the slice.
struct HeatmapCtx {
    // image pipeline
    vtkNew<vtkImageData>         img;
    vtkNew<vtkLookupTable>       lut;
    vtkNew<vtkImageMapToColors>  colorMap;
    vtkNew<vtkImageActor>        imgActor;
    // colour scale legend (text, avoids needing RenderingAnnotation module)
    vtkNew<vtkTextActor>         cbar;
    // hotsite cross-hair (2D overlay actors)
    vtkNew<vtkPolyData>          crossPD;
    vtkNew<vtkPoints>            crossPts;
    vtkNew<vtkCellArray>         crossLines;
    vtkNew<vtkActor2D>           crossActor;
    // title
    vtkNew<vtkTextActor>         titleAct;
    // renderer / window
    vtkNew<vtkRenderer>          ren;
    vtkNew<vtkRenderWindow>      renWin;
    vtkNew<vtkWindowToImageFilter> w2i;
    vtkNew<vtkPNGWriter>         writer;

    int Nx = 0, Ny = 0;           // display dimensions of the slice
    double scale = 0.07;          // colour range: [-scale, +scale]
    int hotX = -1, hotY = -1;     // hotsite projection onto (Nx,Ny) plane
    // pixel coords of the image origin in the render window
    int imgOffX = 60, imgOffY = 40;
    int imgW = 0, imgH = 0;       // rendered pixel size of each lattice cell

    void init(int nx, int ny, double sc, int hx, int hy)
    {
        Nx = nx; Ny = ny; scale = sc;
        hotX = hx; hotY = hy;

        const int WIN_W = 900, WIN_H = 700;
        // Make cells square and as large as possible within the central area
        int cellW = (WIN_W - 160) / Nx;
        int cellH = (WIN_H - 120) / Ny;
        imgW = std::min(cellW, cellH);
        imgH = imgW;
        imgOffX = (WIN_W - Nx * imgW) / 2;
        imgOffY = 60;

        // --- Image data (scalar field, one component) ---
        img->SetDimensions(Nx, Ny, 1);
        img->SetSpacing(imgW, imgH, 1);
        img->SetOrigin(imgOffX, imgOffY, 0);
        img->AllocateScalars(VTK_DOUBLE, 1);

        // --- Diverging LUT: blue(-scale) → white(0) → red(+scale) ---
        lut->SetNumberOfTableValues(512);
        lut->SetRange(-scale, scale);
        lut->SetNanColor(0.2, 0.2, 0.2, 1.0);
        for (int i = 0; i < 512; ++i) {
            double t = i / 511.0;   // 0=blue, 0.5=white, 1=red
            double r = (t > 0.5) ? 1.0 : 2.0 * t;
            double g = (t < 0.5) ? 2.0 * t : 2.0 * (1.0 - t);
            double b = (t < 0.5) ? 1.0 : 2.0 * (1.0 - t);
            lut->SetTableValue(i, r, g, b, 1.0);
        }
        lut->Build();

        // --- Colour map pipeline ---
        colorMap->SetInputData(img);
        colorMap->SetLookupTable(lut);
        colorMap->Update();

        imgActor->GetMapper()->SetInputConnection(colorMap->GetOutputPort());

        // --- Colour scale legend (text) ---
        {
            std::ostringstream ss;
            ss << std::scientific << std::setprecision(2)
               << "blue=-" << sc << "  white=0  red=+" << sc;
            cbar->SetInput(ss.str().c_str());
        }
        cbar->GetTextProperty()->SetFontFamilyToCourier();
        cbar->GetTextProperty()->SetFontSize(13);
        cbar->GetTextProperty()->SetColor(0.9, 0.9, 0.9);
        cbar->SetDisplayPosition(10, 10);

        // --- Cross-hair at hotsite (2D display coords) ---
        if (hotX >= 0 && hotY >= 0) {
            double cx = imgOffX + (hotX + 0.5) * imgW;
            double cy = imgOffY + (hotY + 0.5) * imgH;
            double arm = imgW * 0.8;
            crossPts->InsertNextPoint(cx - arm, cy, 0);
            crossPts->InsertNextPoint(cx + arm, cy, 0);
            crossPts->InsertNextPoint(cx, cy - arm, 0);
            crossPts->InsertNextPoint(cx, cy + arm, 0);
            vtkIdType seg0[2] = {0, 1};
            vtkIdType seg1[2] = {2, 3};
            crossLines->InsertNextCell(2, seg0);
            crossLines->InsertNextCell(2, seg1);
            crossPD->SetPoints(crossPts);
            crossPD->SetLines(crossLines);
            vtkNew<vtkPolyDataMapper2D> crossMap;
            crossMap->SetInputData(crossPD);
            crossActor->SetMapper(crossMap);
            crossActor->GetProperty()->SetColor(1, 1, 1);
            crossActor->GetProperty()->SetLineWidth(2.0);
        }

        // --- Title ---
        titleAct->GetTextProperty()->SetFontFamilyToCourier();
        titleAct->GetTextProperty()->SetFontSize(16);
        titleAct->GetTextProperty()->SetColor(1, 1, 0);
        titleAct->SetDisplayPosition(10, WIN_H - 30);

        // --- Renderer (2D parallel projection so image fills correctly) ---
        ren->SetBackground(0.08, 0.08, 0.12);
        ren->AddActor(imgActor);
        ren->AddActor2D(cbar);
        ren->AddActor2D(crossActor);
        ren->AddActor2D(titleAct);
        ren->GetActiveCamera()->ParallelProjectionOn();
        // Set up camera to look straight down at the image plane
        ren->GetActiveCamera()->SetPosition(WIN_W/2.0, WIN_H/2.0, 1000);
        ren->GetActiveCamera()->SetFocalPoint(WIN_W/2.0, WIN_H/2.0, 0);
        ren->GetActiveCamera()->SetViewUp(0, 1, 0);
        ren->GetActiveCamera()->SetParallelScale(WIN_H / 2.0);
        ren->ResetCameraClippingRange();

        renWin->AddRenderer(ren);
        renWin->SetSize(WIN_W, WIN_H);
        renWin->SetOffScreenRendering(1);
        renWin->SetMultiSamples(0);

        w2i->SetInput(renWin);
        w2i->SetInputBufferTypeToRGB();
        w2i->ReadFrontBufferOff();
    }
};

void renderHeatmap(HeatmapCtx& ctx,
                   LatticeComplexD& field,
                   const Coordinate& latt_size,
                   int dim1, int val1,   // first fixed dimension
                   int dim2, int val2,   // second fixed dimension
                   const std::string& title,
                   const std::string& outpath)
{
    // Display dimensions: the two dims that are NOT fixed
    std::vector<int> disp;
    for (int d = 0; d < 4; d++)
        if (d != dim1 && d != dim2) disp.push_back(d);

    int Nx = latt_size[disp[0]];
    int Ny = latt_size[disp[1]];

    // Fill image data
    for (int ix = 0; ix < Nx; ix++) {
        for (int iy = 0; iy < Ny; iy++) {
            Coordinate site(4);
            site[disp[0]] = ix;
            site[disp[1]] = iy;
            site[dim1]    = val1;
            site[dim2]    = val2;
            RealD val = real(TensorRemove(peekSite(field, site)));
            ctx.img->SetScalarComponentFromDouble(ix, iy, 0, 0, val);
        }
    }
    ctx.img->Modified();
    ctx.colorMap->Update();

    ctx.titleAct->SetInput(title.c_str());

    ctx.renWin->Render();
    ctx.w2i->Modified();
    ctx.w2i->Update();

    ctx.writer->SetFileName(outpath.c_str());
    ctx.writer->SetInputConnection(ctx.w2i->GetOutputPort());
    ctx.writer->Write();
}

// ─── Persistent rendering context (created once, reused every frame) ──────────
// Avoids Metal GPU context exhaustion on macOS when rendering hundreds of frames.
struct RenderCtx {
    vtkNew<vtkNamedColors>       colors;
    vtkNew<vtkImageData>         imageData;
    vtkNew<isosurface>           posEx, negEx;
    vtkNew<vtkStripper>          posSt, negSt;
    vtkNew<vtkPolyDataMapper>    posMap, negMap, outMap;
    vtkNew<vtkActor>             posAct, negAct, outAct;
    vtkNew<vtkOutlineFilter>     outF;
    vtkNew<vtkTextActor>         label;
    vtkNew<vtkRenderer>          ren;
    vtkNew<vtkCamera>            cam;
    vtkNew<vtkRenderWindow>      renWin;
    vtkNew<vtkWindowToImageFilter> w2i;
    vtkNew<vtkPNGWriter>         writer;

    void init(int Nx, int Ny, int Nz)
    {
        std::array<unsigned char,4> posColor{{240,184,160,255}};
        colors->SetColor("posColor", posColor.data());
        std::array<unsigned char,4> bkg{{51,77,102,255}};
        colors->SetColor("BkgColor", bkg.data());

        imageData->SetDimensions(Nx, Ny, Nz);
        imageData->AllocateScalars(VTK_DOUBLE, 1);

        posEx->SetInputData(imageData); posEx->SetValue(0,  1.0);
        posSt->SetInputConnection(posEx->GetOutputPort());
        posMap->SetInputConnection(posSt->GetOutputPort());
        posMap->ScalarVisibilityOff();
        posAct->SetMapper(posMap);
        posAct->GetProperty()->SetDiffuseColor(colors->GetColor3d("posColor").GetData());
        posAct->GetProperty()->SetSpecular(0.3);
        posAct->GetProperty()->SetSpecularPower(20);
        posAct->GetProperty()->SetOpacity(0.6);

        negEx->SetInputData(imageData); negEx->SetValue(0, -1.0);
        negSt->SetInputConnection(negEx->GetOutputPort());
        negMap->SetInputConnection(negSt->GetOutputPort());
        negMap->ScalarVisibilityOff();
        negAct->SetMapper(negMap);
        negAct->GetProperty()->SetDiffuseColor(colors->GetColor3d("Ivory").GetData());
        negAct->GetProperty()->SetOpacity(0.6);

        outF->SetInputData(imageData);
        outMap->SetInputConnection(outF->GetOutputPort());
        outAct->SetMapper(outMap);
        outAct->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());

        label->SetPosition(10, 10);
        label->GetTextProperty()->SetFontFamilyToCourier();
        label->GetTextProperty()->SetFontSize(18);
        label->GetTextProperty()->SetColor(colors->GetColor3d("Gold").GetData());

        ren->AddActor(posAct);
        ren->AddActor(negAct);
        ren->AddActor(outAct);
        ren->AddActor2D(label);
        ren->SetBackground(colors->GetColor3d("BkgColor").GetData());

        cam->SetViewUp(0,0,-1);
        cam->SetPosition(0,-1000,0);
        cam->SetFocalPoint(0,0,0);
        cam->ComputeViewPlaneNormal();
        cam->Azimuth(30.0);
        cam->Elevation(30.0);
        ren->SetActiveCamera(cam);

        renWin->AddRenderer(ren);
        renWin->SetSize(800, 600);
        renWin->SetOffScreenRendering(1);
        renWin->SetMultiSamples(0);

        w2i->SetInput(renWin);
        w2i->SetInputBufferTypeToRGB();
        w2i->ReadFrontBufferOff();
    }
};

// ─── Render one force field snapshot to a PNG (reuses existing RenderCtx) ─────
void renderPNG(RenderCtx& ctx,
               LatticeComplexD& field,
               const Coordinate& latt_size,
               int slice_dim, int sliceval,
               double contour,
               const std::string& title,
               const std::string& outpath)
{
    // Update image data
    fillImageData(ctx.imageData, field, latt_size, slice_dim, sliceval);

    // Update isosurface levels
    ctx.posEx->SetValue(0,  contour);
    ctx.negEx->SetValue(0, -contour);

    // Update label
    ctx.label->SetInput(title.c_str());

    // Reset camera to fit the (possibly new) data bounds
    ctx.ren->ResetCamera();
    ctx.cam->Dolly(1.2);
    ctx.ren->ResetCameraClippingRange();

    ctx.renWin->Render();

    ctx.w2i->Modified();
    ctx.w2i->Update();

    ctx.writer->SetFileName(outpath.c_str());
    ctx.writer->SetInputConnection(ctx.w2i->GetOutputPort());
    ctx.writer->Write();
}

// ─── main ─────────────────────────────────────────────────────────────────────
int main(int argc, char* argv[])
{
    Grid_init(&argc, &argv);
    GridLogMessage.Active(0);
    GridLogIterative.Active(0);
    GridLogDebug.Active(0);
    GridLogPerformance.Active(0);
    GridLogComms.Active(0);
    GridLogDslash.Active(0);
    GridLogMemory.Active(0);

    // ── CLI ──────────────────────────────────────────────────────────────────
    std::string snapdir = ".";
    std::string pngdir  = "";
    int first = 0, last = 1920, step = 1;
    int slice_dim = 3, sliceval = 0;   // default: fix t=0, display xyz
    double iso_rms    = 1.0;
    double fixed_iso  = -1.0;   // if >0, use this absolute contour
    double tau_start  = -1.0;   // if >=0, display MD time tau = tau_start + idx*tau_step
    double tau_step   =  0.0;
    // Heatmap mode: fix two dimensions, show 2D colour map
    bool   do_heatmap  = false;
    int    slice_dim2  = -1, sliceval2 = 0;
    double heat_scale  = -1.0;  // if >0, fixed symmetric colour scale; else auto
    Coordinate hotsite({0,0,0,0});
    bool has_hotsite = false;

    std::string arg;
    if (GridCmdOptionExists(argv, argv+argc, "--snapdir"))
        snapdir = GridCmdOptionPayload(argv, argv+argc, "--snapdir");
    if (GridCmdOptionExists(argv, argv+argc, "--pngdir"))
        pngdir = GridCmdOptionPayload(argv, argv+argc, "--pngdir");
    if (GridCmdOptionExists(argv, argv+argc, "--first")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--first");
        GridCmdOptionInt(arg, first);
    }
    if (GridCmdOptionExists(argv, argv+argc, "--last")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--last");
        GridCmdOptionInt(arg, last);
    }
    if (GridCmdOptionExists(argv, argv+argc, "--step")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--step");
        GridCmdOptionInt(arg, step);
    }
    if (GridCmdOptionExists(argv, argv+argc, "--slicedim")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--slicedim");
        GridCmdOptionInt(arg, slice_dim);
    }
    if (GridCmdOptionExists(argv, argv+argc, "--sliceval")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--sliceval");
        GridCmdOptionInt(arg, sliceval);
    }
    if (GridCmdOptionExists(argv, argv+argc, "--isosurface")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--isosurface");
        GridCmdOptionFloat(arg, iso_rms);
    }
    if (GridCmdOptionExists(argv, argv+argc, "--fixediso")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--fixediso");
        GridCmdOptionFloat(arg, fixed_iso);
    }
    if (GridCmdOptionExists(argv, argv+argc, "--taustart")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--taustart");
        GridCmdOptionFloat(arg, tau_start);
    }
    if (GridCmdOptionExists(argv, argv+argc, "--taustep")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--taustep");
        GridCmdOptionFloat(arg, tau_step);
    }
    if (GridCmdOptionExists(argv, argv+argc, "--hotsite")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--hotsite");
        GridCmdOptionIntVector(arg, hotsite);
        has_hotsite = true;
    }
    if (GridCmdOptionExists(argv, argv+argc, "--heatmap"))
        do_heatmap = true;
    if (GridCmdOptionExists(argv, argv+argc, "--slicedim2")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--slicedim2");
        GridCmdOptionInt(arg, slice_dim2);
    }
    if (GridCmdOptionExists(argv, argv+argc, "--sliceval2")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--sliceval2");
        GridCmdOptionInt(arg, sliceval2);
    }
    if (GridCmdOptionExists(argv, argv+argc, "--heatscale")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--heatscale");
        GridCmdOptionFloat(arg, heat_scale);
    }

    bool do_png = !pngdir.empty();
    if (do_png) mkdir(pngdir.c_str(), 0755);

    // ── Grid setup ───────────────────────────────────────────────────────────
    auto latt_size   = GridDefaultLatt();
    auto simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
    auto mpi_layout  = GridDefaultMpi();
    GridCartesian grid(latt_size, simd_layout, mpi_layout);
    LatticeComplexD field(&grid);

    // Force components
    struct ForceSpec { std::string prefix; std::string label; };
    std::vector<ForceSpec> forces = {
        { "F_IwasakiGaugeAction_lat",                                             "Gauge_lat"    },
        { "F_IwasakiGaugeAction_smr",                                             "Gauge_smr"    },
        { "F_JacobianAction_smr",                                                 "Jacobian"     },
        { "F_TwoFlavourEvenOddRatioPseudoFermionActiondet_0.0047_det_0.05_lat",  "Ferm0047_lat" },
        { "F_TwoFlavourEvenOddRatioPseudoFermionActiondet_0.0047_det_0.05_smr",  "Ferm0047_smr" },
        { "F_TwoFlavourEvenOddRatioPseudoFermionActiondet_0.05_det_0.1_lat",     "Ferm005_lat"  },
        { "F_TwoFlavourEvenOddRatioPseudoFermionActiondet_0.1_det_0.25_lat",     "Ferm01_lat"   },
        { "F_TwoFlavourEvenOddRatioPseudoFermionActiondet_0.25_det_0.5_lat",     "Ferm025_lat"  },
        { "F_TwoFlavourEvenOddRatioPseudoFermionActiondet_0.5_det_1_lat",        "Ferm05_lat"   },
    };

    // ── Stdout header ─────────────────────────────────────────────────────────
    std::cerr << "idx";
    for (auto& fs : forces) {
        std::cerr << "\t" << fs.label << "_rms";
        if (has_hotsite) std::cerr << "\t" << fs.label << "_hot";
    }
    std::cerr << "\n";

    // ── Persistent render contexts (one GPU context for all frames) ──────────
    std::unique_ptr<RenderCtx>    ctx;      // isosurface mode
    std::unique_ptr<HeatmapCtx>   hctx;    // heatmap mode

    // ── Main loop ─────────────────────────────────────────────────────────────
    for (int idx = first; idx <= last; idx += step) {
        std::cerr << idx;

        for (auto& fs : forces) {
            std::string fname = snapdir + "/" + fs.prefix + "." + std::to_string(idx);

            if (!tryReadFile(field, fname)) {
                std::cerr << "\t-";
                if (has_hotsite) std::cerr << "\t-";
                continue;
            }

            // RMS (real part)
            RealD sumsq = 0.0;
            for (int i = 0; i < grid.gSites(); i++) {
                Coordinate site;
                Lexicographic::CoorFromIndex(site, i, latt_size);
                RealD v = real(peekSite(field, site));
                sumsq += v * v;
            }
            RealD rms = std::sqrt(sumsq / grid.gSites());
            std::cerr << "\t" << rms;

            if (has_hotsite) {
                RealD hval = real(TensorRemove(peekSite(field, hotsite)));
                std::cerr << "\t" << hval;
            }

            // PNG output (isosurface or heatmap)
            if (do_png) {
                // Build title string
                std::ostringstream title;
                title << fs.label << "  ";
                if (tau_start >= 0.0 && tau_step > 0.0) {
                    double tau = tau_start + idx * tau_step;
                    title << std::fixed << std::setprecision(6) << "tau=" << tau;
                } else {
                    title << "idx=" << idx;
                }
                title << "  rms=" << std::scientific << std::setprecision(3) << rms;

                std::ostringstream outpath;
                outpath << pngdir << "/" << fs.label
                        << "_" << std::setfill('0') << std::setw(6) << idx << ".png";

                if (do_heatmap && slice_dim2 >= 0) {
                    // ── Heatmap mode ────────────────────────────────────────
                    // Display dims = the two that are NOT fixed
                    std::vector<int> disp;
                    for (int d = 0; d < 4; d++)
                        if (d != slice_dim && d != slice_dim2) disp.push_back(d);

                    if (!hctx) {
                        double sc = (heat_scale > 0) ? heat_scale : rms * 20.0;
                        // Hotsite projection onto display plane
                        int hx = -1, hy = -1;
                        if (has_hotsite) {
                            hx = hotsite[disp[0]];
                            hy = hotsite[disp[1]];
                        }
                        hctx = std::make_unique<HeatmapCtx>();
                        hctx->init(latt_size[disp[0]], latt_size[disp[1]], sc, hx, hy);
                    }

                    title << "  scale=+-" << std::fixed << std::setprecision(4) << hctx->scale;
                    renderHeatmap(*hctx, field, latt_size,
                                  slice_dim,  sliceval,
                                  slice_dim2, sliceval2,
                                  title.str(), outpath.str());
                } else {
                    // ── Isosurface mode ─────────────────────────────────────
                    double contour = (fixed_iso > 0) ? fixed_iso : iso_rms * rms;
                    title << "  iso=" << contour;

                    if (!ctx) {
                        std::vector<int> disp;
                        for (int d = 0; d < 4; d++) if (d != slice_dim) disp.push_back(d);
                        ctx = std::make_unique<RenderCtx>();
                        ctx->init(latt_size[disp[0]], latt_size[disp[1]], latt_size[disp[2]]);
                    }
                    renderPNG(*ctx, field, latt_size, slice_dim, sliceval,
                              contour, title.str(), outpath.str());
                }
            }
        }
        std::cerr << "\n";
        std::cerr.flush();
    }

    Grid_finalize();
    return EXIT_SUCCESS;
}
