/*************************************************************************************

  fft5d.cc — Fourier analysis of a time series of 4-D lattice scalar fields.

  Designed for force-norm files from FTHMC (one RealD per site, SCIDAC or binary).
  Assembles a (4+1)-D structure [trajectory][t][z][y][x] using Grid lattice fields
  and performs three FFT analyses:

    --fft spatial
        4-D spatial FFT using Grid's FFT class (MPI+GPU parallel via Cshift/cufft).
        Outputs shell-averaged P(|k|^2) averaged over trajectories, and a
        per-mode table P(kt,kz,ky,kx).

    --fft traj
        1-D trajectory-axis FFT at each lattice site using FFTW3 locally at each
        MPI rank.  Outputs site-averaged P(f_traj) and autocorrelation C(lag).
        Results are combined with GlobalSumVector.

    --fft all
        Both of the above, plus a 2-D cross spectrum P(f_traj, |k_spatial|^2):
        spatial FFT per trajectory, then trajectory FFTW on the momentum-space series.

  Normalisations:
    Spatial FFT  (volume V): P_spatial(k) = |F(k)|^2 / V^2
      => (1/V) * sum_k P(k) = site mean-square per trajectory
    Trajectory FFT (Ntraj):  P_traj(f) = |F(f)|^2 / Ntraj^2
      => sum_f P(f) = site mean-square over trajectories

  Build: add to examples/Make.inc (see bottom of this file), then make.

  Usage (follows Grid conventions):
    fft5d --grid Lx.Ly.Lz.Lt [--mpi Px.Py.Pz.Pt] \
          [--fft spatial|traj|all] [--format scidac|binary] \
          [--output PREFIX] file1 file2 ...

*************************************************************************************/

#include <Grid/Grid.h>
// Grid's FFT.h uses cufft on CUDA builds; for the trajectory axis we also need
// CPU-side FFTW3 (already linked via -lfftw3 in GRID_LIBS).
#include <fftw3.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace Grid;

// ─────────────────────────────────────────────────────────────────────────────
// Global coordinate of a local site (osite, lane) on this MPI rank
// ─────────────────────────────────────────────────────────────────────────────
static void globalCoor(int osite, int lane, GridCartesian* g, Coordinate& gc)
{
    Coordinate oc, ic;
    Lexicographic::CoorFromIndex(oc, osite, g->_rdimensions);
    Lexicographic::CoorFromIndex(ic, lane,  g->_simd_layout);
    for (int d = 0; d < Nd; d++)
        gc[d] = g->_processor_coor[d] * g->_ldimensions[d]
              + oc[d] * g->_simd_layout[d] + ic[d];
}

// ─────────────────────────────────────────────────────────────────────────────
// Shell-average a power LatticeRealD over |k|^2 bins (MPI-aware via GlobalSumVector)
// ─────────────────────────────────────────────────────────────────────────────
static void shellAverage(const LatticeRealD& power, GridCartesian* g,
                         double norm, std::ostream& ofs)
{
    const int Nsimd = vRealD::Nsimd();
    Coordinate fdims = g->_fdimensions;
    int maxk2 = 0;
    for (int d = 0; d < Nd; d++) { int h = fdims[d]/2; maxk2 += h*h; }

    std::vector<RealD> psum(maxk2+1, 0.0), cnt(maxk2+1, 0.0);
    {
        auto pv = power.View(CpuRead);
        Coordinate gc(Nd);
        for (int os = 0; os < (int)g->oSites(); os++) {
            for (int lane = 0; lane < Nsimd; lane++) {
                globalCoor(os, lane, g, gc);
                int k2 = 0;
                for (int d = 0; d < Nd; d++) {
                    int kd = std::min(gc[d], fdims[d] - gc[d]);
                    k2 += kd*kd;
                }
                psum[k2] += (RealD)extractLane(lane, pv[os]);
                cnt [k2] += 1.0;
            }
        }
    }
    g->GlobalSumVector(psum.data(), (int)psum.size());
    g->GlobalSumVector(cnt .data(), (int)cnt .size());

    for (int k2 = 0; k2 <= maxk2; k2++)
        if (cnt[k2] > 0)
            ofs << k2 << "  " << std::sqrt((double)k2) << "  "
                << (psum[k2] / cnt[k2]) * norm << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// Convert LatticeRealD → LatticeComplexD (zero imaginary part)
// ─────────────────────────────────────────────────────────────────────────────
static LatticeComplexD toComplex(const LatticeRealD& r)
{
    std::vector<RealD>    lr;
    unvectorizeToLexOrdArray(lr, r);
    std::vector<ComplexD> lc(lr.size());
    for (size_t i = 0; i < lr.size(); i++) lc[i] = ComplexD(lr[i], 0.0);
    LatticeComplexD c(r.Grid());
    vectorizeFromLexOrdArray(lc, c);
    return c;
}

// ─────────────────────────────────────────────────────────────────────────────
// File readers
// ─────────────────────────────────────────────────────────────────────────────
static LatticeComplexD readScidac(const std::string& fname, GridCartesian* g)
{
    LatticeRealD field(g);
    emptyUserRecord rec;
    ScidacReader RD;
    RD.open(fname);
    RD.readScidacFieldRecord(field, rec);
    RD.close();
    return toComplex(field);
}

static LatticeComplexD readBinary(const std::string& fname, GridCartesian* g)
{
    // Raw IEEE doubles in lex order [x][y][z][t] written serially.
    // Rank 0 reads the file, broadcasts to all ranks via GlobalSumVector.
    Coordinate fdims = g->_fdimensions;
    long vol4 = 1; for (int d = 0; d < Nd; d++) vol4 *= fdims[d];

    std::vector<RealD> buf(vol4, 0.0);
    if (g->IsBoss()) {
        std::ifstream f(fname, std::ios::binary);
        if (!f) throw std::runtime_error("Cannot open: " + fname);
        f.read(reinterpret_cast<char*>(buf.data()), vol4 * sizeof(RealD));
    }
    g->GlobalSumVector(buf.data(), (int)vol4);   // broadcast from rank 0

    LatticeRealD field(g);
    Coordinate ldims = g->_ldimensions, pcoor = g->_processor_coor;
    for (long ls = 0; ls < g->lSites(); ls++) {
        Coordinate lc;
        g->LocalIndexToLocalCoor(ls, lc);
        long glex = 0, stride = 1;
        for (int d = 0; d < Nd; d++) {
            glex += (pcoor[d]*ldims[d] + lc[d]) * stride;
            stride *= fdims[d];
        }
        pokeLocalSite(buf[glex], field, lc);
    }
    return toComplex(field);
}

// ─────────────────────────────────────────────────────────────────────────────
// Load all files into a trajectory vector
// ─────────────────────────────────────────────────────────────────────────────
static std::vector<LatticeComplexD>
loadFiles(const std::vector<std::string>& files, GridCartesian* g,
          const std::string& fmt)
{
    std::vector<LatticeComplexD> traj;
    traj.reserve(files.size());
    for (int n = 0; n < (int)files.size(); n++) {
        std::cout << GridLogMessage << "[" << n << "] reading " << files[n] << "\n";
        if (fmt == "scidac") traj.push_back(readScidac(files[n], g));
        else                 traj.push_back(readBinary (files[n], g));
    }
    return traj;
}

// ─────────────────────────────────────────────────────────────────────────────
// Spatial FFT analysis
// ─────────────────────────────────────────────────────────────────────────────
static void analyzeSpatial(const std::vector<LatticeComplexD>& traj,
                            GridCartesian* g, const std::string& pfx)
{
    int  Ntraj = (int)traj.size();
    long vol4  = 1; for (int d = 0; d < Nd; d++) vol4 *= g->_fdimensions[d];

    FFT theFFT(g);
    LatticeRealD pavg(g); pavg = Zero();
    for (int n = 0; n < Ntraj; n++) {
        LatticeComplexD fk(g);
        theFFT.FFT_all_dim(fk, traj[n], FFT::forward);
        // Evaluate product before applying real() — real() is not defined for
        // unevaluated LatticeBinaryExpression.
        LatticeComplexD fk_sq(g); fk_sq = conjugate(fk) * fk;
        std::vector<ComplexD> lc; unvectorizeToLexOrdArray(lc, fk_sq);
        std::vector<RealD> lr(lc.size());
        for (size_t i = 0; i < lc.size(); i++) lr[i] = lc[i].real();
        LatticeRealD pk(g); vectorizeFromLexOrdArray(lr, pk);
        pavg += pk;
    }
    pavg *= (1.0 / Ntraj);

    // Shell-averaged spectrum
    {
        std::ofstream fs(pfx + "_spatial_shell.dat");
        fs << "# k2  |k|  P_shell_avg   (P = |F|^2 / V^2 / shell_count)\n";
        shellAverage(pavg, g, 1.0 / ((double)vol4 * vol4), fs);
        std::cout << GridLogMessage << "Written: " << pfx << "_spatial_shell.dat\n";
    }

    // Per-mode table (rank 0 only, via peekSite)
    if (g->IsBoss()) {
        std::ofstream fm(pfx + "_spatial_modes.dat");
        fm << "# kt  kz  ky  kx  k2  |k|  P_traj_avg\n";
        Coordinate fdims = g->_fdimensions, site(Nd);
        double norm = 1.0 / ((double)vol4 * vol4);
        for (site[3]=0; site[3]<fdims[3]; site[3]++) { int kt=std::min(site[3],fdims[3]-site[3]);
        for (site[2]=0; site[2]<fdims[2]; site[2]++) { int kz=std::min(site[2],fdims[2]-site[2]);
        for (site[1]=0; site[1]<fdims[1]; site[1]++) { int ky=std::min(site[1],fdims[1]-site[1]);
        for (site[0]=0; site[0]<fdims[0]; site[0]++) { int kx=std::min(site[0],fdims[0]-site[0]);
            RealD p; peekSite(p, pavg, site);
            int k2 = kt*kt+kz*kz+ky*ky+kx*kx;
            fm << kt<<" "<<kz<<" "<<ky<<" "<<kx<<" "<<k2<<" "
               <<std::sqrt((double)k2)<<" "<<p*norm<<"\n";
        }}}}
        std::cout << GridLogMessage << "Written: " << pfx << "_spatial_modes.dat\n";
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Local FFTW trajectory FFT: each rank FFTs its own local sites independently
// ─────────────────────────────────────────────────────────────────────────────
static void trajFFT(const std::vector<LatticeComplexD>& in,
                    std::vector<LatticeComplexD>& out)
{
    int    Ntraj  = (int)in.size();
    GridBase* g   = in[0].Grid();
    long   lsites = g->lSites();

    // Pack: buf[traj * lsites + lsite]  (traj varies slowly, site varies fast)
    std::vector<fftw_complex> ibuf((long)Ntraj * lsites);
    for (int n = 0; n < Ntraj; n++) {
        std::vector<ComplexD> lc;
        unvectorizeToLexOrdArray(lc, in[n]);
        for (long s = 0; s < lsites; s++) {
            ibuf[(long)n*lsites + s][0] = lc[s].real();
            ibuf[(long)n*lsites + s][1] = lc[s].imag();
        }
    }

    // lsites transforms of length Ntraj, stride=lsites, dist=1
    std::vector<fftw_complex> obuf((long)Ntraj * lsites);
    int n1[1] = {Ntraj};
    fftw_plan p = fftw_plan_many_dft(
        1, n1, (int)lsites,
        ibuf.data(), nullptr, (int)lsites, 1,
        obuf.data(), nullptr, (int)lsites, 1,
        FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    // vector::assign triggers _M_fill_assign which needs a default constructor;
    // Lattice has none.  Use explicit push_back instead.
    out.clear(); out.reserve(Ntraj);
    for (int k = 0; k < Ntraj; k++) out.emplace_back(g);
    for (int k = 0; k < Ntraj; k++) {
        std::vector<ComplexD> lc(lsites);
        for (long s = 0; s < lsites; s++)
            lc[s] = ComplexD(obuf[(long)k*lsites+s][0], obuf[(long)k*lsites+s][1]);
        vectorizeFromLexOrdArray(lc, out[k]);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Trajectory-axis FFT analysis
// ─────────────────────────────────────────────────────────────────────────────
static void analyzeTraj(const std::vector<LatticeComplexD>& traj,
                         GridCartesian* g, const std::string& pfx)
{
    int  Ntraj = (int)traj.size();
    int  Nf    = Ntraj / 2 + 1;
    long vol4  = 1; for (int d = 0; d < Nd; d++) vol4 *= g->_fdimensions[d];

    std::vector<LatticeComplexD> ftraj;
    trajFFT(traj, ftraj);

    // P(k) = (1/vol4) * sum_sites |F_traj(k)|^2 / Ntraj^2
    std::vector<RealD> Pavg(Ntraj, 0.0);
    for (int k = 0; k < Ntraj; k++) {
        LatticeComplexD tmp(g); tmp = conjugate(ftraj[k]) * ftraj[k];
        std::vector<ComplexD> lc; unvectorizeToLexOrdArray(lc, tmp);
        std::vector<RealD> lr(lc.size());
        for (size_t i = 0; i < lc.size(); i++) lr[i] = lc[i].real();
        LatticeRealD pk(g); vectorizeFromLexOrdArray(lr, pk);
        RealD s = 0.0;
        for (long ls = 0; ls < g->lSites(); ls++) {
            Coordinate lc; g->LocalIndexToLocalCoor(ls, lc);
            RealD v; peekLocalSite(v, pk, lc);
            s += v;
        }
        g->GlobalSum(s);
        Pavg[k] = s / ((double)vol4 * Ntraj * Ntraj);
    }

    {
        std::ofstream f(pfx + "_traj_power.dat");
        f << "# k  freq  P_avg   (sum_k P = site mean-square)\n";
        for (int k = 0; k < Nf; k++)
            f << k << "  " << (double)k/Ntraj << "  " << Pavg[k] << "\n";
        std::cout << GridLogMessage << "Written: " << pfx << "_traj_power.dat\n";
    }

    // Autocorrelation via IFFT of the power spectrum
    std::vector<fftw_complex> Pc(Nf);
    for (int k = 0; k < Nf; k++) { Pc[k][0] = Pavg[k]; Pc[k][1] = 0.0; }
    std::vector<double> acorr(Ntraj, 0.0);
    fftw_plan ip = fftw_plan_dft_c2r(1, &Ntraj, Pc.data(), acorr.data(), FFTW_ESTIMATE);
    fftw_execute(ip); fftw_destroy_plan(ip);
    double c0 = acorr[0] / Ntraj;
    {
        std::ofstream f(pfx + "_traj_autocorr.dat");
        f << "# lag  C(lag)  C(lag)/C(0)\n";
        for (int lag = 0; lag < Ntraj/2; lag++) {
            double c = acorr[lag] / Ntraj;
            f << lag << "  " << c << "  " << (c0 > 0 ? c/c0 : 0.0) << "\n";
        }
        std::cout << GridLogMessage << "Written: " << pfx << "_traj_autocorr.dat\n";
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Full 5-D analysis: spatial FFT → trajectory FFT → 2-D cross spectrum
// ─────────────────────────────────────────────────────────────────────────────
static void analyzeAll5D(const std::vector<LatticeComplexD>& traj,
                          GridCartesian* g, const std::string& pfx)
{
    int  Ntraj = (int)traj.size();
    int  Nf    = Ntraj / 2 + 1;
    long vol4  = 1; for (int d = 0; d < Nd; d++) vol4 *= g->_fdimensions[d];
    Coordinate fdims = g->_fdimensions;

    // Step 1: spatial FFT for each trajectory
    FFT theFFT(g);
    std::vector<LatticeComplexD> sfft(Ntraj, LatticeComplexD(g));
    for (int n = 0; n < Ntraj; n++)
        theFFT.FFT_all_dim(sfft[n], traj[n], FFT::forward);

    // Step 2: trajectory FFTW on the momentum-space series
    std::vector<LatticeComplexD> tfft;
    trajFFT(sfft, tfft);

    // Step 3: P(f_traj, |k_spatial|^2) shell-averaged
    int maxk2 = 0;
    for (int d = 0; d < Nd; d++) { int h = fdims[d]/2; maxk2 += h*h; }

    long nb = (long)Nf * (maxk2+1);
    std::vector<RealD> p2d(nb, 0.0), cnt2d(nb, 0.0);
    const int Nsimd = vComplexD::Nsimd();
    double norm = 1.0 / ((double)Ntraj * Ntraj * vol4 * vol4);

    for (int kf = 0; kf < Nf; kf++) {
        LatticeComplexD tmp(g); tmp = conjugate(tfft[kf]) * tfft[kf];
        std::vector<ComplexD> lc2; unvectorizeToLexOrdArray(lc2, tmp);
        std::vector<RealD> lr2(lc2.size());
        for (size_t i = 0; i < lc2.size(); i++) lr2[i] = lc2[i].real();
        LatticeRealD pk(g); vectorizeFromLexOrdArray(lr2, pk);
        auto pv = pk.View(CpuRead);
        Coordinate gc(Nd);
        for (int os = 0; os < (int)g->oSites(); os++) {
            for (int lane = 0; lane < Nsimd; lane++) {
                globalCoor(os, lane, g, gc);
                int k2 = 0;
                for (int d = 0; d < Nd; d++) {
                    int kd = std::min(gc[d], fdims[d] - gc[d]);
                    k2 += kd*kd;
                }
                long b = (long)kf*(maxk2+1) + k2;
                p2d [b] += (RealD)extractLane(lane, pv[os]);
                cnt2d[b] += 1.0;
            }
        }
    }
    g->GlobalSumVector(p2d .data(), (int)nb);
    g->GlobalSumVector(cnt2d.data(), (int)nb);

    if (g->IsBoss()) {
        std::ofstream f(pfx + "_5d_cross.dat");
        f << "# k_traj  freq_traj  k2_spatial  |k_spatial|  P_avg\n";
        for (int kf = 0; kf < Nf; kf++) {
            double freq = (double)kf / Ntraj;
            for (int k2 = 0; k2 <= maxk2; k2++) {
                long b = (long)kf*(maxk2+1) + k2;
                if (cnt2d[b] > 0)
                    f << kf << "  " << freq << "  " << k2 << "  "
                      << std::sqrt((double)k2) << "  "
                      << (p2d[b]/cnt2d[b]) * norm << "\n";
            }
            f << "\n";   // blank line between kf slices for gnuplot pm3d
        }
        std::cout << GridLogMessage << "Written: " << pfx << "_5d_cross.dat\n";
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Main
// ─────────────────────────────────────────────────────────────────────────────
static void usage(const char* prog)
{
    std::cerr
        << "Usage: " << prog << " --grid Lx.Ly.Lz.Lt [--mpi Px.Py.Pz.Pt]\n"
        << "       [--fft spatial|traj|all] [--format scidac|binary]\n"
        << "       [--output PREFIX] file1 [file2 ...]\n";
    exit(1);
}

int main(int argc, char** argv)
{
    Grid_init(&argc, &argv);   // consumes --grid, --mpi, etc.

    std::string fftMode = "all";
    std::string fmt     = "scidac";
    std::string pfx     = "fft5d";
    std::vector<std::string> files;

    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if      (a == "--fft"    && i+1 < argc) fftMode = argv[++i];
        else if (a == "--format" && i+1 < argc) fmt     = argv[++i];
        else if (a == "--output" && i+1 < argc) pfx     = argv[++i];
        else if (!a.empty() && a[0] != '-')     files.push_back(a);
        // unknown --flags skipped (may be Grid flags already consumed)
    }

    if (files.empty()) usage(argv[0]);

    GridCartesian* grid = SpaceTimeGrid::makeFourDimGrid(
        GridDefaultLatt(),
        GridDefaultSimd(Nd, vComplexD::Nsimd()),
        GridDefaultMpi());

    Coordinate latt = GridDefaultLatt();
    std::cout << GridLogMessage << "Lattice : "
              << latt[0]<<"x"<<latt[1]<<"x"<<latt[2]<<"x"<<latt[3] << "\n"
              << GridLogMessage << "Ntraj   : " << files.size() << "\n"
              << GridLogMessage << "FFT     : " << fftMode << "\n"
              << GridLogMessage << "Format  : " << fmt << "\n"
              << GridLogMessage << "Output  : " << pfx << "\n";

    auto traj = loadFiles(files, grid, fmt);

    if (fftMode == "spatial" || fftMode == "all") analyzeSpatial(traj, grid, pfx);
    if (fftMode == "traj"    || fftMode == "all") analyzeTraj   (traj, grid, pfx);
    if (fftMode == "all")                         analyzeAll5D  (traj, grid, pfx);

    Grid_finalize();
    return 0;
}

/*
  Add to examples/Make.inc to build:

    bin_PROGRAMS += ... fft5d
    fft5d_SOURCES = fft5d.cc
    fft5d_LDADD   = $(top_builddir)/Grid/libGrid.a
*/
