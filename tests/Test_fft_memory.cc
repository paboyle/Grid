/*************************************************************************************

 Test_fft_memory.cc

 Memory growth test for PlannedFFT on a spin-colour matrix (propagator) field.

 The test creates a single PlannedFFT object (which allocates FFTW plans once),
 then repeatedly applies FFT_all_dim to the same propagator 400 times.

 If PlannedFFT is working correctly the RSS should remain flat after the first
 iteration — no new plans, no new deviceVector allocations beyond the per-call
 pencil buffer which is freed at the end of each FFT_dim_execute call.

 Build exactly like any other Grid test, e.g.:
   make Test_fft_memory
 or compile manually:
   $(CXX) $(CXXFLAGS) Test_fft_memory.cc -o Test_fft_memory $(LDFLAGS)

*************************************************************************************/

#include <Grid/Grid.h>
using namespace Grid;

// --------------------------------------------------------------------------
// Helper: read RSS (resident set size) in kB from /proc/self/status.
// Returns 0 on platforms where /proc is unavailable.
// --------------------------------------------------------------------------
static long getCPURSSKb()
{
  long rss = 0;
  FILE *fp = fopen("/proc/self/status", "r");
  if (!fp) return -1;
  char line[256];
  while (fgets(line, sizeof(line), fp)) {
    if (strncmp(line, "VmRSS:", 6) == 0) {
      sscanf(line + 6, "%ld", &rss);
      break;
    }
  }
  fclose(fp);
  return rss;
}

static long getGPUUsedMb()
{
#if defined(GRID_CUDA)
  size_t free_bytes  = 0;
  size_t total_bytes = 0;
  cudaError_t err = cudaMemGetInfo(&free_bytes, &total_bytes);
  if (err != cudaSuccess) return -1;
  return (long)((total_bytes - free_bytes) / (1024 * 1024));
 
#elif defined(GRID_HIP)
  size_t free_bytes  = 0;
  size_t total_bytes = 0;
  hipError_t err = hipMemGetInfo(&free_bytes, &total_bytes);
  if (err != hipSuccess) return -1;
  return (long)((total_bytes - free_bytes) / (1024 * 1024));
 
#else
  return -1;   // CPU-only build: no GPU to query
#endif
}

// ============================================================
//  Convenience struct — one snapshot of both sides
// ============================================================
struct MemSnapshot {
  long cpu_rss_kb;   // host RSS in kB  (-1 if unavailable)
  long gpu_used_mb;  // device used in MB (-1 if no GPU)
};
 
static MemSnapshot takeSnapshot()
{
  MemSnapshot s;
  s.cpu_rss_kb  = getCPURSSKb();
  s.gpu_used_mb = getGPUUsedMb();
  return s;
}
 
// ============================================================
//  Pretty-print one row of the monitoring table
// ============================================================
static void printRow(int iter,
                     const MemSnapshot &now,
                     const MemSnapshot &prev)
{
  long cpu_delta = (now.cpu_rss_kb  >= 0 && prev.cpu_rss_kb  >= 0)
                 ? now.cpu_rss_kb  - prev.cpu_rss_kb  : 0;
  long gpu_delta = (now.gpu_used_mb >= 0 && prev.gpu_used_mb >= 0)
                 ? now.gpu_used_mb - prev.gpu_used_mb : 0;
 
  // Sign prefix so deltas are unambiguous
  auto sign = [](long v) -> const char* { return v >= 0 ? "+" : ""; };
 
  std::cout << GridLogMessage
            << std::setw(6)  << iter
            << "  CPU: " << std::setw(10) << now.cpu_rss_kb  << " kB"
            << "  (" << sign(cpu_delta) << std::setw(7) << cpu_delta << " kB)"
            << "  GPU: " << std::setw(7)  << now.gpu_used_mb << " MB"
            << "  (" << sign(gpu_delta) << std::setw(5) << gpu_delta << " MB)"
            << "\n";
}

// ============================================================

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage
            << "Grid is setup to use " << threads << " threads" << std::endl;

  // ------------------------------------------------------------------
  // Grid setup — use whatever lattice/mpi/simd was passed on the CLI,
  // e.g.  --grid 8.8.8.8 --mpi 1.1.1.1
  // ------------------------------------------------------------------
  Coordinate latt_size   = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(Nd, vComplexD::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();

  GridCartesian GRID(latt_size, simd_layout, mpi_layout);

  int vol = 1;
  for (int d = 0; d < (int)latt_size.size(); d++) vol *= latt_size[d];

  std::cout << GridLogMessage << "Lattice : ";
  for (int d = 0; d < Nd; d++) std::cout << latt_size[d] << " ";
  std::cout << std::endl;

  // ------------------------------------------------------------------
  // Propagator field: SpinColourMatrix = 12x12 complex, i.e.
  // LatticePropagatorD (= Lattice<iSpinColourMatrix<vComplexD>>).
  // This is the standard QCD quark propagator type.
  // ------------------------------------------------------------------
  LatticePropagatorD prop(&GRID);
  
  // ------------------------------------------------------------------
  // Fill the propagator with a momentum-space plane wave,
  // following the pattern from Test_fft.cc.
  // We set each spin-colour component (a,b) to exp(i * sum_mu p_mu x_mu)
  // with a fixed momentum p = (1,2,1,2).
  // ------------------------------------------------------------------
  Coordinate pvec({1, 2, 1, 2});

  LatticeComplexD phase(&GRID);
  LatticeComplexD coor(&GRID);
  ComplexD ci(0.0, 1.0);

  phase = Zero();
  for (int mu = 0; mu < Nd; mu++) {
    RealD TwoPiL = M_PI * 2.0 / latt_size[mu];
    LatticeCoordinate(coor, mu);
    phase = phase + (TwoPiL * pvec[mu]) * coor;
  }
  phase = exp(phase * ci);   // e^{i p.x}

  // Broadcast the phase into every spin-colour matrix entry
  prop = Zero();
  prop = prop + phase;

  std::cout << GridLogMessage
            << "Propagator norm2 = " << norm2(prop) << std::endl;

  // ------------------------------------------------------------------
  // Baseline snapshot BEFORE PlannedFFT construction
  // ------------------------------------------------------------------
  MemSnapshot snap_before_plan = takeSnapshot();
  std::cout << GridLogMessage
            << "[mem] Before PlannedFFT construction"
            << "  CPU: " << snap_before_plan.cpu_rss_kb  << " kB"
            << "  GPU: " << snap_before_plan.gpu_used_mb << " MB"
            << std::endl;
 
  // ------------------------------------------------------------------
  // Create the PlannedFFT — plans are allocated here ONCE for all
  // dimensions and stored inside the object.
  // ------------------------------------------------------------------
  PlannedFFT<iSpinColourMatrix<vComplexD>> plannedFFT(&GRID);

  // ------------------------------------------------------------------
  // Snapshot AFTER plan construction — this is the true baseline
  // for the loop, because cufftPlanMany itself grabs device memory.
  // ------------------------------------------------------------------
  MemSnapshot snap_after_plan = takeSnapshot();
  std::cout << GridLogMessage
            << "[mem] After  PlannedFFT construction"
            << "  CPU: " << snap_after_plan.cpu_rss_kb  << " kB"
            << "  GPU: " << snap_after_plan.gpu_used_mb << " MB"
            << "  (plan overhead:"
            << "  CPU +" << snap_after_plan.cpu_rss_kb - snap_before_plan.cpu_rss_kb << " kB"
            << "  GPU +" << snap_after_plan.gpu_used_mb - snap_before_plan.gpu_used_mb << " MB)"
            << std::endl;

  MemoryManager::Print();
  // ------------------------------------------------------------------
  // 400-iteration loop.
  // Each iteration computes the full 4d forward FFT of `prop`.
  // We deliberately do NOT cache the result — we always start from
  // the same `prop` so the FFT is recomputed identically each time.
  // The point is to watch memory, not correctness.
  // ------------------------------------------------------------------
  const int Niter = 40;
  const int Niter2 = 32;

  // Print header for the memory table
  std::cout << GridLogMessage
              << "\n"
              << std::setw(6)  << "iter"
              << "  CPU: " << std::setw(10) << "RSS[kB]"
              << "  (  delta  )"
              << "  GPU: " << std::setw(7)  << "used[MB]"
              << "  (delta)"
              << "\n";
 
  MemSnapshot snap_prev = snap_after_plan;

  for (int i = 0; i < Niter; i++) {
    std::vector<LatticePropagatorD> G;

    for (int j = 0; j < Niter2; j++) {
      LatticePropagatorD prop_fft(&GRID);

      // Full 4d forward FFT using the pre-built plans
      plannedFFT.FFT_all_dim(prop_fft, prop, FFT::forward);

      G.push_back(prop_fft);
    }
 
    // cudaMemGetInfo reflects the state *after* any pooled frees have
    // been committed, so this is accurate without an explicit sync —
    // FFT_dim_execute already calls accelerator_barrier() internally.
    MemSnapshot snap_now = takeSnapshot();
    printRow(i, snap_now, snap_prev);
    MemoryManager::Print();
    snap_prev = snap_now;
  }

  // ------------------------------------------------------------------
  // Summary
  // ------------------------------------------------------------------
  MemSnapshot snap_final = takeSnapshot();
 
  long cpu_growth = snap_final.cpu_rss_kb  - snap_after_plan.cpu_rss_kb;
  long gpu_growth = snap_final.gpu_used_mb - snap_after_plan.gpu_used_mb;
 
  std::cout << GridLogMessage
            << "\n==== Memory summary (baseline = after plan construction) ====\n"
            << "  CPU RSS growth over " << Niter << " FFTs : "
            << cpu_growth << " kB"
            << (cpu_growth == 0 ? "  OK" : "  *** GROWING ***") << "\n"
            << "  GPU used growth over " << Niter << " FFTs : "
            << gpu_growth << " MB"
            << (gpu_growth == 0 ? "  OK" : "  *** GROWING ***") << "\n"
            << "  Note: first-call watermark from pool fill is expected and benign.\n"
            << "        A leak shows as continuous growth beyond iter ~2-3.\n";
 
  Grid_finalize();
  return 0;
}
