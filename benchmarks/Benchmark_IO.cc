#include "Benchmark_IO.hpp"

#ifndef BENCH_IO_LMIN
#define BENCH_IO_LMIN 8
#endif

#ifndef BENCH_IO_LMAX
#define BENCH_IO_LMAX 32
#endif

#ifndef BENCH_IO_NPASS
#define BENCH_IO_NPASS 10
#endif

#ifdef HAVE_LIME
using namespace Grid;

std::string filestem(const int l)
{
  return "iobench_l" + std::to_string(l);
}

int vol(const int i)
{
  return BENCH_IO_LMIN + 2*i;
}

int volInd(const int l)
{
  return (l - BENCH_IO_LMIN)/2;
}

template <typename Mat>
void stats(Mat &mean, Mat &stdDev, const std::vector<Mat> &data)
{
  auto            nr = data[0].rows(), nc = data[0].cols();
  Eigen::MatrixXd sqSum(nr, nc);
  double          n = static_cast<double>(data.size());

  assert(n > 1.);
  mean  = Mat::Zero(nr, nc);
  sqSum = Mat::Zero(nr, nc);
  for (auto &d: data)
  {
    mean  += d;
    sqSum += d.cwiseProduct(d);
  }
  stdDev = ((sqSum - mean.cwiseProduct(mean)/n)/(n - 1.)).cwiseSqrt();
  mean  /= n;
}

#define grid_printf(...) \
{\
  char _buf[1024];\
  sprintf(_buf, __VA_ARGS__);\
  MSG << _buf;\
}

enum {sRead = 0, sWrite = 1, gRead = 2, gWrite = 3};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int64_t                      threads = GridThread::GetThreads();
  auto                         mpi     = GridDefaultMpi();
  unsigned int                 nVol    = (BENCH_IO_LMAX - BENCH_IO_LMIN)/2 + 1;
  unsigned int                 nRelVol = (BENCH_IO_LMAX - 24)/2 + 1;
  std::vector<Eigen::MatrixXd> perf(BENCH_IO_NPASS, Eigen::MatrixXd::Zero(nVol, 4));
  std::vector<Eigen::VectorXd> avPerf(BENCH_IO_NPASS, Eigen::VectorXd::Zero(4));
  std::vector<int>             latt;

  MSG << "Grid is setup to use " << threads << " threads" << std::endl;
  MSG << "MPI partition " << mpi << std::endl;
  for (unsigned int i = 0; i < BENCH_IO_NPASS; ++i)
  {
    MSG << BIGSEP << std::endl;
    MSG << "Pass " << i + 1 << "/" << BENCH_IO_NPASS << std::endl;
    MSG << BIGSEP << std::endl;
    MSG << SEP << std::endl;
    MSG << "Benchmark std write" << std::endl;
    MSG << SEP << std::endl;
    for (int l = BENCH_IO_LMIN; l <= BENCH_IO_LMAX; l += 2)
    {
      latt = {l*mpi[0], l*mpi[1], l*mpi[2], l*mpi[3]};

      MSG << "-- Local volume " << l << "^4" << std::endl;
      writeBenchmark<LatticeFermion>(latt, filestem(l), stdWrite<LatticeFermion>);
      perf[i](volInd(l), sWrite) = BinaryIO::lastPerf.mbytesPerSecond;
    }

    MSG << SEP << std::endl;
    MSG << "Benchmark std read" << std::endl;
    MSG << SEP << std::endl;
    for (int l = BENCH_IO_LMIN; l <= BENCH_IO_LMAX; l += 2)
    {
      latt = {l*mpi[0], l*mpi[1], l*mpi[2], l*mpi[3]};

      MSG << "-- Local volume " << l << "^4" << std::endl;
      readBenchmark<LatticeFermion>(latt, filestem(l), stdRead<LatticeFermion>);
      perf[i](volInd(l), sRead) = BinaryIO::lastPerf.mbytesPerSecond;
    }

  #ifdef HAVE_LIME
    MSG << SEP << std::endl;
    MSG << "Benchmark Grid C-Lime write" << std::endl;
    MSG << SEP << std::endl;
    for (int l = BENCH_IO_LMIN; l <= BENCH_IO_LMAX; l += 2)
    {
      latt = {l*mpi[0], l*mpi[1], l*mpi[2], l*mpi[3]};

      MSG << "-- Local volume " << l << "^4" << std::endl;
      writeBenchmark<LatticeFermion>(latt, filestem(l), limeWrite<LatticeFermion>);
      perf[i](volInd(l), gWrite) = BinaryIO::lastPerf.mbytesPerSecond;
    }

    MSG << SEP << std::endl;
    MSG << "Benchmark Grid C-Lime read" << std::endl;
    MSG << SEP << std::endl;
    for (int l = BENCH_IO_LMIN; l <= BENCH_IO_LMAX; l += 2)
    {
      latt = {l*mpi[0], l*mpi[1], l*mpi[2], l*mpi[3]};

      MSG << "-- Local volume " << l << "^4" << std::endl;
      readBenchmark<LatticeFermion>(latt, filestem(l), limeRead<LatticeFermion>);
      perf[i](volInd(l), gRead) = BinaryIO::lastPerf.mbytesPerSecond;
    }
#endif
    avPerf[i].fill(0.);
    for (int f = 0; f < 4; ++f)
    for (int l = 24; l <= BENCH_IO_LMAX; l += 2)
    {
      avPerf[i](f) += perf[i](volInd(l), f);
    }
    avPerf[i] /= nRelVol;
  }

  Eigen::MatrixXd mean(nVol, 4), stdDev(nVol, 4), rob(nVol, 4);
  Eigen::VectorXd avMean(4), avStdDev(4), avRob(4);
  //  double          n = BENCH_IO_NPASS;

  stats(mean, stdDev, perf);
  stats(avMean, avStdDev, avPerf);
  rob.fill(100.);
  rob -= 100.*stdDev.cwiseQuotient(mean.cwiseAbs());
  avRob.fill(100.);
  avRob -= 100.*avStdDev.cwiseQuotient(avMean.cwiseAbs());

  MSG << BIGSEP << std::endl;
  MSG << "SUMMARY" << std::endl;
  MSG << BIGSEP << std::endl;
  MSG << "Summary of individual results (all results in MB/s)." << std::endl;
  MSG << "Every second colum gives the standard deviation of the previous column." << std::endl;
  MSG << std::endl;
  grid_printf("%4s %12s %12s %12s %12s %12s %12s %12s %12s\n",
              "L", "std read", "std dev", "std write", "std dev",
              "Grid read", "std dev", "Grid write", "std dev");
  for (int l = BENCH_IO_LMIN; l <= BENCH_IO_LMAX; l += 2)
  {
    grid_printf("%4d %12.1f %12.1f %12.1f %12.1f %12.1f %12.1f %12.1f %12.1f\n",
                l, mean(volInd(l), sRead), stdDev(volInd(l), sRead),
                mean(volInd(l), sWrite), stdDev(volInd(l), sWrite),
                mean(volInd(l), gRead), stdDev(volInd(l), gRead),
                mean(volInd(l), gWrite), stdDev(volInd(l), gWrite));
  }
  MSG << std::endl;
  MSG << "Robustness of individual results, in %. (rob = 100% - std dev / mean)" << std::endl;
  MSG << std::endl;
  grid_printf("%4s %12s %12s %12s %12s\n",
              "L", "std read", "std write", "Grid read", "Grid write");
  for (int l = BENCH_IO_LMIN; l <= BENCH_IO_LMAX; l += 2)
  {
    grid_printf("%4d %12.1f %12.1f %12.1f %12.1f\n",
                l, rob(volInd(l), sRead), rob(volInd(l), sWrite),
                rob(volInd(l), gRead), rob(volInd(l), gWrite));
  }
  MSG << std::endl;
  MSG << "Summary of results averaged over local volumes 24^4-" << BENCH_IO_LMAX << "^4 (all results in MB/s)." << std::endl;
  MSG << "Every second colum gives the standard deviation of the previous column." << std::endl;
  MSG << std::endl;
  grid_printf("%12s %12s %12s %12s %12s %12s %12s %12s\n",
              "std read", "std dev", "std write", "std dev",
              "Grid read", "std dev", "Grid write", "std dev");
  grid_printf("%12.1f %12.1f %12.1f %12.1f %12.1f %12.1f %12.1f %12.1f\n",
              avMean(sRead), avStdDev(sRead), avMean(sWrite), avStdDev(sWrite),
              avMean(gRead), avStdDev(gRead), avMean(gWrite), avStdDev(gWrite));
  MSG << std::endl;
  MSG << "Robustness of volume-averaged results, in %. (rob = 100% - std dev / mean)" << std::endl;
  MSG << std::endl;
  grid_printf("%12s %12s %12s %12s\n",
              "std read", "std write", "Grid read", "Grid write");
  grid_printf("%12.1f %12.1f %12.1f %12.1f\n",
              avRob(sRead), avRob(sWrite), avRob(gRead), avRob(gWrite));

  Grid_finalize();

  return EXIT_SUCCESS;
}
#else
int main(int argc,char ** argv){}
#endif
