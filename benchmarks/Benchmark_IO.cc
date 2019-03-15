
#include "Benchmark_IO.hpp"

#ifndef BENCH_IO_LMAX
#define BENCH_IO_LMAX 40
#endif

using namespace Grid;
using namespace QCD;

std::string filestem(const int l)
{
  return "iobench_l" + std::to_string(l);
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int64_t threads = GridThread::GetThreads();
  MSG << "Grid is setup to use " << threads << " threads" << std::endl;
  MSG << SEP << std::endl;
  MSG << "Benchmark Lime write" << std::endl;
  MSG << SEP << std::endl;
  for (int l = 4; l <= BENCH_IO_LMAX; l += 2)
  {
    auto             mpi  = GridDefaultMpi();
    std::vector<int> latt = {l*mpi[0], l*mpi[1], l*mpi[2], l*mpi[3]};

    std::cout << "-- Local volume " << l << "^4" << std::endl;
    writeBenchmark<LatticeFermion>(latt, filestem(l), limeWrite<LatticeFermion>);
  }

  MSG << "Benchmark Lime read" << std::endl;
  MSG << SEP << std::endl;
  for (int l = 4; l <= BENCH_IO_LMAX; l += 2)
  {
    auto             mpi  = GridDefaultMpi();
    std::vector<int> latt = {l*mpi[0], l*mpi[1], l*mpi[2], l*mpi[3]};

    std::cout << "-- Local volume " << l << "^4" << std::endl;
    readBenchmark<LatticeFermion>(latt, filestem(l), limeRead<LatticeFermion>);
  }

  Grid_finalize();

  return EXIT_SUCCESS;
}
