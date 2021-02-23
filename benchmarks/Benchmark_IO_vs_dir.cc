#include "Benchmark_IO.hpp"
#ifdef HAVE_LIME
using namespace Grid;

int main (int argc, char ** argv)
{
  std::vector<std::string> dir;
  unsigned int             Ls;
  bool                     rb;
  if (argc < 4)
  {
    std::cerr << "usage: " << argv[0] << " <Ls> <RB {0|1}> <dir1> [<dir2> ... <dirn>] [Grid options]";
    std::cerr << std::endl;
  }
  Ls = std::stoi(argv[1]);
  rb = (std::string(argv[2]) == "1");
  for (unsigned int i = 3; i < argc; ++i)
  {
    std::string a = argv[i];

    if (a[0] != '-')
    {
      dir.push_back(std::string(argv[i]));
    }
    else
    {
      break;
    }
  }
  Grid_init(&argc,&argv);

  int64_t threads = GridThread::GetThreads();
  auto    mpi     = GridDefaultMpi();

  MSG << "Grid is setup to use " << threads << " threads" << std::endl;
  MSG << "MPI partition " << mpi << std::endl;

  MSG << SEP << std::endl;
  MSG << "Benchmark Grid std write" << std::endl;
  MSG << SEP << std::endl;
  for (auto &d: dir)
  {
    MSG << "-- Directory " << d << std::endl;
    writeBenchmark<LatticeFermion>(GridDefaultLatt(), d + "/ioBench", 
                                   stdWrite<LatticeFermion>, Ls, rb);
  }
  MSG << SEP << std::endl;
  MSG << "Benchmark Grid std read" << std::endl;
  MSG << SEP << std::endl;
  for (auto &d: dir)
  {
    MSG << "-- Directory " << d << std::endl;
    readBenchmark<LatticeFermion>(GridDefaultLatt(), d + "/ioBench", 
                                  stdRead<LatticeFermion>, Ls, rb);
  }

#ifdef HAVE_LIME
  MSG << SEP << std::endl;
  MSG << "Benchmark Grid C-Lime write" << std::endl;
  MSG << SEP << std::endl;
  for (auto &d: dir)
  {
    MSG << "-- Directory " << d << std::endl;
    writeBenchmark<LatticeFermion>(GridDefaultLatt(), d + "/ioBench", 
                                   limeWrite<LatticeFermion>, Ls, rb);
  }
  MSG << SEP << std::endl;
  MSG << "Benchmark Grid C-Lime read" << std::endl;
  MSG << SEP << std::endl;
  for (auto &d: dir)
  {
    MSG << "-- Directory " << d << std::endl;
    readBenchmark<LatticeFermion>(GridDefaultLatt(), d + "/ioBench", 
                                  limeRead<LatticeFermion>, Ls, rb);
  }
#endif

  // MSG << SEP << std::endl;
  // MSG << "Benchmark single precision Lime write" << std::endl;
  // MSG << SEP << std::endl;
  // for (auto &d: dir)
  // {
  //   MSG << "-- Directory " << d << std::endl;
  //   writeBenchmark<LatticeFermionF>(GridDefaultLatt(), d + "/ioBench", limeWrite<LatticeFermionF>, Ls, rb);
  // }

  // MSG << SEP << std::endl;
  // MSG << "Benchmark single precision Lime read" << std::endl;
  // MSG << SEP << std::endl;
  // for (auto &d: dir)
  // {
  //   MSG << "-- Directory " << d << std::endl;
  //   readBenchmark<LatticeFermionF>(GridDefaultLatt(), d + "/ioBench", limeRead<LatticeFermionF>, Ls, rb);
  // }

  Grid_finalize();

  return EXIT_SUCCESS;
}
#else
int main(int argc,char ** argv){}
#endif
