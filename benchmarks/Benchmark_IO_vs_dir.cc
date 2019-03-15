#include "Benchmark_IO.hpp"

#define MSG std::cout << GridLogMessage
#define SEP \
"============================================================================="

using namespace Grid;
using namespace QCD;

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
  MSG << "Grid is setup to use " << threads << " threads" << std::endl;
  MSG << SEP << std::endl;
  MSG << "Benchmark double precision Lime write" << std::endl;
  MSG << SEP << std::endl;
  for (auto &d: dir)
  {
    MSG << "-- Directory " << d << std::endl;
    writeBenchmark<LatticeFermion>(GridDefaultLatt(), d + "/ioBench", limeWrite<LatticeFermion>, Ls, rb);
  }

  MSG << SEP << std::endl;
  MSG << "Benchmark double precision Lime read" << std::endl;
  MSG << SEP << std::endl;
  for (auto &d: dir)
  {
    MSG << "-- Directory " << d << std::endl;
    readBenchmark<LatticeFermion>(GridDefaultLatt(), d + "/ioBench", limeRead<LatticeFermion>, Ls, rb);
  }

  MSG << SEP << std::endl;
  MSG << "Benchmark single precision Lime write" << std::endl;
  MSG << SEP << std::endl;
  for (auto &d: dir)
  {
    MSG << "-- Directory " << d << std::endl;
    writeBenchmark<LatticeFermionF>(GridDefaultLatt(), d + "/ioBench", limeWrite<LatticeFermionF>, Ls, rb);
  }

  MSG << SEP << std::endl;
  MSG << "Benchmark single precision Lime read" << std::endl;
  MSG << SEP << std::endl;
  for (auto &d: dir)
  {
    MSG << "-- Directory " << d << std::endl;
    readBenchmark<LatticeFermionF>(GridDefaultLatt(), d + "/ioBench", limeRead<LatticeFermionF>, Ls, rb);
  }

  Grid_finalize();

  return EXIT_SUCCESS;
}
