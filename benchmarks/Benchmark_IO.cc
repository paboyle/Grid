#include <Grid/Grid.h>
#ifdef HAVE_LIME

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

#define MSG cout << GridLogMessage
#define SEP \
"============================================================================="
#ifndef BENCH_IO_LMAX
#define BENCH_IO_LMAX 40
#endif

typedef function<void(const string, LatticeFermion &)> WriterFn;
typedef function<void(LatticeFermion &, const string)> ReaderFn;

string filestem(const int l)
{
  return "iobench_l" + to_string(l);
}

void limeWrite(const string filestem, LatticeFermion &vec)
{
  emptyUserRecord record;
  ScidacWriter    binWriter(vec._grid->IsBoss());

  binWriter.open(filestem + ".bin");
  binWriter.writeScidacFieldRecord(vec, record);
  binWriter.close();
}

void limeRead(LatticeFermion &vec, const string filestem)
{
  emptyUserRecord record;
  ScidacReader    binReader;

  binReader.open(filestem + ".bin");
  binReader.readScidacFieldRecord(vec, record);
  binReader.close();
}

void writeBenchmark(const int l, const WriterFn &write)
{
  auto                      mpi  = GridDefaultMpi();
  auto                      simd = GridDefaultSimd(Nd, vComplex::Nsimd());
  vector<int>               latt = {l*mpi[0], l*mpi[1], l*mpi[2], l*mpi[3]};
  unique_ptr<GridCartesian> gPt(SpaceTimeGrid::makeFourDimGrid(latt, simd, mpi));
  GridCartesian             *g = gPt.get();
  GridParallelRNG           rng(g);
  LatticeFermion            vec(g);
  emptyUserRecord           record;
  ScidacWriter              binWriter(g->IsBoss());

  cout << "-- Local volume " << l << "^4" << endl;
  random(rng, vec);
  write(filestem(l), vec);
}

void readBenchmark(const int l, const ReaderFn &read)
{
  auto                      mpi  = GridDefaultMpi();
  auto                      simd = GridDefaultSimd(Nd, vComplex::Nsimd());
  vector<int>               latt = {l*mpi[0], l*mpi[1], l*mpi[2], l*mpi[3]};
  unique_ptr<GridCartesian> gPt(SpaceTimeGrid::makeFourDimGrid(latt, simd, mpi));
  GridCartesian             *g = gPt.get();
  LatticeFermion            vec(g);
  emptyUserRecord           record;
  ScidacReader              binReader;

  cout << "-- Local volume " << l << "^4" << endl;
  read(vec, filestem(l));
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  auto simd = GridDefaultSimd(Nd,vComplex::Nsimd());
  auto mpi  = GridDefaultMpi();

  int64_t threads = GridThread::GetThreads();
  MSG << "Grid is setup to use " << threads << " threads" << endl;
  MSG << SEP << endl;
  MSG << "Benchmark Lime write" << endl;
  MSG << SEP << endl;
  for (int l = 4; l <= BENCH_IO_LMAX; l += 2)
  {
    writeBenchmark(l, limeWrite);
  }

  MSG << "Benchmark Lime read" << endl;
  MSG << SEP << endl;
  for (int l = 4; l <= BENCH_IO_LMAX; l += 2)
  {
    readBenchmark(l, limeRead);
  }

  Grid_finalize();

  return EXIT_SUCCESS;
}
#else
int main (int argc, char ** argv)
{
  return EXIT_SUCCESS;
}
#endif
