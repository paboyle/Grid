#ifndef Benchmark_IO_hpp_
#define Benchmark_IO_hpp_

#include <Grid/Grid.h>

#define MSG std::cout << GridLogMessage
#define SEP \
"============================================================================="

namespace Grid {

template <typename Field>
using WriterFn = std::function<void(const std::string, Field &)> ;
template <typename Field>
using ReaderFn = std::function<void(Field &, const std::string)>;

template <typename Field>
void limeWrite(const std::string filestem, Field &vec)
{
  emptyUserRecord   record;
  QCD::ScidacWriter binWriter(vec._grid->IsBoss());

  binWriter.open(filestem + ".bin");
  binWriter.writeScidacFieldRecord(vec, record);
  binWriter.close();
}

template <typename Field>
void limeRead(Field &vec, const std::string filestem)
{
  emptyUserRecord   record;
  QCD::ScidacReader binReader;

  binReader.open(filestem + ".bin");
  binReader.readScidacFieldRecord(vec, record);
  binReader.close();
}

inline void makeGrid(std::shared_ptr<GridBase> &gPt, 
                     const std::shared_ptr<GridCartesian> &gBasePt,
                     const unsigned int Ls = 1, const bool rb = false)
{
  if (rb)
  {
    if (Ls > 1)
    {
      gPt.reset(QCD::SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, gBasePt.get()));
    }
    else
    {
      gPt.reset(QCD::SpaceTimeGrid::makeFourDimRedBlackGrid(gBasePt.get()));
    }
  }
  else
  {
    if (Ls > 1)
    {
        gPt.reset(QCD::SpaceTimeGrid::makeFiveDimGrid(Ls, gBasePt.get()));
    }
    else
    {
        gPt = gBasePt;
    }
  }
}

template <typename Field>
void writeBenchmark(const std::vector<int> &latt, const std::string filename,
                    const WriterFn<Field> &write, 
                    const unsigned int Ls = 1, const bool rb = false)
{
  auto                           mpi  = GridDefaultMpi();
  auto                           simd = GridDefaultSimd(latt.size(), Field::vector_type::Nsimd());
  std::shared_ptr<GridCartesian> gBasePt(QCD::SpaceTimeGrid::makeFourDimGrid(latt, simd, mpi));
  std::shared_ptr<GridBase>      gPt;

  makeGrid(gPt, gBasePt, Ls, rb);

  GridBase                       *g = gPt.get();
  GridParallelRNG                rng(g);
  Field                          vec(g);

  random(rng, vec);
  write(filename, vec);
}

template <typename Field>
void readBenchmark(const std::vector<int> &latt, const std::string filename,
                   const ReaderFn<Field> &read, 
                   const unsigned int Ls = 1, const bool rb = false)
{
  auto                           mpi  = GridDefaultMpi();
  auto                           simd = GridDefaultSimd(latt.size(), Field::vector_type::Nsimd());
  std::shared_ptr<GridCartesian> gBasePt(QCD::SpaceTimeGrid::makeFourDimGrid(latt, simd, mpi));
  std::shared_ptr<GridBase>      gPt;

  makeGrid(gPt, gBasePt, Ls, rb);

  GridBase                       *g = gPt.get();
  Field                          vec(g);

  read(vec, filename);
}

}

#endif // Benchmark_IO_hpp_
