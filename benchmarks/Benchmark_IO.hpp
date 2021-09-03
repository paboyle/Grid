#ifndef Benchmark_IO_hpp_
#define Benchmark_IO_hpp_

#include <Grid/Grid.h>
#define MSG std::cout << GridLogMessage
#define SEP \
"-----------------------------------------------------------------------------"
#define BIGSEP \
"============================================================================="
#ifdef HAVE_LIME

namespace Grid {

template <typename Field>
using WriterFn = std::function<void(const std::string, Field &)> ;
template <typename Field>
using ReaderFn = std::function<void(Field &, const std::string)>;

// AP 06/10/2020: Standard C version in case one is suspicious of the C++ API
// 
// template <typename Field>
// void stdWrite(const std::string filestem, Field &vec)
// {
//   std::string   rankStr = std::to_string(vec.Grid()->ThisRank());
//   std::FILE     *file = std::fopen((filestem + "." + rankStr + ".bin").c_str(), "wb");
//   size_t        size;
//   uint32_t      crc;
//   GridStopWatch ioWatch, crcWatch;

//   size = vec.Grid()->lSites()*sizeof(typename Field::scalar_object);
//   autoView(vec_v, vec, CpuRead);
//   crcWatch.Start();
//   crc = GridChecksum::crc32(vec_v.cpu_ptr, size);
//   std::fwrite(&crc, sizeof(uint32_t), 1, file);
//   crcWatch.Stop();
//   MSG << "Std I/O write: Data CRC32 " << std::hex << crc << std::dec << std::endl;
//   ioWatch.Start();
//   std::fwrite(vec_v.cpu_ptr, sizeof(typename Field::scalar_object), vec.Grid()->lSites(), file);
//   ioWatch.Stop();
//   std::fclose(file);
//   size *= vec.Grid()->ProcessorCount();
//   auto &p = BinaryIO::lastPerf;
//   p.size            = size;
//   p.time            = ioWatch.useconds();
//   p.mbytesPerSecond = size/1024./1024./(ioWatch.useconds()/1.e6);
//   MSG << "Std I/O write: Wrote " << p.size << " bytes in " << ioWatch.Elapsed() 
//       << ", " << p.mbytesPerSecond << " MB/s" << std::endl;
//   MSG << "Std I/O write: checksum overhead " << crcWatch.Elapsed() << std::endl;
// }
//
// template <typename Field>
// void stdRead(Field &vec, const std::string filestem)
// {
//   std::string   rankStr = std::to_string(vec.Grid()->ThisRank());
//   std::FILE     *file = std::fopen((filestem + "." + rankStr + ".bin").c_str(), "rb");
//   size_t        size;
//   uint32_t      crcRead, crcData;
//   GridStopWatch ioWatch, crcWatch;

//   size = vec.Grid()->lSites()*sizeof(typename Field::scalar_object);
//   crcWatch.Start();
//   std::fread(&crcRead, sizeof(uint32_t), 1, file);
//   crcWatch.Stop();
//   {
//     autoView(vec_v, vec, CpuWrite);
//     ioWatch.Start();
//     std::fread(vec_v.cpu_ptr, sizeof(typename Field::scalar_object), vec.Grid()->lSites(), file);
//     ioWatch.Stop();
//     std::fclose(file);
//   }
//   {
//     autoView(vec_v, vec, CpuRead);
//     crcWatch.Start();
//     crcData = GridChecksum::crc32(vec_v.cpu_ptr, size);
//     crcWatch.Stop();
//   }
//   MSG << "Std I/O read: Data CRC32 " << std::hex << crcData << std::dec << std::endl;
//   assert(crcData == crcRead);
//   size *= vec.Grid()->ProcessorCount();
//   auto &p = BinaryIO::lastPerf;
//   p.size            = size;
//   p.time            = ioWatch.useconds();
//   p.mbytesPerSecond = size/1024./1024./(ioWatch.useconds()/1.e6);
//   MSG << "Std I/O read: Read " <<  p.size << " bytes in " << ioWatch.Elapsed() 
//       << ", " << p.mbytesPerSecond << " MB/s" << std::endl;
//   MSG << "Std I/O read: checksum overhead " << crcWatch.Elapsed() << std::endl;
// }

template <typename Field>
void stdWrite(const std::string filestem, Field &vec)
{
  std::string   rankStr = std::to_string(vec.Grid()->ThisRank());
  std::ofstream file(filestem + "." + rankStr + ".bin", std::ios::out | std::ios::binary);
  size_t        size, sizec;
  uint32_t      crc;
  GridStopWatch ioWatch, crcWatch;

  size  = vec.Grid()->lSites()*sizeof(typename Field::scalar_object);
  sizec = size/sizeof(char); // just in case of...
  autoView(vec_v, vec, CpuRead);
  crcWatch.Start();
  crc = GridChecksum::crc32(vec_v.cpu_ptr, size);
  file.write(reinterpret_cast<char *>(&crc), sizeof(uint32_t)/sizeof(char));
  crcWatch.Stop();
  MSG << "Std I/O write: Data CRC32 " << std::hex << crc << std::dec << std::endl;
  ioWatch.Start();
  file.write(reinterpret_cast<char *>(vec_v.cpu_ptr), sizec);
  file.flush();
  ioWatch.Stop();
  size *= vec.Grid()->ProcessorCount();
  auto &p = BinaryIO::lastPerf;
  p.size            = size;
  p.time            = ioWatch.useconds();
  p.mbytesPerSecond = size/1024./1024./(ioWatch.useconds()/1.e6);
  MSG << "Std I/O write: Wrote " << p.size << " bytes in " << ioWatch.Elapsed() 
      << ", " << p.mbytesPerSecond << " MB/s" << std::endl;
  MSG << "Std I/O write: checksum overhead " << crcWatch.Elapsed() << std::endl;
}

template <typename Field>
void stdRead(Field &vec, const std::string filestem)
{
  std::string   rankStr = std::to_string(vec.Grid()->ThisRank());
  std::ifstream file(filestem + "." + rankStr + ".bin", std::ios::in | std::ios::binary);
  size_t        size, sizec;
  uint32_t      crcRead, crcData;
  GridStopWatch ioWatch, crcWatch;

  size  = vec.Grid()->lSites()*sizeof(typename Field::scalar_object);
  sizec = size/sizeof(char); // just in case of...
  crcWatch.Start();
  file.read(reinterpret_cast<char *>(&crcRead), sizeof(uint32_t)/sizeof(char));
  crcWatch.Stop();
  {
    autoView(vec_v, vec, CpuWrite);
    ioWatch.Start();
    file.read(reinterpret_cast<char *>(vec_v.cpu_ptr), sizec);
    ioWatch.Stop();
  }
  {
    autoView(vec_v, vec, CpuRead);
    crcWatch.Start();
    crcData = GridChecksum::crc32(vec_v.cpu_ptr, size);
    crcWatch.Stop();
  }
  MSG << "Std I/O read: Data CRC32 " << std::hex << crcData << std::dec << std::endl;
  assert(crcData == crcRead);
  size *= vec.Grid()->ProcessorCount();
  auto &p = BinaryIO::lastPerf;
  p.size            = size;
  p.time            = ioWatch.useconds();
  p.mbytesPerSecond = size/1024./1024./(ioWatch.useconds()/1.e6);
  MSG << "Std I/O read: Read " <<  p.size << " bytes in " << ioWatch.Elapsed() 
      << ", " << p.mbytesPerSecond << " MB/s" << std::endl;
  MSG << "Std I/O read: checksum overhead " << crcWatch.Elapsed() << std::endl;
}

template <typename Field>
void limeWrite(const std::string filestem, Field &vec)
{
  emptyUserRecord   record;
  ScidacWriter binWriter(vec.Grid()->IsBoss());

  binWriter.open(filestem + ".lime.bin");
  binWriter.writeScidacFieldRecord(vec, record);
  binWriter.close();
}

template <typename Field>
void limeRead(Field &vec, const std::string filestem)
{
  emptyUserRecord   record;
  ScidacReader binReader;

  binReader.open(filestem + ".lime.bin");
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
      gPt.reset(SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, gBasePt.get()));
    }
    else
    {
      gPt.reset(SpaceTimeGrid::makeFourDimRedBlackGrid(gBasePt.get()));
    }
  }
  else
  {
    if (Ls > 1)
    {
        gPt.reset(SpaceTimeGrid::makeFiveDimGrid(Ls, gBasePt.get()));
    }
    else
    {
        gPt = gBasePt;
    }
  }
}

template <typename Field>
void writeBenchmark(const Coordinate &latt, const std::string filename,
                    const WriterFn<Field> &write, 
                    const unsigned int Ls = 1, const bool rb = false)
{
  auto                           mpi  = GridDefaultMpi();
  auto                           simd = GridDefaultSimd(latt.size(), Field::vector_type::Nsimd());
  std::shared_ptr<GridCartesian> gBasePt(SpaceTimeGrid::makeFourDimGrid(latt, simd, mpi));
  std::shared_ptr<GridBase>      gPt;
  std::random_device             rd;

  makeGrid(gPt, gBasePt, Ls, rb);

  GridBase         *g = gPt.get();
  GridParallelRNG  rng(g);
  Field            vec(g);

  rng.SeedFixedIntegers({static_cast<int>(rd()), static_cast<int>(rd()),
                         static_cast<int>(rd()), static_cast<int>(rd()),
                         static_cast<int>(rd()), static_cast<int>(rd()),
                         static_cast<int>(rd()), static_cast<int>(rd())});

  random(rng, vec);
  write(filename, vec);
}

template <typename Field>
void readBenchmark(const Coordinate &latt, const std::string filename,
                   const ReaderFn<Field> &read, 
                   const unsigned int Ls = 1, const bool rb = false)
{
  auto                           mpi  = GridDefaultMpi();
  auto                           simd = GridDefaultSimd(latt.size(), Field::vector_type::Nsimd());
  std::shared_ptr<GridCartesian> gBasePt(SpaceTimeGrid::makeFourDimGrid(latt, simd, mpi));
  std::shared_ptr<GridBase>      gPt;

  makeGrid(gPt, gBasePt, Ls, rb);

  GridBase *g = gPt.get();
  Field    vec(g);

  read(vec, filename);
}

}

#endif //LIME
#endif // Benchmark_IO_hpp_
