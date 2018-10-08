/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/A2AMatrix.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#ifndef A2A_Matrix_hpp_
#define A2A_Matrix_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Grid/Eigen/unsupported/CXX11/Tensor>

#ifndef HADRONS_A2AM_NAME 
#define HADRONS_A2AM_NAME "a2aMatrix"
#endif

#define HADRONS_A2AM_PARALLEL_IO

BEGIN_HADRONS_NAMESPACE

// general A2A matrix set based on Eigen tensors and Grid-allocated memory
// Dimensions:
//   0 - ext - external field (momentum, EM field, ...)
//   1 - str - spin-color structure
//   2 - t   - timeslice
//   3 - i   - left  A2A mode index
//   4 - j   - right A2A mode index
template <typename T>
using A2AMatrixSet = Eigen::TensorMap<Eigen::Tensor<T, 5, Eigen::RowMajor>>;

/******************************************************************************
 *                      Abstract class for A2A kernels                        *
 ******************************************************************************/
template <typename T, typename Field>
class A2AKernel
{
public:
    A2AKernel(void) = default;
    virtual ~A2AKernel(void) = default;
    virtual void operator()(A2AMatrixSet<T> &m, const Field *left, const Field *right,
                          const unsigned int orthogDim, double &time) = 0;
    virtual double flops(const unsigned int blockSizei, const unsigned int blockSizej) = 0;
    virtual double bytes(const unsigned int blockSizei, const unsigned int blockSizej) = 0;
};

/******************************************************************************
 *                  Class to handle A2A matrix block HDF5 I/O                 *
 ******************************************************************************/
template <typename T>
class A2AMatrixIo
{
public:
    // constructors
    A2AMatrixIo(void) = default;
    A2AMatrixIo(std::string filename, std::string dataname, 
                const unsigned int nt, const unsigned int ni,
                const unsigned int nj);
    // destructor
    ~A2AMatrixIo(void) = default;
    // file allocation
    template <typename MetadataType>
    void initFile(const MetadataType &d, const unsigned int chunkSize);
    // block I/O
    void saveBlock(const T *data, const unsigned int i, const unsigned int j,
                   const unsigned int blockSizei, const unsigned int blockSizej);
    void saveBlock(const A2AMatrixSet<T> &m, const unsigned int ext, const unsigned int str,
                   const unsigned int i, const unsigned int j);
private:
    std::string  filename_, dataname_;
    unsigned int nt_, ni_, nj_;
};

/******************************************************************************
 *                  Wrapper for A2A matrix block computation                  *
 ******************************************************************************/
template <typename T, typename Field, typename MetadataType, typename TIo = T>
class A2AMatrixBlockComputation
{
private:
    struct IoHelper
    {
        A2AMatrixIo<TIo> io;
        MetadataType     md;
        unsigned int     e, s, i, j;
    };
    typedef std::function<std::string(const unsigned int, const unsigned int)>  FilenameFn;
    typedef std::function<MetadataType(const unsigned int, const unsigned int)> MetadataFn;
public:
    // constructor
    A2AMatrixBlockComputation(GridBase *grid,
                              const unsigned int orthogDim,
                              const unsigned int next,
                              const unsigned int nstr,
                              const unsigned int blockSize,
                              const unsigned int cacheBlockSize,
                              TimerArray *tArray = nullptr);
    // execution
    void execute(const std::vector<Field> &left, 
                 const std::vector<Field> &right,
                 A2AKernel<T, Field> &kernel,
                 const FilenameFn &ionameFn,
                 const FilenameFn &filenameFn,
                 const MetadataFn &metadataFn);
private:
    // I/O handler
    void saveBlock(const A2AMatrixSet<TIo> &m, IoHelper &h);
private:
    TimerArray            *tArray_;
    GridBase              *grid_;
    unsigned int          orthogDim_, nt_, next_, nstr_, blockSize_, cacheBlockSize_;
    Vector<T>             mCache_;
    Vector<TIo>           mBuf_;
    std::vector<IoHelper> nodeIo_;
};

/******************************************************************************
 *                     A2AMatrixIo template implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename T>
A2AMatrixIo<T>::A2AMatrixIo(std::string filename, std::string dataname, 
                            const unsigned int nt, const unsigned int ni,
                            const unsigned int nj)
: filename_(filename), dataname_(dataname)
, nt_(nt), ni_(ni), nj_(nj)
{}

// file allocation /////////////////////////////////////////////////////////////
template <typename T>
template <typename MetadataType>
void A2AMatrixIo<T>::initFile(const MetadataType &d, const unsigned int chunkSize)
{
#ifdef HAVE_HDF5
    std::vector<hsize_t>    dim = {static_cast<hsize_t>(nt_), 
                                   static_cast<hsize_t>(ni_), 
                                   static_cast<hsize_t>(nj_)},
                            chunk = {static_cast<hsize_t>(nt_), 
                                     static_cast<hsize_t>(chunkSize), 
                                     static_cast<hsize_t>(chunkSize)};
    H5NS::DataSpace         dataspace(dim.size(), dim.data());
    H5NS::DataSet           dataset;
    H5NS::DSetCreatPropList plist;
    
    // create empty file just with metadata
    {
        Hdf5Writer writer(filename_);
        write(writer, dataname_, d);
    }

    // create the dataset
    Hdf5Reader reader(filename_);

    push(reader, dataname_);
    auto &group = reader.getGroup();
    plist.setChunk(chunk.size(), chunk.data());
    dataset = group.createDataSet(HADRONS_A2AM_NAME, Hdf5Type<T>::type(), dataspace, plist);
#else
    HADRONS_ERROR(Implementation, "all-to-all matrix I/O needs HDF5 library");
#endif
}

// block I/O ///////////////////////////////////////////////////////////////////
template <typename T>
void A2AMatrixIo<T>::saveBlock(const T *data, 
                               const unsigned int i, 
                               const unsigned int j,
                               const unsigned int blockSizei,
                               const unsigned int blockSizej)
{
#ifdef HAVE_HDF5
    Hdf5Reader           reader(filename_);
    std::vector<hsize_t> count = {nt_, blockSizei, blockSizej},
                         offset = {0, static_cast<hsize_t>(i),
                                   static_cast<hsize_t>(j)},
                         stride = {1, 1, 1},
                         block  = {1, 1, 1}; 
    H5NS::DataSpace      memspace(count.size(), count.data()), dataspace;
    H5NS::DataSet        dataset;
    size_t               shift;

    push(reader, dataname_);
    auto &group = reader.getGroup();
    dataset     = group.openDataSet(HADRONS_A2AM_NAME);
    dataspace   = dataset.getSpace();
    dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                              stride.data(), block.data());
    dataset.write(data, Hdf5Type<T>::type(), memspace, dataspace);
#else
    HADRONS_ERROR(Implementation, "all-to-all matrix I/O needs HDF5 library");
#endif
}

template <typename T>
void A2AMatrixIo<T>::saveBlock(const A2AMatrixSet<T> &m,
                               const unsigned int ext, const unsigned int str,
                               const unsigned int i, const unsigned int j)
{
    unsigned int blockSizei = m.dimension(3);
    unsigned int blockSizej = m.dimension(4);
    unsigned int nstr       = m.dimension(1);
    size_t       offset     = (ext*nstr + str)*nt_*blockSizei*blockSizej;

    saveBlock(m.data() + offset, i, j, blockSizei, blockSizej);
}

/******************************************************************************
 *               A2AMatrixBlockComputation template implementation            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename T, typename Field, typename MetadataType, typename TIo>
A2AMatrixBlockComputation<T, Field, MetadataType, TIo>
::A2AMatrixBlockComputation(GridBase *grid,
                            const unsigned int orthogDim,
                            const unsigned int next, 
                            const unsigned int nstr,
                            const unsigned int blockSize, 
                            const unsigned int cacheBlockSize,
                            TimerArray *tArray)
: grid_(grid), nt_(grid->GlobalDimensions()[orthogDim]), orthogDim_(orthogDim)
, next_(next), nstr_(nstr), blockSize_(blockSize), cacheBlockSize_(cacheBlockSize)
, tArray_(tArray)
{
    mCache_.resize(nt_*next_*nstr_*cacheBlockSize_*cacheBlockSize_);
    mBuf_.resize(nt_*next_*nstr_*blockSize_*blockSize_);
}

#define START_TIMER(name) if (tArray_) tArray_->startTimer(name)
#define STOP_TIMER(name)  if (tArray_) tArray_->stopTimer(name)
#define GET_TIMER(name)   ((tArray_ != nullptr) ? tArray_->getDTimer(name) : 0.)

// execution ///////////////////////////////////////////////////////////////////
template <typename T, typename Field, typename MetadataType, typename TIo>
void A2AMatrixBlockComputation<T, Field, MetadataType, TIo>
::execute(const std::vector<Field> &left, const std::vector<Field> &right,
          A2AKernel<T, Field> &kernel, const FilenameFn &ionameFn,
          const FilenameFn &filenameFn, const MetadataFn &metadataFn)
{
    //////////////////////////////////////////////////////////////////////////
    // i,j   is first  loop over blockSize_ factors
    // ii,jj is second loop over cacheBlockSize_ factors for high perf contractions
    // iii,jjj are loops within cacheBlock
    // Total index is sum of these  i+ii+iii etc...
    //////////////////////////////////////////////////////////////////////////
    int    N_i = left.size();
    int    N_j = right.size();
    double flops, bytes, t_kernel;
    double nodes = grid_->NodeCount();
    
    int NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0);
    int NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);

    for(int i=0;i<N_i;i+=blockSize_)
    for(int j=0;j<N_j;j+=blockSize_)
    {
        // Get the W and V vectors for this block^2 set of terms
        int N_ii = MIN(N_i-i,blockSize_);
        int N_jj = MIN(N_j-j,blockSize_);
        A2AMatrixSet<TIo> mBlock(mBuf_.data(), next_, nstr_, nt_, N_ii, N_jj);

        LOG(Message) << "All-to-all matrix block " 
                     << j/blockSize_ + NBlock_j*i/blockSize_ + 1 
                     << "/" << NBlock_i*NBlock_j << " [" << i <<" .. " 
                     << i+N_ii-1 << ", " << j <<" .. " << j+N_jj-1 << "]" 
                     << std::endl;
        // Series of cache blocked chunks of the contractions within this block
        flops    = 0.0;
        bytes    = 0.0;
        t_kernel = 0.0;
        for(int ii=0;ii<N_ii;ii+=cacheBlockSize_)
        for(int jj=0;jj<N_jj;jj+=cacheBlockSize_)
        {
            double t;
            int N_iii = MIN(N_ii-ii,cacheBlockSize_);
            int N_jjj = MIN(N_jj-jj,cacheBlockSize_);
            A2AMatrixSet<T> mCacheBlock(mCache_.data(), next_, nstr_, nt_, N_iii, N_jjj);

            START_TIMER("kernel");
            kernel(mCacheBlock, &left[i+ii], &right[j+jj], orthogDim_, t);
            STOP_TIMER("kernel");
            t_kernel += t;
            flops    += kernel.flops(N_iii, N_jjj);
            bytes    += kernel.bytes(N_iii, N_jjj);

            START_TIMER("cache copy");
            parallel_for_nest5(int e =0;e<next_;e++)
            for(int s =0;s< nstr_;s++)
            for(int t =0;t< nt_;t++)
            for(int iii=0;iii< N_iii;iii++)
            for(int jjj=0;jjj< N_jjj;jjj++)
            {
                mBlock(e,s,t,ii+iii,jj+jjj) = mCacheBlock(e,s,t,iii,jjj);
            }
            STOP_TIMER("cache copy");
        }

        // perf
        LOG(Message) << "Kernel perf " << flops/t_kernel/1.0e3/nodes 
                     << " Gflop/s/node " << std::endl;
        LOG(Message) << "Kernel perf " << bytes/t_kernel*1.0e6/1024/1024/1024/nodes 
                     << " GB/s/node "  << std::endl;

        // IO
        double       blockSize, ioTime;
        unsigned int myRank = grid_->ThisRank(), nRank  = grid_->RankCount();
    
        LOG(Message) << "Writing block to disk" << std::endl;
        ioTime = -GET_TIMER("IO: write block");
        START_TIMER("IO: total");
        makeFileDir(filenameFn(0, 0), grid_);
#ifdef HADRONS_A2AM_PARALLEL_IO
        grid_->Barrier();
        // make task list for current node
        nodeIo_.clear();
        for(int f = myRank; f < next_*nstr_; f += nRank)
        {
            IoHelper h;

            h.i  = i;
            h.j  = j;
            h.e  = f/nstr_;
            h.s  = f % nstr_;
            h.io = A2AMatrixIo<TIo>(filenameFn(h.e, h.s), 
                                    ionameFn(h.e, h.s), nt_, N_i, N_j);
            h.md = metadataFn(h.e, h.s);
            nodeIo_.push_back(h);
        }
        // parallel IO
        for (auto &h: nodeIo_)
        {
            saveBlock(mBlock, h);
        }
        grid_->Barrier();
#else
        // serial IO, for testing purposes only
        for(int e = 0; e < next_; e++)
        for(int s = 0; s < nstr_; s++)
        {
            IoHelper h;

            h.i  = i;
            h.j  = j;
            h.e  = e;
            h.s  = s;
            h.io = A2AMatrixIo<TIo>(filenameFn(h.e, h.s), 
                                    ionameFn(h.e, h.s), nt_, N_i, N_j);
            h.md = metadataFn(h.e, h.s);
            saveBlock(mfBlock, h);
        }
#endif
        STOP_TIMER("IO: total");
        blockSize  = static_cast<double>(next_*nstr_*nt_*N_ii*N_jj*sizeof(TIo));
        ioTime    += GET_TIMER("IO: write block");
        LOG(Message) << "HDF5 IO done " << sizeString(blockSize) << " in "
                     << ioTime  << " us (" 
                     << blockSize/ioTime*1.0e6/1024/1024
                     << " MB/s)" << std::endl;
    }
}

// I/O handler /////////////////////////////////////////////////////////////////
template <typename T, typename Field, typename MetadataType, typename TIo>
void A2AMatrixBlockComputation<T, Field, MetadataType, TIo>
::saveBlock(const A2AMatrixSet<TIo> &m, IoHelper &h)
{
    if ((h.i == 0) and (h.j == 0))
    {
        START_TIMER("IO: file creation");
        h.io.initFile(h.md, blockSize_);
        STOP_TIMER("IO: file creation");
    }
    START_TIMER("IO: write block");
    h.io.saveBlock(m, h.e, h.s, h.i, h.j);
    STOP_TIMER("IO: write block");
}

#undef START_TIMER
#undef STOP_TIMER
#undef GET_TIMER

END_HADRONS_NAMESPACE

#endif // A2A_Matrix_hpp_
