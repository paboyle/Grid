/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/DiskVector.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>

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
#ifndef Hadrons_DiskVector_hpp_
#define Hadrons_DiskVector_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <deque>
#include <sys/stat.h>
#include <ftw.h>
#include <unistd.h>

#ifdef DV_DEBUG
#define DV_DEBUG_MSG(dv, stream) LOG(Debug) << "diskvector " << (dv) << ": " << stream << std::endl
#else
#define DV_DEBUG_MSG(dv, stream)
#endif

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                           Abstract base class                              *
 ******************************************************************************/
template <typename T>
class DiskVectorBase
{
public:
    typedef T ObjectType;

    // helper for read/write vector access
    class RwAccessHelper
    {
    public:
        RwAccessHelper(DiskVectorBase<T> &master, const unsigned int i)
        : master_(master), cmaster_(master), i_(i) {}

        // operator=: somebody is trying to store a vector element
        // write to cache and tag as modified
        T &operator=(const T &obj) const
        {
            auto &cache    = *master_.cachePtr_;
            auto &modified = *master_.modifiedPtr_;
            auto &index    = *master_.indexPtr_;

            DV_DEBUG_MSG(&master_, "writing to " << i_);
            master_.cacheInsert(i_, obj);
            modified[index.at(i_)] = true;
            
            return cache[index.at(i_)];
        }

        // implicit cast to const object reference and redirection
        // to the const operator[] for read-only operations
        operator const T&() const
        {
            return cmaster_[i_];
        }
    private:
        DiskVectorBase<T>       &master_;
        const DiskVectorBase<T> &cmaster_;
        const unsigned int      i_;
    };
public:
    DiskVectorBase(const std::string dirname, const unsigned int size = 0,
                   const unsigned int cacheSize = 1, const bool clean = true);
    DiskVectorBase(DiskVectorBase<T> &&v) = default;
    virtual ~DiskVectorBase(void);
    const T & operator[](const unsigned int i) const;
    RwAccessHelper operator[](const unsigned int i);
    double hitRatio(void) const;
    void resetStat(void);
private:
    virtual void load(T &obj, const std::string filename) const = 0;
    virtual void save(const std::string filename, const T &obj) const = 0;
    virtual std::string filename(const unsigned int i) const;
    void evict(void) const;
    void fetch(const unsigned int i) const;
    void cacheInsert(const unsigned int i, const T &obj) const;
    void clean(void);
private:
    std::string                                           dirname_;
    unsigned int                                          size_, cacheSize_;
    double                                                access_{0.}, hit_{0.};
    bool                                                  clean_;
    // using pointers to allow modifications when class is const
    // semantic: const means data unmodified, but cache modification allowed
    std::unique_ptr<std::vector<T>>                       cachePtr_;
    std::unique_ptr<std::vector<bool>>                    modifiedPtr_;
    std::unique_ptr<std::map<unsigned int, unsigned int>> indexPtr_;
    std::unique_ptr<std::stack<unsigned int>>             freePtr_;
    std::unique_ptr<std::deque<unsigned int>>             loadsPtr_;                
};

/******************************************************************************
 *                   Specialisation for serialisable classes                  *
 ******************************************************************************/
template <typename T, typename Reader, typename Writer>
class SerializableDiskVector: public DiskVectorBase<T>
{
public:
    using DiskVectorBase<T>::DiskVectorBase;
private:
    virtual void load(T &obj, const std::string filename) const
    {
        Reader reader(filename);

        read(reader, basename(filename), obj);
    }

    virtual void save(const std::string filename, const T &obj) const
    {
        Writer writer(filename);

        write(writer, basename(filename), obj);
    }
};

/******************************************************************************
 *                      Specialisation for Eigen matrices                     *
 ******************************************************************************/
template <typename T>
using EigenDiskVectorMat = A2AMatrix<T>;

template <typename T>
class EigenDiskVector: public DiskVectorBase<EigenDiskVectorMat<T>>
{
public:
    using DiskVectorBase<EigenDiskVectorMat<T>>::DiskVectorBase;
    typedef EigenDiskVectorMat<T> Matrix;
public:
    T operator()(const unsigned int i, const Eigen::Index j,
                 const Eigen::Index k) const
    {
        return (*this)[i](j, k);
    }
private:
    virtual void load(EigenDiskVectorMat<T> &obj, const std::string filename) const
    {
        std::ifstream f(filename, std::ios::binary);
        uint32_t      crc, check;
        Eigen::Index  nRow, nCol;
        size_t        matSize;
        double        tRead, tHash;

        f.read(reinterpret_cast<char *>(&crc), sizeof(crc));
        f.read(reinterpret_cast<char *>(&nRow), sizeof(nRow));
        f.read(reinterpret_cast<char *>(&nCol), sizeof(nCol));
        obj.resize(nRow, nCol);
        matSize = nRow*nCol*sizeof(T);
        tRead  = -usecond();
        f.read(reinterpret_cast<char *>(obj.data()), matSize);
        tRead += usecond();
        tHash  = -usecond();
#ifdef USE_IPP
        check  = GridChecksum::crc32c(obj.data(), matSize);
#else
        check  = GridChecksum::crc32(obj.data(), matSize);
#endif
        tHash += usecond();
        DV_DEBUG_MSG(this, "Eigen read " << tRead/1.0e6 << " sec " << matSize/tRead*1.0e6/1024/1024 << " MB/s");
        DV_DEBUG_MSG(this, "Eigen crc32 " << std::hex << check << std::dec 
                     << " " << tHash/1.0e6 << " sec " << matSize/tHash*1.0e6/1024/1024 << " MB/s");
        if (crc != check)
        {
            HADRONS_ERROR(Io, "checksum failed")
        }
    }

    virtual void save(const std::string filename, const EigenDiskVectorMat<T> &obj) const
    {
        std::ofstream f(filename, std::ios::binary);
        uint32_t      crc;
        Eigen::Index  nRow, nCol;
        size_t        matSize;
        double        tWrite, tHash;
        
        nRow    = obj.rows();
        nCol    = obj.cols();
        matSize = nRow*nCol*sizeof(T);
        tHash   = -usecond();
#ifdef USE_IPP
        crc     = GridChecksum::crc32c(obj.data(), matSize);
#else
        crc     = GridChecksum::crc32(obj.data(), matSize);
#endif
        tHash  += usecond();
        f.write(reinterpret_cast<char *>(&crc), sizeof(crc));
        f.write(reinterpret_cast<char *>(&nRow), sizeof(nRow));
        f.write(reinterpret_cast<char *>(&nCol), sizeof(nCol));
        tWrite = -usecond();
        f.write(reinterpret_cast<const char *>(obj.data()), matSize);
        tWrite += usecond();
        DV_DEBUG_MSG(this, "Eigen write " << tWrite/1.0e6 << " sec " << matSize/tWrite*1.0e6/1024/1024 << " MB/s");
        DV_DEBUG_MSG(this, "Eigen crc32 " << std::hex << crc << std::dec
                     << " " << tHash/1.0e6 << " sec " << matSize/tHash*1.0e6/1024/1024 << " MB/s");
    }
};

/******************************************************************************
 *                       DiskVectorBase implementation                         *
 ******************************************************************************/
template <typename T>
DiskVectorBase<T>::DiskVectorBase(const std::string dirname, 
                                  const unsigned int size,
                                  const unsigned int cacheSize,
                                  const bool clean)
: dirname_(dirname), size_(size), cacheSize_(cacheSize), clean_(clean)
, cachePtr_(new std::vector<T>(size))
, modifiedPtr_(new std::vector<bool>(size, false))
, indexPtr_(new std::map<unsigned int, unsigned int>())
, freePtr_(new std::stack<unsigned int>)
, loadsPtr_(new std::deque<unsigned int>())
{
    struct stat s;

    if(stat(dirname.c_str(), &s) == 0)
    {
        HADRONS_ERROR(Io, "directory '" + dirname + "' already exists")
    }
    mkdir(dirname);
    for (unsigned int i = 0; i < cacheSize_; ++i)
    {
        freePtr_->push(i);
    }
}

template <typename T>
DiskVectorBase<T>::~DiskVectorBase(void)
{
    if (clean_)
    {
        clean();
    }
}

template <typename T>
const T & DiskVectorBase<T>::operator[](const unsigned int i) const
{
    auto &cache   = *cachePtr_;
    auto &index   = *indexPtr_;
    auto &freeInd = *freePtr_;
    auto &loads   = *loadsPtr_;

    DV_DEBUG_MSG(this, "accessing " << i << " (RO)");

    if (i >= size_)
    {
        HADRONS_ERROR(Size, "index out of range");
    }
    const_cast<double &>(access_)++;
    if (index.find(i) == index.end())
    {
        // cache miss
        DV_DEBUG_MSG(this, "cache miss");
        fetch(i);
    }
    else
    {
        DV_DEBUG_MSG(this, "cache hit");

        auto pos = std::find(loads.begin(), loads.end(), i);

        const_cast<double &>(hit_)++;
        loads.erase(pos);
        loads.push_back(i);
    }

#ifdef DV_DEBUG
    std::string msg;

    for (auto &p: loads)
    {
        msg += std::to_string(p) + " ";
    }
    DV_DEBUG_MSG(this, "in cache: " << msg);
#endif

    return cache[index.at(i)];
}

template <typename T>
typename DiskVectorBase<T>::RwAccessHelper DiskVectorBase<T>::operator[](const unsigned int i)
{
    DV_DEBUG_MSG(this, "accessing " << i << " (RW)");

    if (i >= size_)
    {
        HADRONS_ERROR(Size, "index out of range");
    }

    return RwAccessHelper(*this, i);
}

template <typename T>
double DiskVectorBase<T>::hitRatio(void) const
{
    return hit_/access_;
}

template <typename T>
void DiskVectorBase<T>::resetStat(void)
{
    access_ = 0.;
    hit_    = 0.;
}

template <typename T>
std::string DiskVectorBase<T>::filename(const unsigned int i) const
{
    return dirname_ + "/elem_" + std::to_string(i);
}

template <typename T>
void DiskVectorBase<T>::evict(void) const
{
    auto &cache    = *cachePtr_;
    auto &modified = *modifiedPtr_;
    auto &index    = *indexPtr_;
    auto &freeInd  = *freePtr_;
    auto &loads    = *loadsPtr_;

    if (index.size() >= cacheSize_)
    {
        unsigned int i = loads.front();
        
        DV_DEBUG_MSG(this, "evicting " << i);
        if (modified[index.at(i)])
        {
            DV_DEBUG_MSG(this, "element " << i << " modified, saving to disk");
            save(filename(i), cache[index.at(i)]);
        }
        freeInd.push(index.at(i));
        index.erase(i);
        loads.pop_front();
    }
}

template <typename T>
void DiskVectorBase<T>::fetch(const unsigned int i) const
{
    auto &cache    = *cachePtr_;
    auto &modified = *modifiedPtr_;
    auto &index    = *indexPtr_;
    auto &freeInd  = *freePtr_;
    auto &loads    = *loadsPtr_;

    struct stat s;

    DV_DEBUG_MSG(this, "loading " << i << " from disk");

    evict();
    
    if(stat(filename(i).c_str(), &s) != 0)
    {
        HADRONS_ERROR(Io, "disk vector element " + std::to_string(i) + " uninitialised");
    }
    index[i] = freeInd.top();
    freeInd.pop();
    load(cache[index.at(i)], filename(i));
    loads.push_back(i);
    modified[index.at(i)] = false;
}

template <typename T>
void DiskVectorBase<T>::cacheInsert(const unsigned int i, const T &obj) const
{
    auto &cache    = *cachePtr_;
    auto &modified = *modifiedPtr_;
    auto &index    = *indexPtr_;
    auto &freeInd  = *freePtr_;
    auto &loads    = *loadsPtr_;

    evict();
    index[i] = freeInd.top();
    freeInd.pop();
    cache[index.at(i)] = obj;
    loads.push_back(i);
    modified[index.at(i)] = false;

#ifdef DV_DEBUG
    std::string msg;

    for (auto &p: loads)
    {
        msg += std::to_string(p) + " ";
    }
    DV_DEBUG_MSG(this, "in cache: " << msg);
#endif
}

#ifdef DV_DEBUG
#undef DV_DEBUG_MSG
#endif

template <typename T>
void DiskVectorBase<T>::clean(void)
{
    auto unlink = [](const char *fpath, const struct stat *sb, 
                     int typeflag, struct FTW *ftwbuf)
    {
        int rv = remove(fpath);

        if (rv)
        {
            HADRONS_ERROR(Io, "cannot remove '" + std::string(fpath) + "': "
                          + std::string(std::strerror(errno)));
        }

        return rv;
    };

    nftw(dirname_.c_str(), unlink, 64, FTW_DEPTH | FTW_PHYS);
}

END_HADRONS_NAMESPACE

#endif // Hadrons_DiskVector_hpp_
