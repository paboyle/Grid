/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/DiskVector.hpp

Copyright (C) 2015-2018

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
        // write to disk and cache
        T &operator=(const T &obj) const
        {
            DV_DEBUG_MSG(&master_, "writing to " << i_);
            master_.cacheInsert(i_, obj);
            master_.save(master_.filename(i_), obj);
            
            return master_.cachePtr_->at(i_);
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
    std::string                                dirname_;
    unsigned int                               size_, cacheSize_;
    double                                     access_{0.}, hit_{0.};
    bool                                       clean_;
    // using pointers to allow modifications when class is const
    // semantic: const means data unmodified, but cache modification allowed
    std::unique_ptr<std::map<unsigned int, T>> cachePtr_;
    std::unique_ptr<std::deque<unsigned int>>  loadsPtr_;                
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
using EigenDiskVectorMat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

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
        std::ifstream              f(filename, std::ios::binary);
        std::vector<unsigned char> hash(SHA256_DIGEST_LENGTH);
        Eigen::Index               nRow, nCol;
        size_t                     matSize;
        double                     t;

        f.read(reinterpret_cast<char *>(hash.data()), hash.size()*sizeof(unsigned char));
        f.read(reinterpret_cast<char *>(&nRow), sizeof(Eigen::Index));
        f.read(reinterpret_cast<char *>(&nCol), sizeof(Eigen::Index));
        obj.resize(nRow, nCol);
        matSize = nRow*nCol*sizeof(T);
        t  = -usecond();
        f.read(reinterpret_cast<char *>(obj.data()), matSize);
        t += usecond();
        DV_DEBUG_MSG(this, "Eigen read " << matSize/t*1.0e6/1024/1024 << " MB/s");
        auto check = GridChecksum::sha256(obj.data(), matSize);
        DV_DEBUG_MSG(this, "Eigen sha256 " << GridChecksum::sha256_string(check));
        if (hash != check)
        {
            HADRONS_ERROR(Io, "checksum failed")
        }
    }

    virtual void save(const std::string filename, const EigenDiskVectorMat<T> &obj) const
    {
        std::ofstream              f(filename, std::ios::binary);
        std::vector<unsigned char> hash(SHA256_DIGEST_LENGTH);
        Eigen::Index               nRow, nCol;
        size_t                     matSize;
        double                     t;
        
        nRow    = obj.rows();
        nCol    = obj.cols();
        matSize = nRow*nCol*sizeof(T);
        hash    = GridChecksum::sha256(obj.data(), matSize);
        DV_DEBUG_MSG(this, "Eigen sha256 " << GridChecksum::sha256_string(hash));
        f.write(reinterpret_cast<char *>(hash.data()), hash.size()*sizeof(unsigned char));
        f.write(reinterpret_cast<char *>(&nRow), sizeof(Eigen::Index));
        f.write(reinterpret_cast<char *>(&nCol), sizeof(Eigen::Index));
        t  = -usecond();
        f.write(reinterpret_cast<const char *>(obj.data()), matSize);
        t += usecond();
        DV_DEBUG_MSG(this, "Eigen write " << matSize/t*1.0e6/1024/1024 << " MB/s");
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
, cachePtr_(new std::map<unsigned int, T>())
, loadsPtr_(new std::deque<unsigned int>())
{
    struct stat s;

    if(stat(dirname.c_str(), &s) == 0)
    {
        HADRONS_ERROR(Io, "directory '" + dirname + "' already exists")
    }
    mkdir(dirname);
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
    auto &cache  = *cachePtr_;
    auto &loads  = *loadsPtr_;

    DV_DEBUG_MSG(this, "accessing " << i << " (RO)");

    if (i >= size_)
    {
        HADRONS_ERROR(Size, "index out of range");
    }
    const_cast<double &>(access_)++;
    if (cache.find(i) == cache.end())
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

    return cache.at(i);
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
    auto &cache = *cachePtr_;
    auto &loads = *loadsPtr_;

    if (cache.size() >= cacheSize_)
    {
        DV_DEBUG_MSG(this, "evicting " << loads.front());
        cache.erase(loads.front());
        loads.pop_front();
    }
}

template <typename T>
void DiskVectorBase<T>::fetch(const unsigned int i) const
{
    auto &cache = *cachePtr_;
    auto &loads = *loadsPtr_;
    struct stat s;

    DV_DEBUG_MSG(this, "loading " << i << " from disk");

    evict();
    if(stat(filename(i).c_str(), &s) != 0)
    {
        HADRONS_ERROR(Io, "disk vector element " + std::to_string(i) + " uninitialised");
    }
    load(cache[i], filename(i));
    loads.push_back(i);
}

template <typename T>
void DiskVectorBase<T>::cacheInsert(const unsigned int i, const T &obj) const
{
    auto &cache = *cachePtr_;
    auto &loads = *loadsPtr_;

    evict();
    cache[i] = obj;
    loads.push_back(i);

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
