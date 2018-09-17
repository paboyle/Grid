/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Environment.hpp

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

#ifndef Hadrons_Environment_hpp_
#define Hadrons_Environment_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Global environment                                 *
 ******************************************************************************/
class Object
{
public:
    Object(void) = default;
    virtual ~Object(void) = default;
};

template <typename T>
class Holder: public Object
{
public:
    Holder(void) = default;
    Holder(T *pt);
    virtual ~Holder(void) = default;
    T &       get(void) const;
    T *       getPt(void) const;
    void      reset(T *pt);
private:
    std::unique_ptr<T> objPt_{nullptr};
};

#define DEFINE_ENV_ALIAS \
inline Environment & env(void) const\
{\
    return Environment::getInstance();\
}

#define DEFINE_ENV_LAMBDA \
auto env = [](void)->Environment &{return Environment::getInstance();}

class Environment
{
    SINGLETON(Environment);
public:
    typedef SITE_SIZE_TYPE                         Size;
    typedef std::unique_ptr<GridCartesian>         GridPt;
    typedef std::unique_ptr<GridRedBlackCartesian> GridRbPt;
    typedef std::unique_ptr<GridParallelRNG>       RngPt;
    enum class Storage {object, cache, temporary};
private:
    struct ObjInfo
    {
        Size                    size{0};
        Storage                 storage{Storage::object};
        unsigned int            Ls{0};
        const std::type_info    *type{nullptr}, *derivedType{nullptr};
        std::string             name;
        int                     module{-1};
        std::unique_ptr<Object> data{nullptr};
    };
    typedef std::pair<size_t, unsigned int>     FineGridKey;
    typedef std::pair<size_t, std::vector<int>> CoarseGridKey;
public:
    // grids
    template <typename VType = vComplex>
    void                    createGrid(const unsigned int Ls);
    template <typename VType = vComplex>
    void                    createCoarseGrid(const std::vector<int> &blockSize,
                                             const unsigned int Ls = 1);
    template <typename VType = vComplex>
    GridCartesian *         getGrid(void) const;
    template <typename VType = vComplex>
    GridRedBlackCartesian * getRbGrid(void) const;
    template <typename VType = vComplex>
    GridCartesian *         getCoarseGrid(const std::vector<int> &blockSize) const;
    template <typename VType = vComplex>
    GridCartesian *         getGrid(const unsigned int Ls) const;
    template <typename VType = vComplex>
    GridRedBlackCartesian * getRbGrid(const unsigned int Ls) const;
    template <typename VType = vComplex>
    GridCartesian *         getCoarseGrid(const std::vector<int> &blockSize,
                                          const unsigned int Ls) const;
    std::vector<int>        getDim(void) const;
    int                     getDim(const unsigned int mu) const;
    unsigned int            getNd(void) const;
    double                  getVolume(void) const;
    // random number generator
    GridParallelRNG *       get4dRng(void) const;
    // general memory management
    void                    addObject(const std::string name,
                                      const int moduleAddress = -1);
    template <typename B, typename T, typename ... Ts>
    void                    createDerivedObject(const std::string name,
                                                const Environment::Storage storage,
                                                const unsigned int Ls,
                                                Ts && ... args);
    template <typename T, typename ... Ts>
    void                    createObject(const std::string name,
                                         const Environment::Storage storage,
                                         const unsigned int Ls,
                                         Ts && ... args);
    void                    setObjectModule(const unsigned int objAddress,
                                            const int modAddress);
    template <typename B, typename T>
    T *                     getDerivedObject(const unsigned int address) const;
    template <typename B, typename T>
    T *                     getDerivedObject(const std::string name) const;
    template <typename T>
    T *                     getObject(const unsigned int address) const;
    template <typename T>
    T *                     getObject(const std::string name) const;
    unsigned int            getMaxAddress(void) const;
    unsigned int            getObjectAddress(const std::string name) const;
    std::string             getObjectName(const unsigned int address) const;
    std::string             getObjectType(const unsigned int address) const;
    std::string             getObjectType(const std::string name) const;
    Size                    getObjectSize(const unsigned int address) const;
    Size                    getObjectSize(const std::string name) const;
    Storage                 getObjectStorage(const unsigned int address) const;
    Storage                 getObjectStorage(const std::string name) const;
    int                     getObjectModule(const unsigned int address) const;
    int                     getObjectModule(const std::string name) const;
    unsigned int            getObjectLs(const unsigned int address) const;
    unsigned int            getObjectLs(const std::string name) const;
    bool                    hasObject(const unsigned int address) const;
    bool                    hasObject(const std::string name) const;
    bool                    hasCreatedObject(const unsigned int address) const;
    bool                    hasCreatedObject(const std::string name) const;
    bool                    isObject5d(const unsigned int address) const;
    bool                    isObject5d(const std::string name) const;
    template <typename T>
    bool                    isObjectOfType(const unsigned int address) const;
    template <typename T>
    bool                    isObjectOfType(const std::string name) const;
    Environment::Size       getTotalSize(void) const;
    void                    freeObject(const unsigned int address);
    void                    freeObject(const std::string name);
    void                    freeAll(void);
    void                    protectObjects(const bool protect);
    bool                    objectsProtected(void) const;
    // print environment content
    void                    printContent(void) const;
private:
    // general
    double                              vol_;
    bool                                protect_{true};
    // grids
    std::vector<int>                    dim_;
    FineGridKey                         defaultGrid_;
    std::map<FineGridKey, GridPt>       grid4d_;
    std::map<FineGridKey, GridPt>       grid5d_;
    std::map<FineGridKey, GridRbPt>     gridRb4d_;
    std::map<FineGridKey, GridRbPt>     gridRb5d_;
    std::map<CoarseGridKey, GridPt>     gridCoarse4d_;
    std::map<CoarseGridKey, GridPt>     gridCoarse5d_;
    unsigned int                        nd_;
    // random number generator
    RngPt                               rng4d_;
    // object store
    std::vector<ObjInfo>                object_;
    std::map<std::string, unsigned int> objectAddress_;
};

/******************************************************************************
 *                       Holder template implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename T>
Holder<T>::Holder(T *pt)
: objPt_(pt)
{}

// access //////////////////////////////////////////////////////////////////////
template <typename T>
T & Holder<T>::get(void) const
{
    return *objPt_.get();
}

template <typename T>
T * Holder<T>::getPt(void) const
{
    return objPt_.get();
}

template <typename T>
void Holder<T>::reset(T *pt)
{
    objPt_.reset(pt);
}

/******************************************************************************
 *                     Environment template implementation                    *
 ******************************************************************************/
// grids ///////////////////////////////////////////////////////////////////////
template <typename VType>
void Environment::createGrid(const unsigned int Ls)
{
    size_t hash = typeHash<VType>();

    if (grid4d_.find({hash, 1}) == grid4d_.end())
    {
        grid4d_[{hash, 1}].reset(
            SpaceTimeGrid::makeFourDimGrid(getDim(), 
                                        GridDefaultSimd(getNd(), VType::Nsimd()),
                                        GridDefaultMpi()));
        gridRb4d_[{hash, 1}].reset(
            SpaceTimeGrid::makeFourDimRedBlackGrid(grid4d_[{hash, 1}].get()));
    }
    if (grid5d_.find({hash, Ls}) == grid5d_.end())
    {
        auto g = grid4d_[{hash, 1}].get();
        
        grid5d_[{hash, Ls}].reset(SpaceTimeGrid::makeFiveDimGrid(Ls, g));
        gridRb5d_[{hash, Ls}].reset(SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, g));
    }
}

template <typename VType>
void Environment::createCoarseGrid(const std::vector<int> &blockSize,
                                   const unsigned int Ls)
{
    int              nd      = getNd();
    std::vector<int> fineDim = getDim(), coarseDim(nd);
    unsigned int     cLs;
    auto             key4d = blockSize, key5d = blockSize;
    size_t           hash  = typeHash<VType>();

    createGrid(Ls);
    for (int d = 0; d < coarseDim.size(); d++)
    {
        coarseDim[d] = fineDim[d]/blockSize[d];
        if (coarseDim[d]*blockSize[d] != fineDim[d])
        {
            HADRONS_ERROR(Size, "Fine dimension " + std::to_string(d) 
                         + " (" + std::to_string(fineDim[d]) 
                         + ") not divisible by coarse dimension ("
                         + std::to_string(coarseDim[d]) + ")"); 
        }
    }
    if (blockSize.size() > nd)
    {
        cLs = Ls/blockSize[nd];
        if (cLs*blockSize[nd] != Ls)
        {
            HADRONS_ERROR(Size, "Fine Ls (" + std::to_string(Ls) 
                         + ") not divisible by coarse Ls ("
                         + std::to_string(cLs) + ")");
        }
    }
    else
    {
        cLs = Ls;
    }
    key4d.resize(nd);
    key5d.push_back(Ls);

    CoarseGridKey hkey4d = {hash, key4d}, hkey5d = {hash, key5d};

    if (gridCoarse4d_.find(hkey4d) == gridCoarse4d_.end())
    {
        gridCoarse4d_[hkey4d].reset(
            SpaceTimeGrid::makeFourDimGrid(coarseDim, 
                GridDefaultSimd(nd, VType::Nsimd()), GridDefaultMpi()));
    }
    if (gridCoarse5d_.find(hkey5d) == gridCoarse5d_.end())
    {
        gridCoarse5d_[hkey5d].reset(
            SpaceTimeGrid::makeFiveDimGrid(cLs, gridCoarse4d_[hkey4d].get()));
    }
}

template <typename VType>
GridCartesian * Environment::getGrid(void) const
{
    auto it = grid4d_.find({typeHash<VType>(), 1});

    if (it != grid4d_.end())
    {
        return it->second.get();
    }
    else
    {
        HADRONS_ERROR(Definition, "no 4D grid for SIMD type '" 
                                  + typeName<VType>() + "'");
    }
}

template <typename VType>
GridRedBlackCartesian * Environment::getRbGrid(void) const
{
    auto it = gridRb4d_.find({typeHash<VType>(), 1});

    if (it != gridRb4d_.end())
    {
        return it->second.get();
    }
    else
    {
        HADRONS_ERROR(Definition, "no 4D red-black grid for SIMD type '" 
                                  + typeName<VType>() + "'");
    }
}

template <typename VType>
GridCartesian * Environment::getCoarseGrid(const std::vector<int> &blockSize) const
{
    std::vector<int> s = blockSize;

    s.resize(getNd());
    auto it = gridCoarse4d_.find({typeHash<VType>(), s});
    if (it != gridCoarse4d_.end())
    {
        return it->second.get();
    }
    else
    {
        HADRONS_ERROR(Definition, "no 4D coarse grid for SIMD type '" 
                                  + typeName<VType>() + "' and block size "
                                  + vecToStr(blockSize));
    }
}

template <typename VType>
GridCartesian * Environment::getGrid(const unsigned int Ls) const
{
    auto it = grid5d_.find({typeHash<VType>(), Ls});

    if (it != grid5d_.end())
    {
        return it->second.get();
    }
    else
    {
        HADRONS_ERROR(Definition, "no 5D grid for SIMD type '" 
                                  + typeName<VType>() + "' and Ls = "
                                  + std::to_string(Ls));
    }
}

template <typename VType>
GridRedBlackCartesian * Environment::getRbGrid(const unsigned int Ls) const
{
    auto it = gridRb5d_.find({typeHash<VType>(), Ls});

    if (it != gridRb5d_.end())
    {
        return it->second.get();
    }
    else
    {
        HADRONS_ERROR(Definition, "no 5D red-black grid for SIMD type '" 
                                  + typeName<VType>() + "' and Ls = "
                                  + std::to_string(Ls));
    }
}

template <typename VType>
GridCartesian * Environment::getCoarseGrid(const std::vector<int> &blockSize,
                                           const unsigned int Ls) const
{
    std::vector<int> s = blockSize;

    s.push_back(Ls);
    auto it = gridCoarse5d_.find({typeHash<VType>(), s});
    if (it != gridCoarse5d_.end())
    {
        return it->second.get();
    }
    else
    {
        HADRONS_ERROR(Definition, "no 5D coarse grid for SIMD type '" 
                                  + typeName<VType>() + "', block size "
                                  + vecToStr(blockSize)
                                  + " and Ls = " + std::to_string(Ls));
    }
}


// general memory management ///////////////////////////////////////////////////
template <typename B, typename T, typename ... Ts>
void Environment::createDerivedObject(const std::string name,
                                      const Environment::Storage storage,
                                      const unsigned int Ls,
                                      Ts && ... args)
{
    if (!hasObject(name))
    {
        addObject(name);
    }
    
    unsigned int address = getObjectAddress(name);
    
    if (!object_[address].data or !objectsProtected())
    {
        MemoryStats memStats;
    
        if (!MemoryProfiler::stats)
        {
            MemoryProfiler::stats = &memStats;
        }
        size_t initMem               = MemoryProfiler::stats->currentlyAllocated;
        object_[address].storage     = storage;
        object_[address].Ls          = Ls;
        object_[address].data.reset(new Holder<B>(new T(std::forward<Ts>(args)...)));
        object_[address].size        = MemoryProfiler::stats->maxAllocated - initMem;
        object_[address].type        = typeIdPt<B>();
        object_[address].derivedType = typeIdPt<T>();
        if (MemoryProfiler::stats == &memStats)
        {
            MemoryProfiler::stats = nullptr;
        }
    }
    // object already exists, no error if it is a cache, error otherwise
    else if ((object_[address].storage               != Storage::cache) or 
             (object_[address].storage               != storage)        or
             (object_[address].name                  != name)           or
             (typeHash(object_[address].type)        != typeHash<B>())  or
             (typeHash(object_[address].derivedType) != typeHash<T>()))
    {
        HADRONS_ERROR_REF(ObjectDefinition, "object '" + name + "' already allocated", address);
    }
}

template <typename T, typename ... Ts>
void Environment::createObject(const std::string name, 
                               const Environment::Storage storage,
                               const unsigned int Ls,
                               Ts && ... args)
{
    createDerivedObject<T, T>(name, storage, Ls, std::forward<Ts>(args)...);
}

template <typename B, typename T>
T * Environment::getDerivedObject(const unsigned int address) const
{
    if (hasObject(address))
    {
        if (hasCreatedObject(address))
        {
            if (auto h = dynamic_cast<Holder<B> *>(object_[address].data.get()))
            {
                if (&typeid(T) == &typeid(B))
                {
                    return dynamic_cast<T *>(h->getPt());
                }
                else
                {
                    if (auto hder = dynamic_cast<T *>(h->getPt()))
                    {
                        return hder;
                    }
                    else
                    {
                        HADRONS_ERROR_REF(ObjectType, "object with address " +
                            std::to_string(address) +
                            " cannot be casted to '" + typeName(&typeid(T)) +
                            "' (has type '" + typeName(&typeid(h->get())) + "')", address);
                    }
                }
            }
            else
            {
                HADRONS_ERROR_REF(ObjectType, "object with address " + 
                            std::to_string(address) +
                            " does not have type '" + typeName(&typeid(B)) +
                            "' (has type '" + getObjectType(address) + "')", address);
            }
        }
        else
        {
            HADRONS_ERROR_REF(ObjectDefinition, "object with address " + 
                              std::to_string(address) + " is empty", address);
        }
    }
    else
    {
        HADRONS_ERROR_REF(ObjectDefinition, "no object with address " + 
                          std::to_string(address), address);
    }
}

template <typename B, typename T>
T * Environment::getDerivedObject(const std::string name) const
{
    return getDerivedObject<B, T>(getObjectAddress(name));
}

template <typename T>
T * Environment::getObject(const unsigned int address) const
{
    return getDerivedObject<T, T>(address);
}

template <typename T>
T * Environment::getObject(const std::string name) const
{
    return getObject<T>(getObjectAddress(name));
}

template <typename T>
bool Environment::isObjectOfType(const unsigned int address) const
{
    if (hasObject(address))
    {
        if (auto h = dynamic_cast<Holder<T> *>(object_[address].data.get()))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        HADRONS_ERROR_REF(ObjectDefinition, "no object with address " 
                          + std::to_string(address), address);
    }
}

template <typename T>
bool Environment::isObjectOfType(const std::string name) const
{
    return isObjectOfType<T>(getObjectAddress(name));
}

END_HADRONS_NAMESPACE

#endif // Hadrons_Environment_hpp_
