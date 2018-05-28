/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Environment.hpp

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

#include <Grid/Hadrons/Global.hpp>

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
public:
    // grids
    void                    createGrid(const unsigned int Ls);
    void                    createCoarseGrid(const std::vector<int> &blockSize,
                                             const unsigned int Ls = 1);
    GridCartesian *         getGrid(const unsigned int Ls = 1) const;
    GridRedBlackCartesian * getRbGrid(const unsigned int Ls = 1) const;
    GridCartesian *         getCoarseGrid(const std::vector<int> &blockSize,
                                          const unsigned int Ls = 1) const;
    std::vector<int>        getDim(void) const;
    int                     getDim(const unsigned int mu) const;
    unsigned int            getNd(void) const;
    double                  getVolume(void) const;
    // random number generator
    void                    setSeed(const std::vector<int> &seed);
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
    double                                 vol_;
    bool                                   protect_{true};
    // grids
    std::vector<int>                       dim_;
    GridPt                                 grid4d_;
    std::map<unsigned int, GridPt>         grid5d_;
    GridRbPt                               gridRb4d_;
    std::map<unsigned int, GridRbPt>       gridRb5d_;
    std::map<std::vector<int>, GridPt>     gridCoarse4d_;
    std::map<std::vector<int>, GridPt>     gridCoarse5d_;
    unsigned int                           nd_;
    // random number generator
    RngPt                                  rng4d_;
    // object store
    std::vector<ObjInfo>                   object_;
    std::map<std::string, unsigned int>    objectAddress_;
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
        object_[address].type        = &typeid(B);
        object_[address].derivedType = &typeid(T);
        if (MemoryProfiler::stats == &memStats)
        {
            MemoryProfiler::stats = nullptr;
        }
    }
    // object already exists, no error if it is a cache, error otherwise
    else if ((object_[address].storage     != Storage::cache) or 
             (object_[address].storage     != storage)        or
             (object_[address].name        != name)           or
             (object_[address].type        != &typeid(B))     or
             (object_[address].derivedType != &typeid(T)))
    {
        HADRONS_ERROR(Definition, "object '" + name + "' already allocated");
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
                        HADRONS_ERROR(Definition, "object with address " + std::to_string(address) +
                            " cannot be casted to '" + typeName(&typeid(T)) +
                            "' (has type '" + typeName(&typeid(h->get())) + "')");
                    }
                }
            }
            else
            {
                HADRONS_ERROR(Definition, "object with address " + std::to_string(address) +
                            " does not have type '" + typeName(&typeid(B)) +
                            "' (has type '" + getObjectType(address) + "')");
            }
        }
        else
        {
            HADRONS_ERROR(Definition, "object with address " + std::to_string(address) +
                         " is empty");
        }
    }
    else
    {
        HADRONS_ERROR(Definition, "no object with address " + std::to_string(address));
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
        HADRONS_ERROR(Definition, "no object with address " + std::to_string(address));
    }
}

template <typename T>
bool Environment::isObjectOfType(const std::string name) const
{
    return isObjectOfType<T>(getObjectAddress(name));
}

END_HADRONS_NAMESPACE

#endif // Hadrons_Environment_hpp_
