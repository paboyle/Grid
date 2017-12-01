/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Environment.hpp

Copyright (C) 2015
Copyright (C) 2016

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
#include <Grid/Hadrons/Graph.hpp>

#ifndef SITE_SIZE_TYPE
#define SITE_SIZE_TYPE unsigned int
#endif

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Global environment                                 *
 ******************************************************************************/
// forward declaration of Module
class ModuleBase;

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

class Environment
{
    SINGLETON(Environment);
public:
    typedef SITE_SIZE_TYPE                         Size;
    typedef std::unique_ptr<ModuleBase>            ModPt;
    typedef std::unique_ptr<GridCartesian>         GridPt;
    typedef std::unique_ptr<GridRedBlackCartesian> GridRbPt;
    typedef std::unique_ptr<GridParallelRNG>       RngPt;
    typedef std::unique_ptr<LatticeBase>           LatticePt;
    enum class Storage {object, cache, temporary};
private:
    struct ModuleInfo
    {
        const std::type_info      *type{nullptr};
        std::string               name;
        ModPt                     data{nullptr};
        std::vector<unsigned int> input;
        size_t                    maxAllocated;
    };
    struct ObjInfo
    {
        Size                    size{0};
        Storage                 storage{Storage::object};
        unsigned int            Ls{0};
        const std::type_info    *type{nullptr};
        std::string             name;
        int                     module{-1};
        std::set<unsigned int>  owners, properties;
        std::unique_ptr<Object> data{nullptr};
    };
public:
    // dry run
    void                    dryRun(const bool isDry);
    bool                    isDryRun(void) const;
    void                    memoryProfile(const bool doMemoryProfile);
    bool                    doMemoryProfile(void) const;
    // trajectory number
    void                    setTrajectory(const unsigned int traj);
    unsigned int            getTrajectory(void) const;
    // grids
    void                    createGrid(const unsigned int Ls);
    GridCartesian *         getGrid(const unsigned int Ls = 1) const;
    GridRedBlackCartesian * getRbGrid(const unsigned int Ls = 1) const;
    std::vector<int>        getDim(void) const;
    int                     getDim(const unsigned int mu) const;
    unsigned int            getNd(void) const;
    // random number generator
    void                    setSeed(const std::vector<int> &seed);
    GridParallelRNG *       get4dRng(void) const;
    // module management
    void                    pushModule(ModPt &pt);
    template <typename M>
    void                    createModule(const std::string name);
    template <typename M>
    void                    createModule(const std::string name,
                                         const typename M::Par &par);
    void                    createModule(const std::string name,
                                         const std::string type,
                                         XmlReader &reader);
    unsigned int            getNModule(void) const;
    ModuleBase *            getModule(const unsigned int address) const;
    ModuleBase *            getModule(const std::string name) const;
    template <typename M>
    M *                     getModule(const unsigned int address) const;
    template <typename M>
    M *                     getModule(const std::string name) const;
    unsigned int            getModuleAddress(const std::string name) const;
    std::string             getModuleName(const unsigned int address) const;
    std::string             getModuleType(const unsigned int address) const;
    std::string             getModuleType(const std::string name) const;
    std::string             getModuleNamespace(const unsigned int address) const;
    std::string             getModuleNamespace(const std::string name) const;
    bool                    hasModule(const unsigned int address) const;
    bool                    hasModule(const std::string name) const;
    Graph<unsigned int>     makeModuleGraph(void) const;
    void                    checkGraph(void) const;
    Size                    executeProgram(const std::vector<unsigned int> &p);
    Size                    executeProgram(const std::vector<std::string> &p);
    // general memory management
    void                    addObject(const std::string name,
                                      const int moduleAddress = -1);
    template <typename T, typename P>
    void                    createObject(const std::string name,
                                         const Storage storage,
                                         const unsigned int Ls,
                                         P &&pt);
    template <typename T>
    T *                     getObject(const unsigned int address) const;
    template <typename T>
    T *                     getObject(const std::string name) const;
    unsigned int            getObjectAddress(const std::string name) const;
    std::string             getObjectName(const unsigned int address) const;
    std::string             getObjectType(const unsigned int address) const;
    std::string             getObjectType(const std::string name) const;
    Size                    getObjectSize(const unsigned int address) const;
    Size                    getObjectSize(const std::string name) const;
    unsigned int            getObjectModule(const unsigned int address) const;
    unsigned int            getObjectModule(const std::string name) const;
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
    void                    addOwnership(const unsigned int owner,
                                         const unsigned int property);
    void                    addOwnership(const std::string owner,
                                         const std::string property);
    bool                    hasOwners(const unsigned int address) const;
    bool                    hasOwners(const std::string name) const;
    bool                    freeObject(const unsigned int address);
    bool                    freeObject(const std::string name);
    void                    freeAll(void);
    void                    printContent(void);
private:
    // general
    bool                                   dryRun_{false}, memoryProfile_{false};
    unsigned int                           traj_, locVol_;
    // grids
    std::vector<int>                       dim_;
    GridPt                                 grid4d_;
    std::map<unsigned int, GridPt>         grid5d_;
    GridRbPt                               gridRb4d_;
    std::map<unsigned int, GridRbPt>       gridRb5d_;
    unsigned int                           nd_;
    // random number generator
    RngPt                                  rng4d_;
    // module and related maps
    std::vector<ModuleInfo>                module_;
    std::map<std::string, unsigned int>    moduleAddress_;
    std::string                            currentModule_{""};
    // lattice store
    std::map<unsigned int, LatticePt>      lattice_;
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
    return &objPt_.get();
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
// module management ///////////////////////////////////////////////////////////
template <typename M>
void Environment::createModule(const std::string name)
{
    ModPt pt(new M(name));
    
    pushModule(pt);
}

template <typename M>
void Environment::createModule(const std::string name,
                               const typename M::Par &par)
{
    ModPt pt(new M(name));
    
    static_cast<M *>(pt.get())->setPar(par);
    pushModule(pt);
}

template <typename M>
M * Environment::getModule(const unsigned int address) const
{
    if (auto *pt = dynamic_cast<M *>(getModule(address)))
    {
        return pt;
    }
    else
    {
        HADRON_ERROR("module '" + module_[address].name
                     + "' does not have type " + typeid(M).name()
                     + "(object type: " + getModuleType(address) + ")");
    }
}

template <typename M>
M * Environment::getModule(const std::string name) const
{
    return getModule<M>(getModuleAddress(name));
}

template <typename T, typename P>
void Environment::createObject(const std::string name, 
                               const Environment::Storage storage,
                               const unsigned int Ls,
                               P &&pt)
{
    if (!hasObject(name))
    {
        addObject(name);
    }
    
    unsigned int address = getObjectAddress(name);
    
    if (!object_[address].data)
    {
        MemoryStats  memStats;

        MemoryProfiler::stats    = &memStats;
        object_[address].storage = storage;
        object_[address].Ls      = Ls;
        object_[address].data.reset(new Holder<T>(pt));
        object_[address].size    = memStats.totalAllocated;
        object_[address].type    = &typeid(T);
        MemoryProfiler::stats    = nullptr;
    }
    else
    {
        HADRON_ERROR("object '" + name + "' already allocated");
    }
}

template <typename T>
T * Environment::getObject(const unsigned int address) const
{
    if (hasObject(address))
    {
        if (auto h = dynamic_cast<Holder<T> *>(object_[address].data.get()))
        {
            return h->getPt();
        }
        else
        {
            HADRON_ERROR("object with address " + std::to_string(address) +
                         " does not have type '" + typeName(&typeid(T)) +
                         "' (has type '" + getObjectType(address) + "')");
        }
    }
    else
    {
        HADRON_ERROR("no object with address " + std::to_string(address));
    }
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
        HADRON_ERROR("no object with address " + std::to_string(address));
    }
}

template <typename T>
bool Environment::isObjectOfType(const std::string name) const
{
    return isObjectOfType<T>(getObjectAddress(name));
}

END_HADRONS_NAMESPACE

#endif // Hadrons_Environment_hpp_
