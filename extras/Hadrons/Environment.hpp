/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/Environment.hpp

Copyright (C) 2015

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

See the full license in the file "LICENSE" in the top level distribution 
directory.
*******************************************************************************/

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
private:
    struct ModuleInfo
    {
        const std::type_info      *type{nullptr};
        std::string               name;
        ModPt                     data{nullptr};
        std::vector<unsigned int> input;
    };
    struct ObjInfo
    {
        Size                    size{0};
        unsigned int            Ls{0};
        bool                    isRegistered{false};
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
    // trajectory number
    void                    setTrajectory(const unsigned int traj);
    unsigned int            getTrajectory(void) const;
    // grids
    void                    createGrid(const unsigned int Ls);
    GridCartesian *         getGrid(const unsigned int Ls = 1) const;
    GridRedBlackCartesian * getRbGrid(const unsigned int Ls = 1) const;
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
    bool                    hasModule(const unsigned int address) const;
    bool                    hasModule(const std::string name) const;
    Graph<unsigned int>     makeModuleGraph(void) const;
    Size                    executeProgram(const std::vector<unsigned int> &p);
    Size                    executeProgram(const std::vector<std::string> &p);
    // general memory management
    void                    addObject(const std::string name,
                                      const int moduleAddress);
    void                    registerObject(const unsigned int address,
                                           const unsigned int size,
                                           const unsigned int Ls = 1);
    void                    registerObject(const std::string name,
                                           const unsigned int size,
                                           const unsigned int Ls = 1);
    template <typename T>
    unsigned int            lattice4dSize(void) const;
    template <typename T>
    void                    registerLattice(const unsigned int address,
                                            const unsigned int Ls = 1);
    template <typename T>
    void                    registerLattice(const std::string name,
                                            const unsigned int Ls = 1);
    template <typename T>
    void                    setObject(const unsigned int address, T *object);
    template <typename T>
    void                    setObject(const std::string name, T *object);
    template <typename T>
    T *                     getObject(const unsigned int address) const;
    template <typename T>
    T *                     getObject(const std::string name) const;
    template <typename T>
    T *                     createLattice(const unsigned int address);
    template <typename T>
    T *                     createLattice(const std::string name);
    unsigned int            getObjectAddress(const std::string name) const;
    std::string             getObjectName(const unsigned int address) const;
    std::string             getObjectType(const unsigned int address) const;
    std::string             getObjectType(const std::string name) const;
    Size                    getObjectSize(const unsigned int address) const;
    Size                    getObjectSize(const std::string name) const;
    unsigned int            getObjectLs(const unsigned int address) const;
    unsigned int            getObjectLs(const std::string name) const;
    bool                    hasObject(const unsigned int address) const;
    bool                    hasObject(const std::string name) const;
    bool                    hasRegisteredObject(const unsigned int address) const;
    bool                    hasRegisteredObject(const std::string name) const;
    bool                    isObject5d(const unsigned int address) const;
    bool                    isObject5d(const std::string name) const;
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
    bool                                   dryRun_{false};
    unsigned int                           traj_, locVol_;
    // grids
    GridPt                                 grid4d_;
    std::map<unsigned int, GridPt>         grid5d_;
    GridRbPt                               gridRb4d_;
    std::map<unsigned int, GridRbPt>       gridRb5d_;
    // random number generator
    RngPt                                  rng4d_;
    // module and related maps
    std::vector<ModuleInfo>                module_;
    std::map<std::string, unsigned int>    moduleAddress_;
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

template <typename T>
unsigned int Environment::lattice4dSize(void) const
{
    return sizeof(typename T::vector_object)/getGrid()->Nsimd();
}

template <typename T>
void Environment::registerLattice(const unsigned int address,
                                  const unsigned int Ls)
{
    createGrid(Ls);
    registerObject(address, Ls*lattice4dSize<T>(), Ls);
}

template <typename T>
void Environment::registerLattice(const std::string name, const unsigned int Ls)
{
    createGrid(Ls);
    registerObject(name, Ls*lattice4dSize<T>(), Ls);
}

template <typename T>
void Environment::setObject(const unsigned int address, T *object)
{
    if (hasRegisteredObject(address))
    {
        object_[address].data.reset(new Holder<T>(object));
        object_[address].type = &typeid(T);
    }
    else if (hasObject(address))
    {
        HADRON_ERROR("object with address " + std::to_string(address) +
                     " exists but is not registered");
    }
    else
    {
        HADRON_ERROR("no object with address " + std::to_string(address));
    }
}

template <typename T>
void Environment::setObject(const std::string name, T *object)
{
    setObject(getObjectAddress(name), object);
}

template <typename T>
T * Environment::getObject(const unsigned int address) const
{
    if (hasRegisteredObject(address))
    {
        if (auto h = dynamic_cast<Holder<T> *>(object_[address].data.get()))
        {
            return h->getPt();
        }
        else
        {
            HADRON_ERROR("object with address " + std::to_string(address) +
                         " does not have type '" + typeid(T).name() +
                         "' (has type '" + getObjectType(address) + "')");
        }
    }
    else if (hasObject(address))
    {
        HADRON_ERROR("object with address " + std::to_string(address) +
                     " exists but is not registered");
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
T * Environment::createLattice(const unsigned int address)
{
    GridCartesian *g = getGrid(getObjectLs(address));
    
    setObject(address, new T(g));
    
    return getObject<T>(address);
}

template <typename T>
T * Environment::createLattice(const std::string name)
{
    return createLattice<T>(getObjectAddress(name));
}

END_HADRONS_NAMESPACE

#endif // Hadrons_Environment_hpp_
