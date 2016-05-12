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

#include <Hadrons/Global.hpp>
#include <Hadrons/Graph.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Global environment                                 *
 ******************************************************************************/
// forward declaration of Module
class ModuleBase;

class Environment
{
    SINGLETON(Environment);
public:
    typedef std::unique_ptr<ModuleBase>                 ModPt;
    typedef std::unique_ptr<GridCartesian>              GridPt;
    typedef std::unique_ptr<GridRedBlackCartesian>      GridRbPt;
    typedef FermionOperator<WilsonImplR>                FMat;
    typedef std::unique_ptr<FMat>                       FMatPt;
    typedef std::function<void(LatticeFermion &,
                               const LatticeFermion &)> Solver;
    typedef std::unique_ptr<GridParallelRNG>            RngPt;
    typedef std::unique_ptr<LatticeBase>                LatticePt;
private:
    struct ObjInfo
    {
        unsigned int size{0}, Ls{0};
        bool         isRegistered{false};
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
    void                    createModule(const std::string name,
                                         const std::string type,
                                         XmlReader &reader);
    ModuleBase *            getModule(const unsigned int address) const;
    ModuleBase *            getModule(const std::string name) const;
    unsigned int            getModuleAddress(const std::string name) const;
    std::string             getModuleName(const unsigned int address) const;
    std::string             getModuleType(const unsigned int address) const;
    std::string             getModuleType(const std::string name) const;
    bool                    hasModule(const unsigned int address) const;
    bool                    hasModule(const std::string name) const;
    Graph<unsigned int>     makeModuleGraph(void) const;
    unsigned int            executeProgram(const std::vector<unsigned int> &p);
    unsigned int            executeProgram(const std::vector<std::string> &p);
    // lattice store
    template <typename T>
    T *                     create(const std::string name);
    template <typename T>
    T *                     get(const std::string name) const;
    bool                    hasLattice(const unsigned int address) const;
    bool                    hasLattice(const std::string name) const;
    void                    freeLattice(const unsigned int address);
    void                    freeLattice(const std::string name);
    template <typename T>
    unsigned int            lattice4dSize(void) const;
    // fermion actions
    void                    addFermionMatrix(const std::string name, FMat *mat);
    FMat *                  getFermionMatrix(const std::string name) const;
    bool                    hasFermionMatrix(const unsigned int address) const;
    bool                    hasFermionMatrix(const std::string name) const;
    void                    freeFermionMatrix(const unsigned int address);
    void                    freeFermionMatrix(const std::string name);
    // solvers
    void                    addSolver(const std::string name, Solver s);
    bool                    hasSolver(const unsigned int address) const;
    bool                    hasSolver(const std::string name) const;
    void                    setSolverAction(const std::string name,
                                            const std::string actionName);
    std::string             getSolverAction(const std::string name) const;
    void                    callSolver(const std::string name,
                                       LatticeFermion &sol,
                                       const LatticeFermion &src) const;
    // general memory management
    void                    registerObject(const unsigned int address,
                                           const unsigned int size,
                                           const unsigned int Ls = 1);
    void                    registerObject(const std::string name,
                                           const unsigned int size,
                                           const unsigned int Ls = 1);
    template <typename T>
    void                    registerLattice(const unsigned int address,
                                            const unsigned int Ls = 1);
    template <typename T>
    void                    registerLattice(const std::string name,
                                            const unsigned int Ls = 1);
    unsigned int            getObjectAddress(const std::string name) const;
    std::string             getObjectName(const unsigned int address) const;
    unsigned int            getObjectSize(const unsigned int address) const;
    unsigned int            getObjectSize(const std::string name) const;
    unsigned int            getObjectLs(const unsigned int address) const;
    unsigned int            getObjectLs(const std::string name) const;
    bool                    hasObject(const unsigned int address) const;
    bool                    hasObject(const std::string name) const;
    bool                    hasRegisteredObject(const unsigned int address) const;
    bool                    hasRegisteredObject(const std::string name) const;
    bool                    isObject5d(const unsigned int address) const;
    bool                    isObject5d(const std::string name) const;
    long unsigned int       getTotalSize(void) const;
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
    std::vector<ModPt>                     module_;
    std::vector<std::string>               moduleType_;
    std::vector<std::string>               moduleName_;
    std::map<std::string, unsigned int>    moduleAddress_;
    std::vector<std::vector<unsigned int>> moduleInput_;
    // lattice store
    std::map<unsigned int, LatticePt>      lattice_;
    // fermion matrix store
    std::map<unsigned int, FMatPt>         fMat_;
    // solver store & solver/action map
    std::map<unsigned int, Solver>         solver_;
    std::map<std::string, std::string>     solverAction_;
    // object register
    std::vector<ObjInfo>                   object_;
    std::vector<std::string>               objectName_;
    std::map<std::string, unsigned int>    objectAddress_;
    std::vector<int>                       objectModule_;
    std::vector<std::set<unsigned int>>    owners_;
    std::vector<std::set<unsigned int>>    properties_;
};

/******************************************************************************
 *                        template implementation                             *
 ******************************************************************************/
template <typename T>
unsigned int Environment::lattice4dSize(void) const
{
    return sizeof(typename T::vector_object)/getGrid()->Nsimd();
}

template <typename T>
T * Environment::create(const std::string name)
{
    auto           i = getObjectAddress(name);
    GridCartesian *g = getGrid(getObjectLs(i));
    
    lattice_[i].reset(new T(g));

    return dynamic_cast<T *>(lattice_[i].get());
}

template <typename T>
T * Environment::get(const std::string name) const
{
    if (hasLattice(name))
    {
        auto i = getObjectAddress(name);
        
        if (auto pt = dynamic_cast<T *>(lattice_.at(i).get()))
        {
            return pt;
        }
        else
        {
            HADRON_ERROR("object '" + name + "' does not have type "
                         + typeName<T>() + "(object type: "
                         + typeName(*lattice_.at(i).get()) + ")");
        }
    }
    else
    {
        HADRON_ERROR("no lattice with name '" + name + "'");
        
        return nullptr;
    }
}

template <typename T>
void Environment::registerLattice(const unsigned int address,
                                  const unsigned int Ls)
{
    createGrid(Ls);
    registerObject(address, Ls*lattice4dSize<T>());
}

template <typename T>
void Environment::registerLattice(const std::string name, const unsigned int Ls)
{
    createGrid(Ls);
    registerObject(name, Ls*lattice4dSize<T>());
}

END_HADRONS_NAMESPACE

#endif // Hadrons_Environment_hpp_
