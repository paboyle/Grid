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

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Global environment                                 *
 ******************************************************************************/
class Environment
{
    SINGLETON(Environment);
public:
    typedef FermionOperator<WilsonImplR>                FMat;
    typedef std::function<void(LatticeFermion &,
                               const LatticeFermion &)> Solver;
    typedef std::unique_ptr<GridCartesian>              GridPt;
    typedef std::unique_ptr<GridRedBlackCartesian>      GridRbPt;
    typedef std::unique_ptr<GridParallelRNG>            RngPt;
    typedef std::unique_ptr<FMat>                       FMatPt;
    typedef std::unique_ptr<LatticeBase>                LatticePt;
private:
    struct ObjInfo
    {
        unsigned int size, Ls;
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
    // fermion actions
    void                    addFermionMatrix(const std::string name, FMat *mat);
    FMat *                  getFermionMatrix(const std::string name) const;
    void                    freeFermionMatrix(const std::string name);
    bool                    hasFermionMatrix(const std::string name) const;
    // solvers
    void                    addSolver(const std::string name, Solver s);
    bool                    hasSolver(const std::string name) const;
    void                    setSolverAction(const std::string name,
                                            const std::string actionName);
    std::string             getSolverAction(const std::string name) const;
    void                    callSolver(const std::string name,
                                       LatticeFermion &sol,
                                       const LatticeFermion &src) const;
    // random number generator
    void                    setSeed(const std::vector<int> &seed);
    GridParallelRNG *       get4dRng(void) const;
    // lattice store
    template <typename T>
    T *                     create(const std::string name);
    template <typename T>
    T *                     get(const std::string name) const;
    bool                    hasLattice(const std::string name) const;
    void                    freeLattice(const std::string name);
    template <typename T>
    unsigned int            lattice4dSize(void) const;
    // general memory management
    bool                    hasObject(const std::string name) const;
    void                    registerObject(const std::string name,
                                           const unsigned int size,
                                           const unsigned int Ls = 1);
    template <typename T>
    void                    registerLattice(const std::string name,
                                            const unsigned int Ls = 1);
    unsigned int            getObjectSize(const std::string name) const;
    long unsigned int       getTotalSize(void) const;
    unsigned int            getObjectLs(const std::string name) const;
    bool                    isObject5d(const std::string name) const;
    void                    addOwnership(const std::string owner,
                                         const std::string property);
    bool                    hasOwners(const std::string name) const;
    bool                    freeObject(const std::string name);
    void                    freeAll(void);
private:
    
private:
    bool                                         dryRun_{false};
    unsigned int                                 traj_;
    GridPt                                       grid4d_;
    std::map<unsigned int, GridPt>               grid5d_;
    GridRbPt                                     gridRb4d_;
    std::map<unsigned int, GridRbPt>             gridRb5d_;
    RngPt                                        rng4d_;
    std::map<std::string, ObjInfo>               object_;
    std::map<std::string, LatticePt>             lattice_;
    std::map<std::string, FMatPt>                fMat_;
    std::map<std::string, Solver>                solver_;
    std::map<std::string, std::string>           solverAction_;
    std::map<std::string, std::set<std::string>> owners_;
    std::map<std::string, std::set<std::string>> properties_;
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
    GridCartesian *g = getGrid(getObjectLs(name));
    
    lattice_[name].reset(new T(g));

    return dynamic_cast<T *>(lattice_[name].get());
}

template <typename T>
T * Environment::get(const std::string name) const
{
    if (hasLattice(name))
    {
        if (auto pt = dynamic_cast<T *>(lattice_.at(name).get()))
        {
            return pt;
        }
        else
        {
            HADRON_ERROR("object '" + name + "' does not have type "
                         + typeid(T *).name() + "(object type: "
                         + typeid(lattice_.at(name).get()).name() + ")");
        }
    }
    else
    {
        HADRON_ERROR("no lattice name '" + name + "'");
        
        return nullptr;
    }
}

template <typename T>
void Environment::registerLattice(const std::string name, const unsigned int Ls)
{
    createGrid(Ls);
    registerObject(name, Ls*lattice4dSize<T>());
}

END_HADRONS_NAMESPACE

#endif // Hadrons_Environment_hpp_
