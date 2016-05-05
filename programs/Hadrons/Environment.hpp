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
public:
    // dry run
    void                    dryRun(const bool isDry);
    bool                    isDryRun(void) const;
    // trajectory number
    void                    setTrajectory(const unsigned int traj);
    unsigned int            getTrajectory(void) const;
    // grids
    GridCartesian *         getGrid(const unsigned int Ls = 1) const;
    GridRedBlackCartesian * getRbGrid(const unsigned int Ls = 1) const;
    // fermion actions
    void                    addFermionMatrix(const std::string name, FMat *mat,
                                             const unsigned int size);
    FMat *                  getFermionMatrix(const std::string name) const;
    void                    freeFermionMatrix(const std::string name);
    bool                    hasFermionMatrix(const std::string name) const;
    // solvers
    void                    addSolver(const std::string name, Solver s,
                                      const std::string actionName);
    bool                    hasSolver(const std::string name) const;
    std::string             getSolverAction(const std::string name) const;
    void                    callSolver(const std::string name,
                                       LatticeFermion &sol,
                                       const LatticeFermion &src) const;
    // random number generator
    void                    setSeed(const std::vector<int> &seed);
    GridParallelRNG *       get4dRng(void) const;
    // lattice store
    template <typename T>
    unsigned int            lattice4dSize(void) const;
    template <typename T>
    void                    create(const std::string name,
                                   const unsigned int Ls = 1);
    template <typename T>
    T *                     get(const std::string name) const;
    void                    freeLattice(const std::string name);
    bool                    hasLattice(const std::string name) const;
    bool                    isLattice5d(const std::string name) const;
    unsigned int            getLatticeLs(const std::string name) const;
    // general memory management
    void                    addOwnership(const std::string owner,
                                         const std::string property);
    bool                    hasOwners(const std::string name) const;
    bool                    free(const std::string name);
    void                    freeAll(void);
    unsigned int            getSize(const std::string name) const;
    long unsigned int       getTotalSize(void) const;
private:
    void addSize(const std::string name, const unsigned int size);
private:
    bool                                dryRun_{false};
    unsigned int                        traj_;
    GridPt                              grid4d_;
    std::map<unsigned int, GridPt>      grid5d_;
    GridRbPt                            gridRb4d_;
    std::map<unsigned int, GridRbPt>    gridRb5d_;
    RngPt                               rng4d_;
    std::map<std::string, FMatPt>       fMat_;
    std::map<std::string, Solver>       solver_;
    std::map<std::string, std::string>  solverAction_;
    std::map<std::string, LatticePt>    lattice_;
    std::map<std::string, unsigned int> objectSize_;
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
void Environment::create(const std::string name, const unsigned int Ls)
{
    GridCartesian *g4 = getGrid();
    GridCartesian *g;
    
    if (hasLattice(name))
    {
        HADRON_ERROR("object '" + name + "' already exists");
    }
    if (Ls > 1)
    {
        try
        {
            g = grid5d_.at(Ls).get();
        }
        catch(std::out_of_range &)
        {
            grid5d_[Ls].reset(SpaceTimeGrid::makeFiveDimGrid(Ls, g4));
            gridRb5d_[Ls].reset(SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, g4));
            g = grid5d_[Ls].get();
        }
    }
    else
    {
        g = g4;
    }
    if (!isDryRun())
    {
        lattice_[name].reset(new T(g));
    }
    else
    {
        lattice_[name].reset(nullptr);
    }
    addSize(name, lattice4dSize<T>()*Ls);
}

template <typename T>
T * Environment::get(const std::string name) const
{
    if (hasLattice(name))
    {
        try
        {
            return dynamic_cast<T *>(lattice_.at(name).get());
        }
        catch (std::bad_cast &)
        {
            HADRON_ERROR("object '" + name + "' does not have type "
                         + typeid(T *).name() + "(object type: "
                         + typeid(lattice_.at(name).get()).name() + ")");
        }
    }
    else
    {
        HADRON_ERROR("object '" + name + "' undefined");
        
        return nullptr;
    }
}

END_HADRONS_NAMESPACE

#endif // Hadrons_Environment_hpp_
