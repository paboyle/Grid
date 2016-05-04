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
    typedef std::unique_ptr<LatticePropagator>          PropPt;
    typedef std::unique_ptr<LatticeGaugeField>          GaugePt;
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
    void                    addFermionMatrix(const std::string name, FMat *mat);
    FMat *                  getFermionMatrix(const std::string name) const;
    // solvers
    void                    addSolver(const std::string name, Solver s,
                                      const std::string actionName);
    void                    callSolver(const std::string name,
                                       LatticeFermion &sol,
                                       const LatticeFermion &src) const;
    // quark propagators
    void                    createProp(const std::string name,
                                       const unsigned int Ls = 1);
    void                    freeProp(const std::string name);
    bool                    isProp5d(const std::string name) const;
    unsigned int            getPropLs(const std::string name) const;
    LatticePropagator *     getProp(const std::string name) const;
    bool                    propExists(const std::string name) const;
    unsigned int            nProp(void) const;
    // gauge configurations
    void                    createGauge(const std::string name);
    void                    freeGauge(const std::string name);
    LatticeGaugeField *     getGauge(const std::string name) const;
    bool                    gaugeExists(const std::string name) const;
    // random number generator
    void                    setSeed(const std::vector<int> &seed);
    GridParallelRNG *       get4dRng(void) const;
    // general free
    void                    free(const std::string name);
    void                    freeAll(void);
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
    std::map<std::string, PropPt>       prop_;
    std::map<std::string, unsigned int> propSize_;
    std::map<std::string, GaugePt>      gauge_;
};

END_HADRONS_NAMESPACE

#endif // Hadrons_Environment_hpp_
