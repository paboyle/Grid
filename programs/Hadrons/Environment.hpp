/*
 * Environment.hpp, part of Grid
 *
 * Copyright (C) 2015 Antonin Portelli
 *
 * Grid is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Grid is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Grid.  If not, see <http://www.gnu.org/licenses/>.
 */

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
    typedef std::unique_ptr<GridCartesian>         GridPt;
    typedef std::unique_ptr<GridRedBlackCartesian> GridRbPt;
    typedef std::unique_ptr<GridParallelRNG>       RngPt;
    typedef std::unique_ptr<LatticePropagator>     PropPt;
    typedef std::unique_ptr<LatticeGaugeField>     GaugePt;
public:
    // dry run
    void                    dryRun(const bool isDry);
    bool                    isDryRun(void);
    // grids
    GridCartesian *         get4dGrid(void);
    GridRedBlackCartesian * getRb4dGrid(void);
    GridCartesian *         get5dGrid(const unsigned int Ls);
    GridRedBlackCartesian * getRb5dGrid(const unsigned int Ls);
    // quark propagators
    void                    addProp(const std::string name,
                                    const unsigned int Ls = 1);
    void                    freeProp(const std::string name);
    LatticePropagator *     getProp(const std::string name);
    bool                    propExists(const std::string name);
    unsigned int            nProp(void);
    // gauge configuration
    LatticeGaugeField *     getGauge(void);
    void                    loadUnitGauge(void);
    void                    loadRandomGauge(void);
    // general free
    void                    free(const std::string name);
    void                    freeAll(void);
private:
    bool                                dryRun_{false};
    GridPt                              grid4d_;
    std::map<unsigned int, GridPt>      grid5d_;
    GridRbPt                            gridRb4d_;
    std::map<unsigned int, GridRbPt>    gridRb5d_;
    RngPt                               rng4d_;
    std::map<std::string, PropPt>       prop_;
    std::map<std::string, unsigned int> propSize_;
    GaugePt                             gauge_;
    
};

END_HADRONS_NAMESPACE

#endif // Hadrons_Environment_hpp_
