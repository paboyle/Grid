/*
 * MQuark.cc, part of Grid
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

#include <Hadrons/MQuark.hpp>

using namespace Grid;
using namespace Hadrons;

/******************************************************************************
 *                          MQuark implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
MQuark::MQuark(const std::string &name)
: Module(name)
{}

// parse parameters
void MQuark::parseParameters(XmlReader &reader, const std::string &name)
{
    read(reader, name, par_);
}

// dependency relation
std::vector<std::string> MQuark::getInput(void)
{
    return std::vector<std::string>();
}

std::vector<std::string> MQuark::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_5d"};
    
    return out;
}

// allocation //////////////////////////////////////////////////////////////////
void MQuark::allocate(Environment &env)
{
    env.addProp(getName());
    quark_ = env.getProp(getName());
    if (par_.Ls > 1)
    {
        env.addProp(getName() + "_5d", par_.Ls);
        quark5d_ = env.getProp(getName() + "_5d");
    }
}

// execution
void MQuark::execute(Environment &env)
{
    LOG(Message) << "computing quark propagator '" << getName() << "'"
                 << std::endl;
    
    GridCartesian         *g4d   = env.get4dGrid(),
                          *g5d   = env.get5dGrid(par_.Ls);
    GridRedBlackCartesian *gRb4d = env.getRb4dGrid(),
                          *gRb5d = env.getRb5dGrid(par_.Ls);
    LatticeGaugeField     &Umu   = *env.getGauge();
    LatticeFermion        src(g5d); src=zero;
    LatticeFermion        result(g5d); result=zero;
    
    RealD mass=0.1;
    RealD M5=1.8;
    DomainWallFermionR Ddwf(Umu, *g5d, *gRb5d, *g4d, *gRb4d, mass, M5);
    
    ConjugateGradient<LatticeFermion> CG(1.0e-8,10000);
    SchurRedBlackDiagMooeeSolve<LatticeFermion> SchurSolver(CG);
    SchurSolver(Ddwf,src,result);
}
