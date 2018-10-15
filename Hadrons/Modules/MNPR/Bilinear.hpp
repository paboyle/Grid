/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MNPR/Bilinear.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Julia Kettle J.R.Kettle-2@sms.ed.ac.uk
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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

#ifndef Hadrons_Bilinear_hpp_
#define Hadrons_Bilinear_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/ModuleFactory.hpp>
//#include <Grid/qcd/utils/PropagatorUtils.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                TBilinear                                       *
        Performs bilinear contractions of the type tr[g5*adj(Sout)*g5*G*Sin]
        Suitable for non exceptional momenta in Rome-Southampton NPR
******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNPR)

class BilinearPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BilinearPar,
                                    std::string,    Sin,
                                    std::string,    Sout,
                                    std::string,    pin,
                                    std::string,    pout,
                                    std::string,    output);
};

template <typename FImpl1, typename FImpl2>
class TBilinear: public Module<BilinearPar>
{
public:
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result, 
                                        std::vector<SpinColourMatrix>, bilinear);
    };
public:
    // constructor
    TBilinear(const std::string name);
    // destructor
    virtual ~TBilinear(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    //LatticeSpinColourMatrix PhaseProps(LatticeSpinColourMatrix S, std::vector<Real> p);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Bilinear, ARG(TBilinear<FIMPL, FIMPL>), MNPR);

/******************************************************************************
 *                           TBilinear implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
TBilinear<FImpl1, FImpl2>::TBilinear(const std::string name)
: Module<BilinearPar>(name)
{}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TBilinear<FImpl1, FImpl2>::setup(void)
{
    //env().template registerLattice<LatticeSpinColourMatrix>(getName());
    //env().template registerObject<SpinColourMatrix>(getName());
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> TBilinear<FImpl1, FImpl2>::getInput(void)
{
    std::vector<std::string> input = {par().Sin, par().Sout};
    
    return input;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TBilinear<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

/*
/////Phase propagators//////////////////////////
template <typename FImpl1, typename FImpl2>
LatticeSpinColourMatrix TBilinear<FImpl1, FImpl2>::PhaseProps(LatticeSpinColourMatrix S, std::vector<Real> p)
{
    GridBase *grid = S._grid;
    LatticeComplex      pdotx(grid),  coor(grid);
    std::vector<int>   latt_size = grid->_fdimensions; 
    Complex             Ci(0.0,1.0);
    pdotx=zero;
    for (unsigned int mu = 0; mu < 4; ++mu)
    {
        Real TwoPiL =  M_PI * 2.0/ latt_size[mu];
        LatticeCoordinate(coor,mu);
        pdotx = pdotx +(TwoPiL * p[mu]) * coor;
    }
    S = S*exp(-Ci*pdotx);
    return S;
}
*/
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TBilinear<FImpl1, FImpl2>::execute(void)
{
/**************************************************************************

Compute the bilinear vertex needed for the NPR.
V(G) = sum_x  [ g5 * adj(S'(x,p2)) * g5 * G * S'(x,p1) ]_{si,sj,ci,cj}
G is one of the 16 gamma vertices [I,gmu,g5,g5gmu,sig(mu,nu)]

        * G
       / \
    p1/   \p2
     /     \
    /       \

Returns a spin-colour matrix, with indices si,sj, ci,cj

Conventions:
p1 - incoming momenta
p2 - outgoing momenta
q = (p1-p2)
**************************************************************************/

    LOG(Message) << "Computing bilinear contractions '" << getName() << "' using"
                 << " momentum '" << par().Sin << "' and '" << par().Sout << "'"
                 << std::endl;
     
    BinaryWriter             writer(par().output);
    

    // Propogators
    LatticeSpinColourMatrix     &Sin = *env().template getObject<LatticeSpinColourMatrix>(par().Sin);
    LatticeSpinColourMatrix     &Sout = *env().template getObject<LatticeSpinColourMatrix>(par().Sout);
    LatticeComplex              pdotxin(env().getGrid()), pdotxout(env().getGrid()), coor(env().getGrid());
    // momentum on legs
    std::vector<Real>           pin  = strToVec<Real>(par().pin), pout = strToVec<Real>(par().pout);
    std::vector<Real>           latt_size(pin.begin(), pin.end()); 
    //bilinears
    LatticeSpinColourMatrix     bilinear_x(env().getGrid());
    SpinColourMatrix            bilinear;
    Gamma                       g5(Gamma::Algebra::Gamma5);
    Result                      result;
    Complex                     Ci(0.0,1.0);

    //

    pdotxin=zero;
    pdotxout=zero;
    for (unsigned int mu = 0; mu < 4; ++mu)
    {
        Real TwoPiL =  M_PI * 2.0/ latt_size[mu];
        LatticeCoordinate(coor,mu);
        pdotxin = pdotxin +(TwoPiL * pin[mu]) * coor;
        pdotxout= pdotxout +(TwoPiL * pout[mu]) * coor;
    }
    Sin = Sin*exp(-Ci*pdotxin); //phase corrections
    Sout = Sout*exp(-Ci*pdotxout);
    
    ////Set up gamma vector//////////////////////////
    std::vector<Gamma> gammavector;
    for( int i=0; i<Gamma::nGamma; i++){
        Gamma::Algebra gam = i;
        gammavector.push_back(Gamma(gam));
    }
    result.bilinear.resize(Gamma::nGamma);
    /////////////////////////////////////////////////
    //LatticeSpinMatrix temp = g5*Sout;
    ////////Form Vertex//////////////////////////////
    for (int i=0; i < Gamma::nGamma; i++){
        bilinear_x = g5*adj(Sout)*g5*gammavector[i]*Sin; 
        result.bilinear[i] = sum(bilinear_x); //sum over lattice sites
    }
    //////////////////////////////////////////////////
    write(writer, par().output, result.bilinear);
    LOG(Message) << "Complete. Writing results to " << par().output << std:: endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_Bilinear_hpp_
