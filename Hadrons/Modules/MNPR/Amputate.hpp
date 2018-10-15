/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MNPR/Amputate.hpp

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

#ifndef Hadrons_Amputate_hpp_
#define Hadrons_Amputate_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Grid/Eigen/LU>
//#include <Grid/qcd/utils/PropagatorUtils.h>
//#include <Grid/serialisation/Serialisation.h>
BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                TAmputate                                       *
        Performs bilinear contractions of the type tr[g5*adj(Sout)*g5*G*Sin]
        Suitable for non exceptional momenta
******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNPR)

class AmputatePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(AmputatePar,
                                    std::string,    Sin, //need to make this a propogator type?
                                    std::string,    Sout, //same
                                    std::string,    vertex,
                                    std::string,    pin,
                                    std::string,    pout,
                                    std::string,    output,
                                    std::string,    input);
};

template <typename FImpl1, typename FImpl2>
class TAmputate: public Module<AmputatePar>
{
public:
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::vector<Complex>, Vamp,
                                        ); 
    };
public:
    // constructor
    TAmputate(const std::string name);
    // destructor
    virtual ~TAmputate(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual SpinColourMatrix invertspincolmat(SpinColourMatrix &scmat);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Amputate, ARG(TAmputate<FIMPL, FIMPL>), MNPR);

/******************************************************************************
 *                           TAmputate implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
TAmputate<FImpl1, FImpl2>::TAmputate(const std::string name)
: Module<AmputatePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> TAmputate<FImpl1, FImpl2>::getInput(void)
{
    std::vector<std::string> input = {par().Sin, par().Sout, par().vertex};
    
    return input;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TAmputate<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> output = {getName()};
    
    
    return output;
}

// Invert spin colour matrix using Eigen
template <typename Fimpl1, typename Fimpl2>
SpinColourMatrix TAmputate<Fimpl1, Fimpl2>::invertspincolmat(SpinColourMatrix &scmat)
{
    Eigen::MatrixXcf scmat_2d(Ns*Nc,Ns*Nc);
    for(int ic=0; ic<Nc; ic++){
    for(int jc=0; jc<Nc; jc++){
        for(int is=0; is<Ns; is++){
        for(int js=0; js<Ns; js++){
            scmat_2d(Ns*ic+is,Ns*jc+js) = scmat()(is,js)(ic,jc);
        }}
    }}      
    Eigen::MatrixXcf scmat_2d_inv = scmat_2d.inverse();
    SpinColourMatrix scmat_inv;
    for(int ic=0; ic<Nc; ic++){
    for(int jc=0; jc<Nc; jc++){
        for(int is=0; is<Ns; is++){
        for(int js=0; js<Ns; js++){
            scmat_inv()(is,js)(ic,jc) = scmat_2d_inv(Ns*ic+is,Ns*jc+js);
        }}
    }}      
    return scmat_inv;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TAmputate<FImpl1, FImpl2>::execute(void)
{
    LOG(Message) << "Computing bilinear amputations '" << getName() << "' using"
                 << " momentum '" << par().Sin << "' and '" << par().Sout << "'"
                 << std::endl;
    BinaryWriter                    writer(par().output);
    PropagatorField1                &Sin = *env().template getObject<PropagatorField1>(par().Sin); //Do these have the phases taken into account?? Don't think so. FIX
    PropagatorField2                &Sout = *env().template getObject<PropagatorField2>(par().Sout);
    std::vector<int>                pin  = strToVec<int>(par().pin), pout = strToVec<int>(par().pout);
    std::vector<Real>               latt_size(pin.begin(), pin.end()); 
    LatticeComplex                  pdotxin(env().getGrid()), pdotxout(env().getGrid()), coor(env().getGrid());
    LOG(Message) << "Propagators set up " << std::endl;
    std::vector<SpinColourMatrix>   vertex; // Let's read from file here
    Gamma                           g5(Gamma::Algebra::Gamma5);
    Result                          result;
    LOG(Message) << "reading file - "  << par().input << std::endl;
    BinaryReader                    reader(par().input); 
    Complex                         Ci(0.0,1.0);

    std::string svertex;
    read(reader,"vertex", vertex);
    LOG(Message) << "vertex read" << std::endl;

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

    SpinColourMatrix Sin_mom = sum(Sin);
    SpinColourMatrix Sout_mom = sum(Sout);
    LOG(Message) << "summed over lattice" << std::endl;
   
    LOG(Message) << "Lattice -> spincolourmatrix conversion" << std::endl;

    SpinColourMatrix Sin_inv = invertspincolmat(Sin_mom);
    SpinColourMatrix Sout_inv = invertspincolmat(Sout_mom);
    LOG(Message) << "Inversions done" << std::endl;

    result.Vamp.resize(Gamma::nGamma/2);
    for( int mu=0; mu < Gamma::nGamma/2; mu++){
        Gamma::Algebra gam = mu;
        result.Vamp[mu] = 1/12.0*trace(adj(Gamma(mu*2+1))*g5*Sout_inv*g5*vertex[mu]*Sin_inv);
        LOG(Message) << "Vamp[" << mu << "] - " << result.Vamp[mu] << std::endl;
        }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_Amputate_hpp_
