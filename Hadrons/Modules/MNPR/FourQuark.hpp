/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MNPR/FourQuark.hpp

Copyright (C) 2015-2019

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

#ifndef Hadrons_FourQuark_hpp_
#define Hadrons_FourQuark_hpp_

#include <typeinfo>
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Grid/serialisation/Serialisation.h>
BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                TFourQuark                                       *
        Performs fourquark contractions of the type tr[g5*adj(Sout)*g5*G*Sin]
        Suitable for non exceptional momenta
******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNPR)

class FourQuarkPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FourQuarkPar,
                                    std::string,    Sin, //need to make this a propogator type?
                                    std::string,    Sout, //same
                                    std::string,    pin,
                                    std::string,    pout,
                                    bool,           fullbasis,
                                    std::string,    output);
};

template <typename FImpl1, typename FImpl2>
class TFourQuark: public Module<FourQuarkPar>
{
public:
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::vector<SpinColourSpinColourMatrix>, fourquark);
    };
public:
    // constructor
    TFourQuark(const std::string name);
    // destructor
    virtual ~TFourQuark(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void tensorprod(LatticeSpinColourSpinColourMatrix &lret, LatticeSpinColourMatrix a, LatticeSpinColourMatrix b);
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(FourQuark, ARG(TFourQuark<FIMPL, FIMPL>), MNPR);

/******************************************************************************
 *                           TFourQuark implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
TFourQuark<FImpl1, FImpl2>::TFourQuark(const std::string name)
: Module<FourQuarkPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> TFourQuark<FImpl1, FImpl2>::getInput(void)
{
    std::vector<std::string> input = {par().Sin, par().Sout};
    
    return input;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TFourQuark<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> output = {getName()};
    
    return output;
}


template <typename FImpl1, typename FImpl2>
void TFourQuark<FImpl1, FImpl2>::tensorprod(LatticeSpinColourSpinColourMatrix &lret, LatticeSpinColourMatrix a, LatticeSpinColourMatrix b)
{
#if 0
            parallel_for(auto site=lret.begin();site<lret.end();site++) {
                for (int si; si < 4; ++si){
                for(int sj; sj <4; ++sj){
                    for (int ci; ci < 3; ++ci){
                    for (int cj; cj < 3; ++cj){
                        for (int sk; sk < 4; ++sk){
                        for(int sl; sl <4; ++sl){
                            for (int ck; ck < 3; ++ck){
                            for (int cl; cl < 3; ++cl){
                        lret[site]()(si,sj)(ci,cj)(sk,sl)(ck,cl)=a[site]()(si,sj)(ci,cj)*b[site]()(sk,sl)(ck,cl);
                            }}
                        }}
                    }}
                }}
        }
#else 
            // FIXME ; is there a general need for this construct ? In which case we should encapsulate the
            //         below loops in a helper function.
            //LOG(Message) << "sp co mat a is - " << a << std::endl;
            //LOG(Message) << "sp co mat b is - " << b << std::endl;
            parallel_for(auto site=lret.begin();site<lret.end();site++) {
            vTComplex left;
                for(int si=0; si < Ns; ++si){
                for(int sj=0; sj < Ns; ++sj){
                    for (int ci=0; ci < Nc; ++ci){
                    for (int cj=0; cj < Nc; ++cj){
                      //LOG(Message) << "si, sj, ci, cj -  " << si << ", " << sj  << ", "<< ci  << ", "<< cj << std::endl;
                      left()()() = a[site]()(si,sj)(ci,cj);
                      //LOG(Message) << left << std::endl;
                      lret[site]()(si,sj)(ci,cj)=left()*b[site]();
                    }}
                }}
            }
#endif      
}





// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TFourQuark<FImpl1, FImpl2>::setup(void)
{
    envCreateLat(LatticeSpinColourMatrix, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TFourQuark<FImpl1, FImpl2>::execute(void)
{

/*********************************************************************************

TFourQuark : Creates the four quark vertex required for the NPR of four-quark ops

V_{Gamma_1,Gamma_2} = sum_x [ ( g5 * adj(S'(x,p2)) * g5 * G1 * S'(x,p1) )_ci,cj;si,sj x ( g5 * adj(S'(x,p2)) * g5 * G2 S'(x,p1) )_ck,cl;sk,cl ]

Create a bilinear vertex for G1 and G2  the spin and colour indices are kept free. Where there are 16 potential Gs.
We then find the outer product of V1 and V2, keeping the spin and colour indices uncontracted
Then this is summed over the lattice coordinate
Result is a SpinColourSpinColourMatrix - with 4 colour and 4 spin indices. 
We have up to 256 of these including the offdiag (G1 != G2).

        \         /
         \p1   p1/
          \     /
           \   /
         G1 * * G2
           /   \
          /     \
         /p2   p2\
        /         \

*********************************************************************************/




    LOG(Message) << "Computing fourquark contractions '" << getName() << "' using"
                 << " momentum '" << par().Sin << "' and '" << par().Sout << "'"
                 << std::endl;
    
    BinaryWriter             writer(par().output);
    
    PropagatorField1                            &Sin = *env().template getObject<PropagatorField1>(par().Sin);
    PropagatorField2                            &Sout = *env().template getObject<PropagatorField2>(par().Sout);
    std::vector<Real>                           pin  = strToVec<Real>(par().pin), pout = strToVec<Real>(par().pout);
    bool                                        fullbasis = par().fullbasis;
    Gamma                                       g5(Gamma::Algebra::Gamma5);
    Result                                      result;
    std::vector<Real>                           latt_size(pin.begin(), pin.end());
    LatticeComplex                              pdotxin(env().getGrid()), pdotxout(env().getGrid()), coor(env().getGrid());
    LatticeSpinColourMatrix                     bilinear_mu(env().getGrid()), bilinear_nu(env().getGrid());
    LatticeSpinColourSpinColourMatrix           lret(env().getGrid()); 
    Complex                         Ci(0.0,1.0);

    //Phase propagators
    //Sin = Grid::QCD::PropUtils::PhaseProps(Sin,pin);
    //Sout = Grid::QCD::PropUtils::PhaseProps(Sout,pout);
    
    //find p.x for in and out so phase can be accounted for in propagators
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


    //Set up Gammas 
    std::vector<Gamma> gammavector;
     for( int i=1; i<Gamma::nGamma; i+=2){
         Gamma::Algebra gam = i;
         gammavector.push_back(Gamma(gam));
       }
    
    lret = zero;
    if (fullbasis == true){ // all combinations of mu and nu
        result.fourquark.resize(Gamma::nGamma/2*Gamma::nGamma/2);
        for( int mu=0; mu<Gamma::nGamma/2; mu++){ 
            bilinear_mu = g5*adj(Sout)*g5*gammavector[mu]*Sin;
            for ( int nu=0; nu<Gamma::nGamma; nu++){
                LatticeSpinColourMatrix     bilinear_nu(env().getGrid());
                bilinear_nu = g5*adj(Sout)*g5*gammavector[nu]*Sin;
                LOG(Message) << "bilinear_nu for nu = " << nu << " is - " << bilinear_mu << std::endl;
                result.fourquark[mu*Gamma::nGamma/2 + nu] = zero;
                tensorprod(lret,bilinear_mu,bilinear_nu);
                result.fourquark[mu*Gamma::nGamma/2 + nu] = sum(lret);
            }
        }
    } else {
        result.fourquark.resize(Gamma::nGamma/2);
        for ( int mu=0; mu<1; mu++){
        //for( int mu=0; mu<Gamma::nGamma/2; mu++ ){
            bilinear_mu = g5*adj(Sout)*g5*gammavector[mu]*Sin;
            //LOG(Message) << "bilinear_mu for mu = " << mu << " is - " << bilinear_mu << std::endl;
            result.fourquark[mu] = zero;
            tensorprod(lret,bilinear_mu,bilinear_mu); //tensor outer product
            result.fourquark[mu] = sum(lret);
        }
    }
    write(writer, "fourquark", result.fourquark);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_FourQuark_hpp_
