/*************************************************************************************
Grid physics library, www.github.com/paboyle/Grid 
Source file: Hadrons/Modules/MNPR/Bilinear.hpp
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

#ifndef Hadrons_Bilinear_hpp_
#define Hadrons_Bilinear_hpp_


#include <typeinfo>
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Grid/serialisation/Serialisation.h>
#include <Grid/Eigen/Core>
#include <Grid/Eigen/Dense>
#include <Grid/Eigen/SVD>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                TBilinear                                       *
        Performs bilinear contractions of the type tr[g5*adj(Sout)*g5*G*Sin]
        Suitable for non exceptional momenta in Rome-Southampton Bilinear
        Also outputs the four quark contractions
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
                                    bool,           fullbasis,
                                    std::string,    fourQuark_output,
                                    std::string,    propInOutput,
                                    std::string,    propOutOutput,
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
                                        std::vector<SpinColourMatrix>,  bilinear,
                                        std::vector<SpinColourSpinColourMatrix>, fourquark
                                       );
    };
public:
    // constructor
    TBilinear(const std::string name);
    // destructor
    virtual ~TBilinear(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void tensorprod(LatticeSpinColourSpinColourMatrix &lret, LatticeSpinColourMatrix a, LatticeSpinColourMatrix b);
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
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> TBilinear<FImpl1, FImpl2>::getInput(void)
{
    std::vector<std::string> input = {par().Sin, par().Sout};
    
    return input;
}



template <typename FImpl1, typename FImpl2>
void TBilinear<FImpl1, FImpl2>::tensorprod(LatticeSpinColourSpinColourMatrix &lret, LatticeSpinColourMatrix a, LatticeSpinColourMatrix b)
{
            // Tensor product of 2 Lattice Spin Colour Matrices
            parallel_for(auto site=lret.begin();site<lret.end();site++) {
            vTComplex left;
                for(int si=0; si < Ns; ++si){
                for(int sj=0; sj < Ns; ++sj){
                    for (int ci=0; ci < Nc; ++ci){
                    for (int cj=0; cj < Nc; ++cj){
                      left()()() = a[site]()(si,sj)(ci,cj);
                      lret[site]()(si,sj)(ci,cj)=left()*b[site]();
                    }}
                }}
            }
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TBilinear<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    return out;
}
    

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TBilinear<FImpl1, FImpl2>::execute(void)
{
/**************************************************************************
Compute the bilinear vertex needed for the Bilinear.
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

Also computes and saves the four quarks
V4q(G) = sum_x ( [ g5 * adj(S'(x,p2)) * g5 * G * S'(x,p1) ]_{si,sj,ci,cj} [ g5 * adj(S'(x,p2)) * g5 * G * S'(x,p1) ]_{sk,sl,ck,cl} )
**************************************************************************/

    LOG(Message) << "Computing bilinear and fourquark contractions '" << getName() << "' using"
                 << " momentum '" << par().pin << "' and '" << par().pout << "'"
                 << std::endl;
     
    // Propogators
    LatticeSpinColourMatrix                     &Sin_unphased = *env().template getObject<LatticeSpinColourMatrix>(par().Sin);
    LatticeSpinColourMatrix                     &Sout_unphased = *env().template getObject<LatticeSpinColourMatrix>(par().Sout);
    // Averages and inverses
    SpinColourMatrix                            SinAve, SoutAve;
    SpinColourMatrix                            SinInv, SoutInv;
    
    LatticeComplex                              pdotxin(env().getGrid()), pdotxout(env().getGrid()), coor(env().getGrid());
    Complex                                     Ci(0.0,1.0);
    std::vector<Real>                           pin  = strToVec<Real>(par().pin), pout = strToVec<Real>(par().pout);
    std::vector<int>   		                    latt_size(env().getGrid()->_fdimensions); 
    
    bool                                        fullbasis = par().fullbasis;

    //bilinear_x holds bilinear before summing
    LatticeSpinColourMatrix                     bilinear_x(env().getGrid());
    // vector to hold all lattice bilinear vertices required for 4q calculation
    std::vector<LatticeSpinColourMatrix>        bilinear_x_allGamma;
    SpinColourMatrix                            bilinear;
    SpinColourSpinColourMatrix                  fourquark;
    LatticeSpinColourSpinColourMatrix           lret(env().getGrid()); 

    Gamma                                       g5(Gamma::Algebra::Gamma5);
    Result                                      result;

    

    ////Find volume factor for normalisation/////////////////////////
    RealD vol=1;
    for (unsigned int mu = 0; mu < 4; ++mu)
    { 
      vol=vol*latt_size[mu];
    }
    
    ////Find the phases sun_mu( p_mu x^mu ) for both legs////////////
    pdotxin=zero;
    pdotxout=zero;
    for (unsigned int mu = 0; mu < 4; ++mu)
    {
        Real TwoPiL =  M_PI * 2.0/ latt_size[mu];
        LatticeCoordinate(coor,mu);
        pdotxin = pdotxin  + (TwoPiL * pin[mu] ) * coor;
        pdotxout= pdotxout + (TwoPiL * pout[mu]) * coor;
    }
    // Then phase the props 
    // new variable so code can be called repeatedly
    LatticeSpinColourMatrix Sin = Sin_unphased*exp(-Ci*pdotxin); 
    LatticeSpinColourMatrix Sout = Sout_unphased*exp(-Ci*pdotxout);
    
    ////Set up gamma vector//////////////////////////
    std::vector<Gamma> gammavector;
    for( int i=0; i<Gamma::nGamma; i++)
    {
        Gamma::Algebra gam = i;
        gammavector.push_back(Gamma(gam));
    }
    // ensure results vector can hold all gamma
    result.bilinear.resize(Gamma::nGamma);

    /////////////////////////////////////////////////
    // sum_x ( (g5*Sout^dag*g5) * Gamma * Sin )
    ////////Form Vertex//////////////////////////////
    for (int i=0; i < Gamma::nGamma; i++)
    {
        bilinear_x = g5*adj(Sout)*g5*gammavector[i]*Sin;
        bilinear_x_allGamma.push_back(bilinear_x); 
        result.bilinear[i] = sum(bilinear_x); 
        result.bilinear[i] = result.bilinear[i]*(1.0/vol);
    }
    /////////////////////////////////////////////////
    // amputate on a single config and print as a check
    /////////////////////////////////////////////////
    
    ////Find the propogator averages ( normalising by vol factor )
    SinAve = sum(Sin); SinAve = SinAve*(1.0/vol); 
    SoutAve= sum(Sout);SoutAve= SoutAve*(1.0/vol);
    SoutAve= g5*adj(SoutAve)*g5; // quark flow reversal

    ////////////////////////////////////////////////////
    // Invert the volume averaged propogators
    // Rewrite as 12x12 matrix for inverting 
    // then restructure as spin-colour matrix after
    ////////////////////////////////////////////////////
    Eigen::MatrixXcd Sout_ei = Eigen::MatrixXcd::Zero(12,12);
    Eigen::MatrixXcd Sin_ei  = Eigen::MatrixXcd::Zero(12,12);
    for(int c1=0;c1<Nc;c1++)
    for(int c2=0;c2<Nc;c2++)
    for(int s1=0;s1<Ns;s1++)
    for(int s2=0;s2<Ns;s2++)
    {
        int i1=s1*Nc+c1;
        int i2=s2*Nc+c2;
        Sout_ei(i1,i2) = SoutAve()(s1,s2)(c1,c2);
        Sin_ei(i1,i2)  =  SinAve()(s1,s2)(c1,c2);
    }
    
    // SVD
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd_o(Sout_ei,Eigen::ComputeThinU|Eigen::ComputeThinV);
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd_i( Sin_ei,Eigen::ComputeThinU|Eigen::ComputeThinV);

    // Inversion
    Eigen::MatrixXcd Ident = Eigen::MatrixXcd::Identity(12,12);
    Eigen::MatrixXcd SoutInv_ei = Sout_ei.inverse();
    Eigen::MatrixXcd SinInv_ei  = Sin_ei.inverse();
    
    // recover Spin-colour matrix form
    for(int c1=0;c1<Nc;c1++)
    for(int c2=0;c2<Nc;c2++)
    for(int s1=0;s1<Ns;s1++)
    for(int s2=0;s2<Ns;s2++)
    {
        int i1=s1*Nc+c1;
        int i2=s2*Nc+c2;
        SoutInv()(s1,s2)(c1,c2) = SoutInv_ei(i1,i2);
        SinInv()(s1,s2)(c1,c2)  = SinInv_ei(i1,i2) ;
    }
    
    /////////////////////////////////////////////////////
    // Amputate then project the bilinear results for every gamma
    // Tr ( Sout^-1 * Vertex_i * Sin^-1 * Gamma_i ) / 12
    // By looking at the single config the volume average approximates the gauge average
    // This is only the gamma scheme. Need qslash scheme too
    /////////////////////////////////////////////////////
    for (int i=0; i < Gamma::nGamma; i++)
    {
      auto tr = trace (SoutInv * result.bilinear[i] * SinInv * gammavector[i] ) *( 1.0 / 12.0) ;
      std::cout << "Gamma "<< i << " " << tr <<std::endl;
    }
    /////////////////////////////////////////////////////
    // Find Lambda_S and Lambda_A and print to screen
    /////////////////////////////////////////////////////
    {
      int S = Gamma::Algebra::Identity;
      int AX= Gamma::Algebra::GammaXGamma5;
      int AY= Gamma::Algebra::GammaYGamma5;
      int AZ= Gamma::Algebra::GammaZGamma5;
      int AT= Gamma::Algebra::GammaTGamma5;
      auto tr_S = trace (SoutInv * result.bilinear[S] * SinInv * gammavector[S] ) * (1.0 / 12.0);
      auto tr_A = trace (SoutInv * result.bilinear[AX]* SinInv * gammavector[AX] ) * (1.0 / 48.0);
           tr_A+= trace (SoutInv * result.bilinear[AY]* SinInv * gammavector[AY] ) * (1.0 / 48.0);
           tr_A+= trace (SoutInv * result.bilinear[AZ]* SinInv * gammavector[AZ] ) * (1.0 / 48.0);
           tr_A+= trace (SoutInv * result.bilinear[AT]* SinInv * gammavector[AT] ) * (1.0 / 48.0);

      std::cout << "Lambda_S "<< " " << tr_S <<std::endl;
      std::cout << "Lambda_A "<< " " << tr_A <<std::endl;
      std::cout << "Lambda_S/Lambda_A "<< " " << tr_S/tr_A <<std::endl;
      std::cout << "Lambda_A/Lambda_S "<< " " << tr_A/tr_S <<std::endl;
	   
    }


    ////////////////////////////////
    // Find FourQuark
    if (fullbasis == true)  
    {
        // Compute all vertes Gamma_i Gamma_j  eg  VA or SP
        for ( int mu=0; mu<Gamma::nGamma; mu++) 
        for ( int nu=0; nu<Gamma::nGamma; nu++)
        {
            lret = zero;
            tensorprod(lret,bilinear_x_allGamma[mu],bilinear_x_allGamma[nu]);
            result.fourquark.push_back(sum(lret));
        }
    } 
    else 
    {
        // only diagonal gammas 
        //  only VV, AA etc
        for( int mu=0; mu<Gamma::nGamma; mu++ )
        {
            lret = zero;
            tensorprod(lret,bilinear_x_allGamma[mu],bilinear_x_allGamma[mu]); //tensor outer product
            result.fourquark.push_back(sum(lret)*(1.0/vol));
        }
    }


    /////////////////////////////////////////////////////
    // write the results, props, bilinear and fourquark
    saveResult(par().propInOutput , "SinAve",   SinAve  );
    saveResult(par().propOutOutput, "SoutAve",  SoutAve );
    saveResult(par().output, "bilinear", result.bilinear);
    saveResult(par().fourQuark_output, "fourquark", result.fourquark);
    LOG(Message) << "Complete. Writing results to " << par().output << std:: endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_Bilinear_hpp_
