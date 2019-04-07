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
                                    std::string,    bilinear_lattice,
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
    std::vector<std::string> input = {par().Sin, par().Sout, par().bilinear_lattice};
    
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
    
    // Propogators
    LatticeSpinColourMatrix                     &Sin = *env().template getObject<LatticeSpinColourMatrix>(par().Sin);
    LatticeSpinColourMatrix                     &Sout = *env().template getObject<LatticeSpinColourMatrix>(par().Sout);
    
    // new code
    std::vector<LatticeSpinColourMatrix>        &bilinear_lattice = *env().template getObject<std::vector<LatticeSpinColourMatrix>>(par().bilinear_lattice);
   
    // Averages and inverses
    SpinColourMatrix                            SinAve, SoutAve;
    SpinColourMatrix                            SinInv, SoutInv;
    // g5 for adjoint
    Gamma                                       g5(Gamma::Algebra::Gamma5);
    std::vector<int>   		                    latt_size(env().getGrid()->_fdimensions); 
    Complex                                     Ci(0.0,1.0);
    // momenta and p.x
    std::vector<Real>                           pin  = strToVec<Real>(par().pin), pout = strToVec<Real>(par().pout);
    LatticeComplex                              pdotxin(env().getGrid()), pdotxout(env().getGrid()), coor(env().getGrid());
    //bilinear
    LatticeSpinColourMatrix                     bilinear_mu(env().getGrid());
    std::vector<LatticeSpinColourMatrix>        bilinear;
    //boolean for full or diagonal only basis
    bool                                        fullbasis = par().fullbasis;
    // tensor product results + fourquark result
    LatticeSpinColourSpinColourMatrix           lret(env().getGrid()); 
    Result                                      result;
    
    ////Find volume factor for normalisation/////////////////////////
    RealD vol=1;
    for (unsigned int mu = 0; mu < 4; ++mu)
    { 
      vol=vol*latt_size[mu];
    }
   

    /* 
    pdotxin=zero;
    pdotxout=zero;
    ////Find the phases sun_mu( p_mu x^mu ) for both legs////////////
    for (unsigned int mu = 0; mu < 4; ++mu)
    {
        Real TwoPiL =  M_PI * 2.0/ latt_size[mu];
        LatticeCoordinate(coor,mu);
        pdotxin  = pdotxin  + (TwoPiL * pin[mu] ) * coor;
        pdotxout = pdotxout + (TwoPiL * pout[mu]) * coor;
    }
    // But these are already phased in the bilinear step?
    // Make this all one step?
    // Or let's output the lattice spin colour bilinear from the bilinear code?
    // Then phase the props
    Sin = Sin*exp(-Ci*pdotxin); 
    Sout = Sout*exp(-Ci*pdotxout);
    


    /////////////////////////////////////////////////
    // sum_x ( (g5*Sout^dag*g5) * Gamma * Sin )
    ////////Form Vertex//////////////////////////////
    for (int i=0; i < Gamma::nGamma; i++)
    {
        //bilinear[i] = zero;
        bilinear_mu  = g5*adj(Sout)*g5*gammavector[i]*Sin; 
        bilinear.push_back(bilinear_mu); 
    }
    */
    
    ////Set up gamma vector//////////////////////////
    std::vector<Gamma> gammavector;
    for( int i=0; i<Gamma::nGamma; i++)
    {
        Gamma::Algebra gam = i;
        gammavector.push_back(Gamma(gam));
    }
    // ensure results vector can hold all gamma
    //bilinear.resize(Gamma::nGamma);

    /////////////////////////////////////////////////
    // Debug code to amputate on a single config
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
    std::cout << "Setting up eigen matrices" << std::endl;
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
    
    std::cout << "SVD" << std::endl;
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd_o(Sout_ei,Eigen::ComputeThinU|Eigen::ComputeThinV);
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd_i( Sin_ei,Eigen::ComputeThinU|Eigen::ComputeThinV);

    // Inversion
    std::cout << "inverting" << std::endl;
    Eigen::MatrixXcd Ident = Eigen::MatrixXcd::Identity(12,12);
    Eigen::MatrixXcd SoutInv_ei = Sout_ei.inverse();
    Eigen::MatrixXcd SinInv_ei  = Sin_ei.inverse();
    

    std::cout << "recover spin matrix" << std::endl;
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
    // This is only the gamma scheme for bilinears as a debug check
    /////////////////////////////////////////////////////
    std::cout << "amputate" << std::endl;
    for (int i=0; i < Gamma::nGamma; i++)
    {
      /*
      auto tr = trace (SoutInv * sum(bilinear[i]) * (1.0/vol) * SinInv * gammavector[i] ) / 12.0 ;
      std::cout << "Gamma "<< i << " " << tr <<std::endl;
      */
      auto tr = trace (SoutInv * sum(bilinear_lattice[i]) * (1.0/vol) * SinInv * gammavector[i] ) / 12.0 ;
      std::cout << "Gamma "<< i << " " << tr <<std::endl;
    }


    /*
    ////////////////////////////////
    // Find FourQuark
    if (fullbasis == true)  
    {
        for ( int mu=0; mu<Gamma::nGamma; mu++) 
        for ( int nu=0; nu<Gamma::nGamma; nu++)
        {
            lret = zero;
            tensorprod(lret,bilinear[mu],bilinear[nu]);
            result.fourquark.push_back(sum(lret));
        }
    } 
    else 
    {
        // only diagonal gammas 
        for( int mu=0; mu<Gamma::nGamma; mu++ )
        {
            lret = zero;
            tensorprod(lret,bilinear[mu],bilinear[mu]); //tensor outer product
            result.fourquark.push_back(sum(lret)*(1.0/vol));
        }
    }
    */

    ////////////////////////////////
    // Find FourQuark
    if (fullbasis == true)  
    {
        for ( int mu=0; mu<Gamma::nGamma; mu++) 
        for ( int nu=0; nu<Gamma::nGamma; nu++)
        {
            lret = zero;
            tensorprod(lret,bilinear_lattice[mu],bilinear_lattice[nu]);
            result.fourquark.push_back(sum(lret));
        }
    } 
    else 
    {
        // only diagonal gammas 
        for( int mu=0; mu<Gamma::nGamma; mu++ )
        {
            lret = zero;
            tensorprod(lret,bilinear_lattice[mu],bilinear_lattice[mu]); //tensor outer product
            result.fourquark.push_back(sum(lret)*(1.0/vol));
        }
    }

    saveResult(par().output, "fourquark", result.fourquark); 
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_FourQuark_hpp_
