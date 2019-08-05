/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/Baryon.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Felix Erben <felix.erben@ed.ac.uk>

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

#ifndef Hadrons_MContraction_Baryon_hpp_
#define Hadrons_MContraction_Baryon_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Grid/qcd/utils/BaryonUtils.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                               Baryon                                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class BaryonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BaryonPar,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, q3,
                                    std::string, gamma,
                                    std::string, output);
};

template <typename FImpl1, typename FImpl2, typename FImpl3>
class TBaryon: public Module<BaryonPar>
{
public:
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    FERM_TYPE_ALIASES(FImpl3, 3);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::vector<Complex>, corr);
    };
public:
    // constructor
    TBaryon(const std::string name);
    // destructor
    virtual ~TBaryon(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
    // Which gamma algebra was specified
    Gamma::Algebra  al;
};

MODULE_REGISTER_TMP(Baryon, ARG(TBaryon<FIMPL, FIMPL, FIMPL>), MContraction);

/******************************************************************************
 *                         TBaryon implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
TBaryon<FImpl1, FImpl2, FImpl3>::TBaryon(const std::string name)
: Module<BaryonPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TBaryon<FImpl1, FImpl2, FImpl3>::getInput(void)
{
    std::vector<std::string> input = {par().q1, par().q2, par().q3};
    
    return input;
}

template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TBaryon<FImpl1, FImpl2, FImpl3>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TBaryon<FImpl1, FImpl2, FImpl3>::setup(void)
{
    envTmpLat(LatticeComplex, "c");
    envTmpLat(LatticeComplex, "c1");
    envTmpLat(LatticeComplex, "c2");
    envTmpLat(LatticeComplex, "c3");
    envTmpLat(LatticeComplex, "c4");
    envTmpLat(LatticeComplex, "c5");
    envTmpLat(LatticeComplex, "c6");
    envTmpLat(LatticeComplex, "diquark");
  // Translate the full string naming the desired gamma structure into the one we need to use
  const std::string gamma{ par().gamma };
  int iGamma = 0;
  do
  {
    const char * pGammaName = Gamma::name[iGamma];
    int iLen = 0;
    while( pGammaName[iLen] && pGammaName[iLen] != ' ' )
      iLen++;
    if( !gamma.compare( 0, gamma.size(), pGammaName, iLen ) )
      break;
  }
  while( ++iGamma < Gamma::nGamma );
  if( iGamma >= Gamma::nGamma ) {
    LOG(Message) << "Unrecognised gamma structure \"" << gamma << "\"" << std::endl;
    assert( 0 && "Invalid gamma structure specified" );
  }
  switch( iGamma ) {
    case Gamma::Algebra::GammaX:
      std::cout << "using interpolator C gamma_X" << std::endl;
      al = Gamma::Algebra::GammaZGamma5; //Still hardcoded CgX = i gamma_3 gamma_5
      break;
    case Gamma::Algebra::GammaY:
      std::cout << "using interpolator C gamma_Y" << std::endl;
      al = Gamma::Algebra::GammaT; //Still hardcoded CgX = - gamma_4
      break;
    case Gamma::Algebra::GammaZ:
      std::cout << "using interpolator C gamma_Z" << std::endl;
      al = Gamma::Algebra::GammaXGamma5; //Still hardcoded CgX = i gamma_1 gamma_5
      break;
    default:
    {
      LOG(Message) << "Unsupported gamma structure " << gamma << " = " << iGamma << std::endl;
      assert( 0 && "Unsupported gamma structure" );
      // or you could do something like
      al = static_cast<Gamma::Algebra>( iGamma );
      break;
    }
  }
  LOG(Message) << "Gamma structure " << gamma << " = " << iGamma
               << " translated to " << Gamma::name[al] << std::endl;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TBaryon<FImpl1, FImpl2, FImpl3>::execute(void)
{
    LOG(Message) << "Computing baryon contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "', and a diquark formed of ('" << par().q2 << "', and '"
                 << par().q3 << "')" << std::endl;
    
    auto       &q1 = envGet(PropagatorField1, par().q1);
    auto       &q2 = envGet(PropagatorField2, par().q2);
    auto       &q3 = envGet(PropagatorField3, par().q3);
    envGetTmp(LatticeComplex, c);
    envGetTmp(LatticeComplex, c1);
    envGetTmp(LatticeComplex, c2);
    envGetTmp(LatticeComplex, c3);
    envGetTmp(LatticeComplex, c4);
    envGetTmp(LatticeComplex, c5);
    envGetTmp(LatticeComplex, c6);
    envGetTmp(LatticeComplex, diquark);
    Result     result;
    int nt = env().getDim(Tp);
    result.corr.resize(nt);
    const std::string gamma{ par().gamma };
    std::vector<TComplex> buf;
    
    Result     result1;
    Result     result2;
    Result     result3;
    Result     result4;
    Result     result5;
    Result     result6;
    result1.corr.resize(nt);
    result2.corr.resize(nt);
    result3.corr.resize(nt);
    result4.corr.resize(nt);
    result5.corr.resize(nt);
    result6.corr.resize(nt);
    std::vector<TComplex> buf1;
    std::vector<TComplex> buf2;
    std::vector<TComplex> buf3;
    std::vector<TComplex> buf4;
    std::vector<TComplex> buf5;
    std::vector<TComplex> buf6;

    const Gamma GammaA{ Gamma::Algebra::Identity };
    const Gamma GammaB{ al };

    //BaryonUtils<FIMPL>::ContractBaryons(q1,q2,q3,GammaA,GammaB,c);
    BaryonUtils<FIMPL>::ContractBaryons_debug(q1,q2,q3,GammaA,GammaB,c1,c2,c3,c4,c5,c6,c);

    sliceSum(c,buf,Tp);
    sliceSum(c1,buf1,Tp);
    sliceSum(c2,buf2,Tp);
    sliceSum(c3,buf3,Tp);
    sliceSum(c4,buf4,Tp);
    sliceSum(c5,buf5,Tp);
    sliceSum(c6,buf6,Tp);

    for (unsigned int t = 0; t < buf.size(); ++t)
    {
        result.corr[t] = TensorRemove(buf[t]);
        result1.corr[t] = TensorRemove(buf1[t]);
        result2.corr[t] = TensorRemove(buf2[t]);
        result3.corr[t] = TensorRemove(buf3[t]);
        result4.corr[t] = TensorRemove(buf4[t]);
        result5.corr[t] = TensorRemove(buf5[t]);
        result6.corr[t] = TensorRemove(buf6[t]);
    }

    std::string ostr1{ par().output + "_1"};
    std::string ostr2{ par().output + "_2"};
    std::string ostr3{ par().output + "_3"};
    std::string ostr4{ par().output + "_4"};
    std::string ostr5{ par().output + "_5"};
    std::string ostr6{ par().output + "_6"};
    saveResult(par().output, "baryon", result);
    saveResult(ostr1, "baryon1", result1);
    saveResult(ostr2, "baryon2", result2);
    saveResult(ostr3, "baryon3", result3);
    saveResult(ostr4, "baryon4", result4);
    saveResult(ostr5, "baryon5", result5);
    saveResult(ostr6, "baryon6", result6);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_Baryon_hpp_
