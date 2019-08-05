/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/Baryon2.hpp

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

#ifndef Hadrons_MContraction_Baryon2_hpp_
#define Hadrons_MContraction_Baryon2_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Grid/qcd/utils/BaryonUtils.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                               Baryon2                                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class Baryon2Par: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(Baryon2Par,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, q3,
                                    std::string, gamma,
                                    std::string, output);
};

template <typename FImpl1, typename FImpl2, typename FImpl3>
class TBaryon2: public Module<Baryon2Par>
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
    TBaryon2(const std::string name);
    // destructor
    virtual ~TBaryon2(void) {};
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

MODULE_REGISTER_TMP(Baryon2, ARG(TBaryon2<FIMPL, FIMPL, FIMPL>), MContraction);

/******************************************************************************
 *                         TBaryon2 implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
TBaryon2<FImpl1, FImpl2, FImpl3>::TBaryon2(const std::string name)
: Module<Baryon2Par>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TBaryon2<FImpl1, FImpl2, FImpl3>::getInput(void)
{
    std::vector<std::string> input = {par().q1, par().q2, par().q3};
    
    return input;
}

template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TBaryon2<FImpl1, FImpl2, FImpl3>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TBaryon2<FImpl1, FImpl2, FImpl3>::setup(void)
{
    envTmpLat(LatticeComplex, "c");
  //  envTmpLat(PropagatorField1, "diquark");
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
void TBaryon2<FImpl1, FImpl2, FImpl3>::execute(void)
{
    LOG(Message) << "Computing baryon contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "', and a diquark formed of ('" << par().q2 << "', and '"
                 << par().q3 << "')" << std::endl;
    
    auto       &q1 = envGet(PropagatorField1, par().q1);
    auto       &q2 = envGet(PropagatorField2, par().q2);
    auto       &q3 = envGet(PropagatorField3, par().q3);
    envGetTmp(LatticeComplex, c);
//    envGetTmp(PropagatorField1, diquark); //ACTUALLY MIX OF 2 AND 3!!!!
    Result     result;
    int nt = env().getDim(Tp);
    result.corr.resize(nt);
    const std::string gamma{ par().gamma };
    std::vector<TComplex> buf;

    const Gamma GammaA{ Gamma::Algebra::Identity };
    const Gamma GammaB{ al };

    LatticeSpinColourMatrix diquark( q1._grid ); // TODO: Felix, I added "q1._grid". I presume this is correct?

    diquark = BaryonUtils<FIMPL>::quarkContract13(q2*GammaB,GammaB*q3);

    //result = trace(GammaA*GammaA * traceColour(q1*traceSpin(diquark))) + 2.0 * trace(GammaA*GammaA*traceColour(q1*diquark));
  //result = trace(q1*diquark); // TODO: Apologies, Felix - compiler errors
  assert( 0 && "TODO: Felix, please fix prior line - compiler errors" );

    sliceSum(c,buf,Tp);

    for (unsigned int t = 0; t < buf.size(); ++t)
    {
        result.corr[t] = TensorRemove(buf[t]);
    }

    saveResult(par().output, "baryon", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_Baryon2_hpp_
