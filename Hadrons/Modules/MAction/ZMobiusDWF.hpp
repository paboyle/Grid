/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MAction/ZMobiusDWF.hpp

Copyright (C) 2015-2019

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

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#ifndef Hadrons_MAction_ZMobiusDWF_hpp_
#define Hadrons_MAction_ZMobiusDWF_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                      z-Mobius domain-wall fermion action                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MAction)

class ZMobiusDWFPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ZMobiusDWFPar,
                                    std::string                      , gauge,
                                    unsigned int                     , Ls,
                                    double                           , mass,
                                    double                           , M5,
                                    double                           , b,
                                    double                           , c,
                                    std::vector<std::complex<double> >, omega,
                                    std::string                      , boundary,
                                    std::string                      , twist);
};

template <typename FImpl>
class TZMobiusDWF: public Module<ZMobiusDWFPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TZMobiusDWF(const std::string name);
    // destructor
    virtual ~TZMobiusDWF(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ZMobiusDWF, TZMobiusDWF<ZFIMPL>, MAction);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(ZMobiusDWFF, TZMobiusDWF<ZFIMPLF>, MAction);
#endif

/******************************************************************************
 *                     TZMobiusDWF implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TZMobiusDWF<FImpl>::TZMobiusDWF(const std::string name)
: Module<ZMobiusDWFPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TZMobiusDWF<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TZMobiusDWF<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TZMobiusDWF<FImpl>::setup(void)
{
    LOG(Message) << "Setting up z-Mobius domain wall fermion matrix with m= "
                 << par().mass << ", M5= " << par().M5 << ", Ls= " << par().Ls 
                 << ", b= " << par().b << ", c= " << par().c
                 << " using gauge field '" << par().gauge << "'"
                 << std::endl;
    LOG(Message) << "Omegas: " << std::endl;
    for (unsigned int i = 0; i < par().omega.size(); ++i)
    {
        LOG(Message) << "  omega[" << i << "]= " << par().omega[i] << std::endl;
    }

    auto &U    = envGet(GaugeField, par().gauge);
    auto &g4   = *envGetGrid(FermionField);
    auto &grb4 = *envGetRbGrid(FermionField);
    auto &g5   = *envGetGrid(FermionField, par().Ls);
    auto &grb5 = *envGetRbGrid(FermionField, par().Ls);
    auto omega = par().omega;
    typename ZMobiusFermion<FImpl>::ImplParams implParams;
    

    if (!par().boundary.empty())
    {
        implParams.boundary_phases = strToVec<Complex>(par().boundary);
    }
    if (!par().twist.empty())
    {
        implParams.twist_n_2pi_L   = strToVec<Real>(par().twist);
    }
    LOG(Message) << "Fermion boundary conditions: " << implParams.boundary_phases << std::endl;
    LOG(Message) << "Twists: " << implParams.twist_n_2pi_L << std::endl;
    if (implParams.boundary_phases.size() != env().getNd())
    {
        HADRONS_ERROR(Size, "Wrong number of boundary phase");
    }
    if (implParams.twist_n_2pi_L.size() != env().getNd())
    {
        HADRONS_ERROR(Size, "Wrong number of twist");
    }

    assert(par().Ls==omega.size());
    int Ls=par().Ls;
    std::vector<ComplexD> _omega(Ls);
    for(int i=0;i<Ls;i++){
      _omega[i] = omega[i];
    }
    envCreateDerived(FMat, ZMobiusFermion<FImpl>, getName(), par().Ls,
		     U, g5, grb5, g4, grb4, 
		     par().mass, par().M5, 
		     _omega, par().b, par().c, implParams);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TZMobiusDWF<FImpl>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MAction_ZMobiusDWF_hpp_
