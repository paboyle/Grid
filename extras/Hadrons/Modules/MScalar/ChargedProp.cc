#include <Grid/Hadrons/Modules/MScalar/ChargedProp.hpp>
#include <Grid/Hadrons/Modules/MScalar/Scalar.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MScalar;

/******************************************************************************
*                  TChargedProp implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TChargedProp::TChargedProp(const std::string name)
: Module<ChargedPropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TChargedProp::getInput(void)
{
    std::vector<std::string> in = {par().source, par().emField};
    
    return in;
}

std::vector<std::string> TChargedProp::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TChargedProp::setup(void)
{
    freeMomPropName_ = FREEMOMPROP(par().mass);
    shiftedMomPropName_.clear();
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        shiftedMomPropName_.push_back(freeMomPropName_ + "_"
                                      + std::to_string(mu));
    }
    if (!env().hasRegisteredObject(freeMomPropName_))
    {
        env().registerLattice<ScalarField>(freeMomPropName_);
    }
    if (!env().hasRegisteredObject(shiftedMomPropName_[0]))
    {
        for (unsigned int mu = 0; mu < env().getNd(); ++mu)
        {
            env().registerLattice<ScalarField>(shiftedMomPropName_[mu]);
        }
    }
    env().registerLattice<ScalarField>(getName());
    
}

// execution ///////////////////////////////////////////////////////////////////
void TChargedProp::execute(void)
{
    ScalarField &prop   = *env().createLattice<ScalarField>(getName());
    ScalarField &source = *env().getObject<ScalarField>(par().source);
    ScalarField *freeMomProp;
    std::vector<ScalarField *> shiftedMomProp;
    Complex                    ci(0.0,1.0);
    
    if (!env().hasCreatedObject(freeMomPropName_))
    {
        LOG(Message) << "Caching momentum space free scalar propagator"
                     << " (mass= " << par().mass << ")..." << std::endl;
        freeMomProp = env().createLattice<ScalarField>(freeMomPropName_);
        Scalar<SIMPL>::MomentumSpacePropagator(*freeMomProp, par().mass);
    }
    else
    {
        freeMomProp = env().getObject<ScalarField>(freeMomPropName_);
    }
    if (!env().hasCreatedObject(shiftedMomPropName_[0]))
    {
        std::vector<int> &l = env().getGrid()->_fdimensions;
        
        LOG(Message) << "Caching shifted momentum space free scalar propagator"
                     << " (mass= " << par().mass << ")..." << std::endl;
        for (unsigned int mu = 0; mu < env().getNd(); ++mu)
        {
            Real    twoPiL = M_PI*2./l[mu];
            
            shiftedMomProp.push_back(
                env().createLattice<ScalarField>(shiftedMomPropName_[mu]));
            LatticeCoordinate(*(shiftedMomProp[mu]), mu);
            *(shiftedMomProp[mu]) = exp(ci*twoPiL*(*(shiftedMomProp[mu])))
                                    *(*freeMomProp);
        }
    }
    else
    {
        for (unsigned int mu = 0; mu < env().getNd(); ++mu)
        {
            shiftedMomProp.push_back(
                env().getObject<ScalarField>(shiftedMomPropName_[mu]));
        }
    }
    
}
