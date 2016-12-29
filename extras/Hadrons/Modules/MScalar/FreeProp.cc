#include <Grid/Hadrons/Modules/MScalar/FreeProp.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MScalar;

/******************************************************************************
*                  TFreeProp implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TFreeProp::TFreeProp(const std::string name)
: Module<FreePropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TFreeProp::getInput(void)
{
    std::vector<std::string> in = {par().source};
    
    return in;
}

std::vector<std::string> TFreeProp::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TFreeProp::setup(void)
{
    env().registerLattice<ScalarField>(getName());
}

// execution ///////////////////////////////////////////////////////////////////
void TFreeProp::execute(void)
{
    ScalarField &prop   = *env().createLattice<ScalarField>(getName());
    ScalarField &source = *env().getObject<ScalarField>(par().source);
    ScalarField *momKernel;
    std::string kerName = "_" + getName() + "_momKernel";
    
    if (!env().hasCreatedObject(kerName))
    {
        LOG(Message) << "Caching momentum space free scalar propagator"
                     << "(mass= " << par().mass << ")..." << std::endl;
        momKernel = env().template createLattice<ScalarField>(kerName);
        Scalar<SIMPL>::MomentumSpacePropagator(*momKernel, par().mass);
    }
    else
    {
        momKernel = env().getObject<ScalarField>(kerName);
    }
    LOG(Message) << "Computing free scalar propagator..." << std::endl;
    Scalar<SIMPL>::FreePropagator(source, prop, *momKernel);
}
