#include <Grid/Hadrons/Modules/MScalar/FreeProp.hpp>

#define KERNAME "_" + getName() + "_momKernel"

using namespace Grid;
using namespace Hadrons;
using namespace MScalar;

/******************************************************************************
*                        TFreeProp implementation                             *
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
    std::string kerName = KERNAME;
    
    if (!env().hasRegisteredObject(kerName))
    {
        env().registerLattice<ScalarField>(kerName);
    }
    env().registerLattice<ScalarField>(getName());
}

// execution ///////////////////////////////////////////////////////////////////
void TFreeProp::execute(void)
{
    ScalarField &prop   = *env().createLattice<ScalarField>(getName());
    ScalarField &source = *env().getObject<ScalarField>(par().source);
    ScalarField *momKernel;
    std::string kerName = KERNAME;

    if (!env().hasCreatedObject(kerName))
    {
        LOG(Message) << "Caching momentum space free scalar propagator"
                     << " (mass= " << par().mass << ")..." << std::endl;
        momKernel = env().createLattice<ScalarField>(kerName);
        Scalar<SIMPL>::MomentumSpacePropagator(*momKernel, par().mass);
    }
    else
    {
        momKernel = env().getObject<ScalarField>(kerName);
    }
    LOG(Message) << "Computing free scalar propagator..." << std::endl;
    Scalar<SIMPL>::FreePropagator(source, prop, *momKernel);
    
    if (!par().output.empty())
    {
        TextWriter            writer(par().output + "." +
                                     std::to_string(env().getTrajectory()));
        std::vector<TComplex> buf;
        std::vector<Complex>  result;
        
        sliceSum(prop, buf, Tp);
        result.resize(buf.size());
        for (unsigned int t = 0; t < buf.size(); ++t)
        {
            result[t] = TensorRemove(buf[t]);
        }
        write(writer, "prop", result);
    }
}
