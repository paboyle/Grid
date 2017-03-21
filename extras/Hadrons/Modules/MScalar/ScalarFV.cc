#include <Grid/Hadrons/Modules/MScalar/ScalarFV.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MScalar;

/******************************************************************************
*                  TScalarFV implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TScalarFV::TScalarFV(const std::string name)
: Module<ScalarFVPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TScalarFV::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

std::vector<std::string> TScalarFV::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TScalarFV::setup(void)
{

}

// execution ///////////////////////////////////////////////////////////////////
void TScalarFV::execute(void)
{

}
