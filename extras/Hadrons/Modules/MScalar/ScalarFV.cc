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
	std::string prop0Name = par().scalarProp + "_0";
	std::string propD1Name = par().scalarProp + "_D1";
    std::string propD1D1Name = par().scalarProp + "_D1D1";
    std::string propD2Name = par().scalarProp + "_D2";
    std::vector<std::string> in = {par().source, par().emField, par().scalarProp,
    								prop0Name, propD1Name, propD1D1Name, propD2Name};
    
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
