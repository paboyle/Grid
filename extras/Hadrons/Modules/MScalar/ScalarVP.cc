#include <Grid/Hadrons/Modules/MScalar/ScalarVP.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MScalar;

/******************************************************************************
*                  TScalarVP implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TScalarVP::TScalarVP(const std::string name)
: Module<ScalarVPPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TScalarVP::getInput(void)
{
    prop0Name_ = par().scalarProp + "_0";
	propD1Name_ = par().scalarProp + "_D1";
    propD1D1Name_ = par().scalarProp + "_D1D1";
    propD2Name_ = par().scalarProp + "_D2";
    std::vector<std::string> in = {par().source, par().emField, par().scalarProp,
    								prop0Name_, propD1Name_, propD1D1Name_,
    								propD2Name_};
    
    return in;
}

std::vector<std::string> TScalarVP::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TScalarVP::setup(void)
{

}

// execution ///////////////////////////////////////////////////////////////////
void TScalarVP::execute(void)
{

}
