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
	phaseName_.clear();
	for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        phaseName_.push_back("_shiftphase_" + std::to_string(mu));
    }

}

// execution ///////////////////////////////////////////////////////////////////
void TScalarVP::execute(void)
{
	for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        phase_.push_back(env().getObject<ScalarField>(phaseName_[mu]));
    }

}

void TScalarVP::momD1(ScalarField &s, EmField &A, FFT &fft)
{
    ScalarField buf(env().getGrid()), result(env().getGrid()),
                Amu(env().getGrid());
    Complex     ci(0.0,1.0);

    result = zero;

    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        Amu = peekLorentz(A, mu);
        buf = (*phase_[mu])*s;
        fft.FFT_all_dim(buf, buf, FFT::backward);
        buf = Amu*buf;
        fft.FFT_all_dim(buf, buf, FFT::forward);
        result = result - ci*buf;
    }
    fft.FFT_all_dim(s, s, FFT::backward);
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        Amu = peekLorentz(A, mu);
        buf = Amu*s;
        fft.FFT_all_dim(buf, buf, FFT::forward);
        result = result + ci*adj(*phase_[mu])*buf;
    }

    s = result;
}

void TScalarVP::momD2(ScalarField &s, EmField &Asquared, FFT &fft)
{
    ScalarField buf(env().getGrid()), result(env().getGrid()),
                A2mu(env().getGrid());

    result = zero;
    
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        A2mu = peekLorentz(Asquared, mu);
        buf = (*phase_[mu])*s;
        fft.FFT_all_dim(buf, buf, FFT::backward);
        buf = A2mu*buf;
        fft.FFT_all_dim(buf, buf, FFT::forward);
        result = result + .5*buf;
    }
    fft.FFT_all_dim(s, s, FFT::backward);
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        A2mu = peekLorentz(Asquared, mu);        
        buf = A2mu*s;
        fft.FFT_all_dim(buf, buf, FFT::forward);
        result = result + .5*adj(*phase_[mu])*buf;
    }

    s = result;
}
