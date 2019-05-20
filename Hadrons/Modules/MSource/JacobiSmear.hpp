#ifndef Hadrons_MSource_JacobiSmear_hpp_
#define Hadrons_MSource_JacobiSmear_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         JacobiSmear                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class JacobiSmearPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(JacobiSmearPar,
                                    std::string, gauge,
                                    double, width,
                                    int, iterations,
                                    int, orthog,
                                    std::string, source);
};

template <typename FImpl>
class TJacobiSmear: public Module<JacobiSmearPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename FImpl::GaugeLinkField GaugeMat;
public:
    // constructor
    TJacobiSmear(const std::string name);
    // destructor
    virtual ~TJacobiSmear(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(JacobiSmear, TJacobiSmear<FIMPL>, MSource);

/******************************************************************************
 *                 TJacobiSmear implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TJacobiSmear<FImpl>::TJacobiSmear(const std::string name)
: Module<JacobiSmearPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TJacobiSmear<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source, par().gauge};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TJacobiSmear<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TJacobiSmear<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
    envTmp(std::vector<GaugeMat>, "Umu", 1, 4, envGetGrid(LatticeColourMatrix));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TJacobiSmear<FImpl>::execute(void)
{
    auto &out = envGet(PropagatorField, getName());
    auto &src = envGet(PropagatorField, par().source);
    auto &U = envGet(GaugeField, par().gauge);
    envGetTmp(std::vector<GaugeMat>, Umu);
    for(int mu=0; mu<4; mu++)
    {
       Umu.at(mu)=peekLorentz(U,mu);
    }
    CovariantSmearing<FImpl> covsmear;
    out=src;
    startTimer("Jacobi iteration");
    covsmear.GaussianSmear(Umu, out, par().width, par().iterations, par().orthog);
    stopTimer("Jacobi iteration");
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_JacobiSmear_hpp_
