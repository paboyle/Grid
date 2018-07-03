#ifndef Hadrons_MGauge_StoutSmearing_hpp_
#define Hadrons_MGauge_StoutSmearing_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                            Stout smearing                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class StoutSmearingPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StoutSmearingPar,
                                    std::string, gauge,
                                    unsigned int, steps,
                                    double, rho);
};

template <typename GImpl>
class TStoutSmearing: public Module<StoutSmearingPar>
{
public:
    typedef typename GImpl::Field GaugeField;
public:
    // constructor
    TStoutSmearing(const std::string name);
    // destructor
    virtual ~TStoutSmearing(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(StoutSmearing, TStoutSmearing<GIMPL>, MGauge);

/******************************************************************************
 *                     TStoutSmearing implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TStoutSmearing<GImpl>::TStoutSmearing(const std::string name)
: Module<StoutSmearingPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TStoutSmearing<GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename GImpl>
std::vector<std::string> TStoutSmearing<GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TStoutSmearing<GImpl>::setup(void)
{
    envCreateLat(GaugeField, getName());
    envTmpLat(GaugeField, "buf");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TStoutSmearing<GImpl>::execute(void)
{
    LOG(Message) << "Smearing '" << par().gauge << "' with " << par().steps
                 << " step" << ((par().steps > 1) ? "s" : "") 
                 << " of stout smearing and rho= " << par().rho << std::endl;

    Smear_Stout<GImpl> smearer(par().rho);
    auto               &U    = envGet(GaugeField, par().gauge);
    auto               &Usmr = envGet(GaugeField, getName());

    envGetTmp(GaugeField, buf);
    buf = U;
    for (unsigned int n = 0; n < par().steps; ++n)
    {
        smearer.smear(Usmr, buf);
        buf = Usmr;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_StoutSmearing_hpp_
