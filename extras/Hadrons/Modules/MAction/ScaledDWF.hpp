#ifndef Hadrons_MAction_ScaledDWF_hpp_
#define Hadrons_MAction_ScaledDWF_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                      Scaled domain wall fermion                            *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MAction)

class ScaledDWFPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ScaledDWFPar,
                                    std::string , gauge,
                                    unsigned int, Ls,
                                    double      , mass,
                                    double      , M5,
                                    double      , scale,
                                    std::string , boundary);
};

template <typename FImpl>
class TScaledDWF: public Module<ScaledDWFPar>
{
public:
    FG_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TScaledDWF(const std::string name);
    // destructor
    virtual ~TScaledDWF(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ScaledDWF, TScaledDWF<FIMPL>, MAction);

/******************************************************************************
 *                      TScaledDWF implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TScaledDWF<FImpl>::TScaledDWF(const std::string name)
: Module<ScaledDWFPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TScaledDWF<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TScaledDWF<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TScaledDWF<FImpl>::setup(void)
{
    LOG(Message) << "Setting up scaled domain wall fermion matrix with m= "
                 << par().mass << ", M5= " << par().M5 << ", Ls= " << par().Ls 
                 << ", scale= " << par().scale
                 << " using gauge field '" << par().gauge << "'"
                 << std::endl;
    LOG(Message) << "Fermion boundary conditions: " << par().boundary
                 << std::endl;

    env().createGrid(par().Ls);
    auto &U    = envGet(LatticeGaugeField, par().gauge);
    auto &g4   = *env().getGrid();
    auto &grb4 = *env().getRbGrid();
    auto &g5   = *env().getGrid(par().Ls);
    auto &grb5 = *env().getRbGrid(par().Ls);
    std::vector<Complex> boundary = strToVec<Complex>(par().boundary);
    typename MobiusFermion<FImpl>::ImplParams implParams(boundary);
    envCreateDerived(FMat, ScaledShamirFermion<FImpl>, getName(), par().Ls, U, g5,
                     grb5, g4, grb4, par().mass, par().M5, par().scale,
                     implParams);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TScaledDWF<FImpl>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MAction_ScaledDWF_hpp_
