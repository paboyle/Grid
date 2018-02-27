#ifndef Hadrons_MAction_ZMobiusDWF_hpp_
#define Hadrons_MAction_ZMobiusDWF_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         ZMobiusDWF                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MAction)

class ZMobiusDWFPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ZMobiusDWFPar,
                                    std::string                      , gauge,
                                    unsigned int                     , Ls,
                                    double                           , mass,
                                    double                           , M5,
                                    double                           , b,
                                    double                           , c,
                                    std::vector<std::complex<double>>, omega,
                                    std::string                      , boundary);
};

template <typename FImpl>
class TZMobiusDWF: public Module<ZMobiusDWFPar>
{
public:
    FGS_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TZMobiusDWF(const std::string name);
    // destructor
    virtual ~TZMobiusDWF(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(ZMobiusDWF, TZMobiusDWF<ZFIMPL>, MAction);

/******************************************************************************
 *                 TZMobiusDWF implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TZMobiusDWF<FImpl>::TZMobiusDWF(const std::string name)
: Module<ZMobiusDWFPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TZMobiusDWF<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TZMobiusDWF<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TZMobiusDWF<FImpl>::setup(void)
{
    LOG(Message) << "Setting up z-Mobius domain wall fermion matrix with m= "
                 << par().mass << ", M5= " << par().M5 << ", Ls= "
                 << ", b= " << par().b << ", c= " << par().c
                 << par().Ls << " using gauge field '" << par().gauge << "'"
                 << std::endl;
    LOG(Message) << "Omegas: " << std::endl;
    for (unsigned int i = 0; i < par().omega.size(); ++i)
    {
        LOG(Message) << "  omega[" << i << "]= " << par().omega[i] << std::endl;
    }
    LOG(Message) << "Fermion boundary conditions: " << par().boundary
                 << std::endl;

    env().createGrid(par().Ls);
    auto &U    = envGet(LatticeGaugeField, par().gauge);
    auto &g4   = *env().getGrid();
    auto &grb4 = *env().getRbGrid();
    auto &g5   = *env().getGrid(par().Ls);
    auto &grb5 = *env().getRbGrid(par().Ls);
    auto omega = par().omega;
    std::vector<Complex> boundary = strToVec<Complex>(par().boundary);
    typename ZMobiusFermion<FImpl>::ImplParams implParams(boundary);
    envCreateDerived(FMat, ZMobiusFermion<FImpl>, getName(), par().Ls, U, g5,
                     grb5, g4, grb4, par().mass, par().M5, omega,
                     par().b, par().c, implParams);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TZMobiusDWF<FImpl>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MAction_ZMobiusDWF_hpp_
