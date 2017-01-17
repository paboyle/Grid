#ifndef Hadrons_NoiseLoop_hpp_
#define Hadrons_NoiseLoop_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         NoiseLoop                                          *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MLoop)

class NoiseLoopPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(NoiseLoopPar,
                                    std::string,  q,
                                    std::string, eta);
};

template <typename FImpl>
class TNoiseLoop: public Module<NoiseLoopPar>
{
public:
    TYPE_ALIASES(FImpl,);
public:
    // constructor
    TNoiseLoop(const std::string name);
    // destructor
    virtual ~TNoiseLoop(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(NoiseLoop, TNoiseLoop<FIMPL>, MLoop);

/******************************************************************************
 *                 TNoiseLoop implementation                                  *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TNoiseLoop<FImpl>::TNoiseLoop(const std::string name)
: Module<NoiseLoopPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TNoiseLoop<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().eta};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TNoiseLoop<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TNoiseLoop<FImpl>::setup(void)
{
    env().template registerLattice<PropagatorField>(getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TNoiseLoop<FImpl>::execute(void)
{
    PropagatorField &loop = *env().template createLattice<PropagatorField>(getName());
    PropagatorField &q    = *env().template getObject<PropagatorField>(par().q);
    PropagatorField &eta  = *env().template getObject<PropagatorField>(par().eta);
    loop = q*adj(eta);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_NoiseLoop_hpp_
