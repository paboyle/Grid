#ifndef Hadrons_Wall_hpp_
#define Hadrons_Wall_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Wall                                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class WallPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WallPar,
                                    unsigned int, tW);
};

template <typename FImpl>
class TWall: public Module<WallPar>
{
public:
    TYPE_ALIASES(FImpl,);
public:
    // constructor
    TWall(const std::string name);
    // destructor
    virtual ~TWall(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(Wall, TWall<FIMPL>, MSource);

/******************************************************************************
 *                 TWall implementation                                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TWall<FImpl>::TWall(const std::string name)
: Module<WallPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TWall<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TWall<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWall<FImpl>::setup(void)
{
    env().template registerLattice<PropagatorField>(getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWall<FImpl>::execute(void)
{
    Lattice<iScalar<vInteger>> t(env().getGrid());
    
    LOG(Message) << "Generating wall source at t = " << par().tW << std::endl;
    PropagatorField &src = *env().template createLattice<PropagatorField>(getName());
    LatticeCoordinate(t, Tp);
    src = 1.;
    src = where((t == par().tW), src, 0.*src);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_Wall_hpp_
