#ifndef Hadrons_Wall_hpp_
#define Hadrons_Wall_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Wall source
 -----------------------------
 * src_x = theta(x_3 - tA) * theta(tB - x_3) * exp(i x.mom)
 
 * options:
 - tA: begin timeslice (integer)
 - tB: end timeslice (integer)
 - mom: momentum insertion, space-separated float sequence (e.g ".1 .2 1. 0.")
 
 */

/******************************************************************************
 *                         Wall                                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class WallPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WallPar,
                                    unsigned int, tW,
                                    std::string, mom);
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
    LOG(Message) << "Generating wall source at t = " << par().tW 
                 << " with momentum " << par().mom << std::endl;
    
    PropagatorField &src = *env().template createLattice<PropagatorField>(getName());
    Lattice<iScalar<vInteger>> t(env().getGrid());
    LatticeComplex             ph(env().getGrid()), coor(env().getGrid());
    std::vector<Real>          p;
    Complex                    i(0.0,1.0);
    
    p  = strToVec<Real>(par().mom);
    ph = zero;
    for(unsigned int mu = 0; mu < Nd; mu++)
    {
        LatticeCoordinate(coor, mu);
        ph = ph + p[mu]*coor;
    }
    ph = exp(i*ph);
    LatticeCoordinate(t, Tp);
    src = 1.;
    src = where((t == par().tW), src*ph, 0.*src);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_Wall_hpp_
