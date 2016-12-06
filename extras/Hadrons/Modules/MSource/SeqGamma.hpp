#ifndef Hadrons_SeqGamma_hpp_
#define Hadrons_SeqGamma_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Sequential source
 -----------------------------
 * src_x = q_x * theta(x_3 - tA) * theta(tB - x_3) * gamma * exp(i x.mom)
 
 * options:
 - q: input propagator (string)
 - tA: begin timeslice (integer)
 - tB: end timesilce (integer)
 - gamma: gamma product to insert (integer)
 - mom: momentum insertion, space-separated float sequence (e.g ".1 .2 1. 0.")
 
 */

/******************************************************************************
 *                         SeqGamma                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class SeqGammaPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SeqGammaPar,
                                    std::string,  q,
                                    unsigned int, tA,
                                    unsigned int, tB,
                                    unsigned int, gamma,
                                    std::string,  mom);
};

template <typename FImpl>
class TSeqGamma: public Module<SeqGammaPar>
{
public:
    TYPE_ALIASES(FImpl,);
public:
    // constructor
    TSeqGamma(const std::string name);
    // destructor
    virtual ~TSeqGamma(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

typedef TSeqGamma<FIMPL> SeqGamma;

/******************************************************************************
 *                 TSeqGamma implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSeqGamma<FImpl>::TSeqGamma(const std::string name)
: Module<SeqGammaPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSeqGamma<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSeqGamma<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqGamma<FImpl>::setup(void)
{
    env().template registerLattice<PropagatorField>(getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqGamma<FImpl>::execute(void)
{
    if (par().tA == par().tB)
    {
        LOG(Message) << "Generating gamma_" << par().gamma
                     << " sequential source at t= " << par().tA << std::endl;
    }
    else
    {
        LOG(Message) << "Generating gamma_" << par().gamma
                     << " sequential source for "
                     << par().tA << " <= t <= " << par().tB << std::endl;
    }
    PropagatorField &src = *env().template createLattice<PropagatorField>(getName());
    PropagatorField &q   = *env().template getObject<PropagatorField>(par().q);
    Lattice<iScalar<vInteger>> t(env().getGrid());
    LatticeComplex             ph(env().getGrid()), coor(env().getGrid());
    SpinMatrix                 g;
    std::vector<Real>          p;
    Complex                    i(0.0,1.0);
    
    g  = makeGammaProd(par().gamma);
    p  = strToVec<Real>(par().mom);
    ph = zero;
    for(unsigned int mu = 0; mu < Nd; mu++)
    {
        LatticeCoordinate(coor, mu);
        ph = ph + p[mu]*coor;
    }
    ph = exp(i*ph);
    LatticeCoordinate(t, Tp);
    src = where((t >= par().tA) and (t <= par().tB), g*ph*q, 0.*q);
}

END_MODULE_NAMESPACE

MODULE_REGISTER_NS(SeqGamma, MSource);

END_HADRONS_NAMESPACE

#endif // Hadrons_SeqGamma_hpp_
