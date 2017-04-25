#ifndef Hadrons_SeqConserved_hpp_
#define Hadrons_SeqConserved_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Sequential source
 -----------------------------
 * src_x = q_x * theta(x_3 - tA) * theta(tB - x_3) * J_mu * exp(i x.mom)
 
 * options:
 - q: input propagator (string)
 - action: fermion action used for propagator q (string)
 - tA: begin timeslice (integer)
 - tB: end timesilce (integer)
 - curr_type: type of conserved current to insert (Current)
 - mu: Lorentz index of current to insert (integer)
 - mom: momentum insertion, space-separated float sequence (e.g ".1 .2 1. 0.")
 
 */

/******************************************************************************
 *                              SeqConserved                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class SeqConservedPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SeqConservedPar,
                                    std::string,  q,
                                    std::string,  action,
                                    unsigned int, tA,
                                    unsigned int, tB,
                                    Current,      curr_type,
                                    unsigned int, mu,
                                    std::string,  mom);
};

template <typename FImpl>
class TSeqConserved: public Module<SeqConservedPar>
{
public:
    TYPE_ALIASES(FImpl,);
public:
    // constructor
    TSeqConserved(const std::string name);
    // destructor
    virtual ~TSeqConserved(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(SeqConserved, TSeqConserved<FIMPL>, MSource);

/******************************************************************************
 *                      TSeqConserved implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSeqConserved<FImpl>::TSeqConserved(const std::string name)
: Module<SeqConservedPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSeqConserved<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSeqConserved<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqConserved<FImpl>::setup(void)
{
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqConserved<FImpl>::execute(void)
{
    if (par().tA == par().tB)
    {
        LOG(Message) << "Generating sequential source with conserved "
                     << par().curr_type << " current insertion (mu = " 
                     << par().mu << ") at " << "t = " << par().tA << std::endl;
    }
    else
    {
        LOG(Message) << "Generating sequential source with conserved "
                     << par().curr_type << " current insertion (mu = " 
                     << par().mu << ") for " << par().tA << " <= t <= " 
                     << par().tB << std::endl;
    }
    PropagatorField &src = *env().template createLattice<PropagatorField>(getName());
    PropagatorField &q   = *env().template getObject<PropagatorField>(par().q);
    FMat            &mat = *(env().template getObject<FMat>(par().action));

    std::vector<Real> mom = strToVec<Real>(par().mom);
    mat.SeqConservedCurrent(q, src, par().curr_type, par().mu, 
                            mom, par().tA, par().tB);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_SeqConserved_hpp_
