#ifndef Hadrons_WardIdentity_hpp_
#define Hadrons_WardIdentity_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
  Ward Identity contractions
 -----------------------------
 
 * options:
 - q:      propagator, 5D if available (string)
 - q4d:    4D propagator, duplicate of q if q is not 5D (string)
 - action: action module used for propagator solution (string)
 - mass:   mass of quark (double)
*/

/******************************************************************************
 *                              WardIdentity                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class WardIdentityPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WardIdentityPar,
                                    std::string, q,
                                    std::string, q4d,
                                    std::string, action,
                                    double,      mass);
};

template <typename FImpl>
class TWardIdentity: public Module<WardIdentityPar>
{
public:
    TYPE_ALIASES(FImpl,);
public:
    // constructor
    TWardIdentity(const std::string name);
    // destructor
    virtual ~TWardIdentity(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Ls_;
};

MODULE_REGISTER_NS(WardIdentity, TWardIdentity<FIMPL>, MContraction);

/******************************************************************************
 *                     TWardIdentity implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TWardIdentity<FImpl>::TWardIdentity(const std::string name)
: Module<WardIdentityPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TWardIdentity<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().q4d, par().action};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TWardIdentity<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWardIdentity<FImpl>::setup(void)
{
    Ls_ = env().getObjectLs(par().q);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWardIdentity<FImpl>::execute(void)
{
    LOG(Message) << "Performing Ward Identity checks for quark '" << par().q
                 << "'." << std::endl;

    PropagatorField psi(env().getGrid()), tmp(env().getGrid());
    PropagatorField q    = *env().template getObject<PropagatorField>(par().q);
    PropagatorField q4d  = *env().template getObject<PropagatorField>(par().q4d);
    FMat            &act = *(env().template getObject<FMat>(par().action));
    Gamma           g5(Gamma::Algebra::Gamma5);
    LatticeComplex  PP(env().getGrid()), PA(env().getGrid()),
                    c(env().getGrid()), PJ5q(env().getGrid()),
                    vector_WI(env().getGrid()), defect(env().getGrid());
    c = zero; PJ5q = zero; vector_WI = zero; defect = zero;
    std::vector<LatticeComplex> Vmu(Nd, c);
    std::vector<LatticeComplex> Amu(Nd, c);

    // Get PP, PA, V_mu, A_mu for 4D.
    PP = trace(adj(q4d)*q4d);
    PA = trace(adj(q4d)*g5*q4d);
    for (unsigned int mu = 0; mu < Nd; ++mu)
    {
        act.ContractConservedCurrent(q, q, tmp, Current::Vector, mu);
        Vmu[mu] = trace(g5*tmp);
        act.ContractConservedCurrent(q, q, tmp, Current::Axial, mu);
        Amu[mu] = trace(g5*tmp);
    }

    // Get PJ5q for 5D (zero for 4D).
    if (Ls_ > 1)
    {
        ExtractSlice(psi, q, 0, Ls_/2 - 1);
        psi  = 0.5 * (psi + g5*psi);
        ExtractSlice(tmp, q, 0, Ls_/2);
        psi += 0.5 * (tmp - g5*tmp);
        PJ5q = trace(adj(psi)*psi);
    }

    // Test ward identities, D_mu V_mu = 0; D_mu A_mu = 2m<PP> + 2 PJ5q
    for (unsigned int mu = 0; mu < Nd; ++mu)
    {
        vector_WI += Vmu[mu] - Cshift(Vmu[mu], mu, -1); 
        defect    += Amu[mu] - Cshift(Amu[mu], mu, -1);
    }
    defect -= 2.*PJ5q;
    defect -= 2.*(par().mass)*PP;

    LOG(Message) << "Vector Ward Identity check Delta_mu V_mu = " 
                 << norm2(vector_WI) << std::endl;
    LOG(Message) << "Axial Ward Identity defect Delta_mu A_mu = "
                 << norm2(defect) << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_WardIdentity_hpp_
