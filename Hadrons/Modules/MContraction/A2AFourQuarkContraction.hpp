#ifndef Hadrons_MContraction_A2AFourQuarkContraction_hpp_
#define Hadrons_MContraction_A2AFourQuarkContraction_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DiskVector.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         A2AFourQuarkContraction                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2AFourQuarkContractionPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AFourQuarkContractionPar,
                                    std::string, v1,
                                    std::string, v2,
                                    std::string, mf12);
};

template <typename FImpl>
class TA2AFourQuarkContraction: public Module<A2AFourQuarkContractionPar>
{
  public:
    FERM_TYPE_ALIASES(FImpl, );
    // constructor
    TA2AFourQuarkContraction(const std::string name);
    // destructor
    virtual ~TA2AFourQuarkContraction(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(A2AFourQuarkContraction, TA2AFourQuarkContraction<FIMPL>, MContraction);

/******************************************************************************
 *                 TA2AFourQuarkContraction implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AFourQuarkContraction<FImpl>::TA2AFourQuarkContraction(const std::string name)
: Module<A2AFourQuarkContractionPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AFourQuarkContraction<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().v1, par().v2, par().mf12};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TA2AFourQuarkContraction<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AFourQuarkContraction<FImpl>::setup(void)
{
    int nt = env().getDim(Tp);

    envCreate(std::vector<PropagatorField>, getName(), 1, nt, envGetGrid(PropagatorField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AFourQuarkContraction<FImpl>::execute(void)
{
    int nt = env().getDim(Tp);

    auto &v1   = envGet(std::vector<FermionField>, par().v1);
    auto &v2   = envGet(std::vector<FermionField>, par().v2);
    auto &mf12 = envGet(EigenDiskVector<Complex>, par().mf12);

    auto &wwvv = envGet(std::vector<PropagatorField>, getName());

    for (int t = 0; t < nt; t++)
    {
        wwvv[t] = zero;
    }

    LOG(Message) << "Computing 4 quark contraction for: " << getName() << std::endl;
    A2Autils<FImpl>::ContractWWVV(wwvv, mf12, &v1[0], &v2[0]);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AFourQuarkContraction_hpp_
