#ifndef Hadrons_MUtilities_RandomVectors_hpp_
#define Hadrons_MUtilities_RandomVectors_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *           Module generating random lattices for testing purposes           *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class RandomVectorsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RandomVectorsPar,
                                    unsigned int, size,
                                    unsigned int, Ls);
};

template <typename Field>
class TRandomVectors: public Module<RandomVectorsPar>
{
public:
    // constructor
    TRandomVectors(const std::string name);
    // destructor
    virtual ~TRandomVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(RandomFermions, TRandomVectors<FIMPL::FermionField>, MUtilities);

/******************************************************************************
 *                      TRandomVectors implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TRandomVectors<Field>::TRandomVectors(const std::string name)
: Module<RandomVectorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TRandomVectors<Field>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename Field>
std::vector<std::string> TRandomVectors<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TRandomVectors<Field>::setup(void)
{
    envCreate(std::vector<Field>, getName(), par().Ls, par().size, 
              env().getGrid(par().Ls));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TRandomVectors<Field>::execute(void)
{
    LOG(Message) << "Generating " << par().size << " random vectors" << std::endl;

    auto &vec = envGet(std::vector<Field>, getName());
    
    for (unsigned int i = 0; i < vec.size(); ++i)
    {
        random(*env().get4dRng(), vec[i]);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_RandomVectors_hpp_
