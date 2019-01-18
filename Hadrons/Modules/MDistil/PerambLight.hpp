#ifndef Hadrons_MDistil_perambulator_l_hpp_
#define Hadrons_MDistil_perambulator_l_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         perambulator_l                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

class perambulator_lPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(perambulator_lPar,
                                    unsigned int, i);
};

template <typename FImpl>
class Tperambulator_l: public Module<perambulator_lPar>
{
public:
    // constructor
    Tperambulator_l(const std::string name);
    // destructor
    virtual ~Tperambulator_l(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(perambulator_l, Tperambulator_l<FIMPL>, MDistil);

/******************************************************************************
 *                 Tperambulator_l implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
Tperambulator_l<FImpl>::Tperambulator_l(const std::string name)
: Module<perambulator_lPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> Tperambulator_l<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> Tperambulator_l<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void Tperambulator_l<FImpl>::setup(void)
{
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void Tperambulator_l<FImpl>::execute(void)
{
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_perambulator_l_hpp_
