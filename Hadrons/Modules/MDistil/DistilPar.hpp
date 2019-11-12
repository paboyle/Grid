#ifndef Hadrons_MDistil_DistilPar_hpp_
#define Hadrons_MDistil_DistilPar_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MDistil/DistilCommon.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         DistilPar                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

template <typename FImpl>
class TDistilPar: public Module<DistilParameters>
{
public:
    // constructor
    TDistilPar(const std::string name);
    // destructor
    virtual ~TDistilPar(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(DistilPar, TDistilPar<FIMPL>, MDistil);

/******************************************************************************
 *                 TDistilPar implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDistilPar<FImpl>::TDistilPar(const std::string name)
: Module<DistilParameters>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDistilPar<FImpl>::getInput(void)
{
    return {};
}

template <typename FImpl>
std::vector<std::string> TDistilPar<FImpl>::getOutput(void)
{
    return {getName()};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilPar<FImpl>::setup(void)
{
    envCreate(DistilParameters, getName(), 1, par() );
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilPar<FImpl>::execute(void)
{
    // Nothing to do. setup() created and initialised the output object
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DistilPar_hpp_
