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


class DistilParPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilParPar,
                                    int, nvec,
                                    int, nnoise,
                                    int, tsrc,
                                    int, TI,
                                    int, LI,
                                    int, SI );
};

template <typename FImpl>
class TDistilPar: public Module<DistilParPar>
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
: Module<DistilParPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDistilPar<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TDistilPar<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilPar<FImpl>::setup(void)
{
 //   envCreate(Hadrons::MDistil::DistilParameters, getName(), 1); //DOES NOT WORK
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilPar<FImpl>::execute(void)
{
    Hadrons::MDistil::DistilParameters &out = envGet(Hadrons::MDistil::DistilParameters, getName());
  /*  out.nvec=par().nvec;
    out.nnoise=par().nnoise;
    out.tsrc=par().tsrc;
    out.TI=par().TI;
    out.LI=par().LI;
    out.SI=par().SI; */
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DistilPar_hpp_
