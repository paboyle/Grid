#ifndef Hadrons_MSink_Smear_hpp_
#define Hadrons_MSink_Smear_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                 Smear                                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSink)

class SmearPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SmearPar,
                                    std::string, q,
                                    std::string, sink);
};

template <typename FImpl>
class TSmear: public Module<SmearPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SINK_TYPE_ALIASES();
public:
    // constructor
    TSmear(const std::string name);
    // destructor
    virtual ~TSmear(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(Smear, TSmear<FIMPL>, MSink);

/******************************************************************************
 *                          TSmear implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSmear<FImpl>::TSmear(const std::string name)
: Module<SmearPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSmear<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().sink};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSmear<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSmear<FImpl>::setup(void)
{
    unsigned int nt = env().getDim(Tp);
    unsigned int size = nt * sizeof(SitePropagator);
    env().registerObject(getName(), size);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSmear<FImpl>::execute(void)
{
    LOG(Message) << "Sink smearing propagator '" << par().q
                 << "' using sink function '" << par().sink << "'."
                 << std::endl;

    SinkFn          &sink = *env().template getObject<SinkFn>(par().sink);
    PropagatorField &q    = *env().template getObject<PropagatorField>(par().q);
    SlicedPropagator *out = new SlicedPropagator(env().getDim(Tp));
    *out  = sink(q);
    env().setObject(getName(), out);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSink_Smear_hpp_
