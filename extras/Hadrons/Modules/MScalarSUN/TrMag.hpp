#ifndef Hadrons_MScalarSUN_TrMag_hpp_
#define Hadrons_MScalarSUN_TrMag_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         TrMag                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class TrMagPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TrMagPar,
                                    std::string,  field,
                                    unsigned int, maxPow,
                                    std::string,  output);
};

template <typename SImpl>
class TTrMag: public Module<TrMagPar>
{
public:
    typedef typename SImpl::Field        Field;
    typedef typename SImpl::ComplexField ComplexField;
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::string, op,
                                        Real,        value);
    };
public:
    // constructor
    TTrMag(const std::string name);
    // destructor
    virtual ~TTrMag(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(TrMagSU2, TTrMag<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_NS(TrMagSU3, TTrMag<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_NS(TrMagSU4, TTrMag<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_NS(TrMagSU5, TTrMag<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_NS(TrMagSU6, TTrMag<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                 TTrMag implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TTrMag<SImpl>::TTrMag(const std::string name)
: Module<TrMagPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TTrMag<SImpl>::getInput(void)
{
    std::vector<std::string> in = {par().field};
    
    return in;
}

template <typename SImpl>
std::vector<std::string> TTrMag<SImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTrMag<SImpl>::setup(void)
{}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTrMag<SImpl>::execute(void)
{
    LOG(Message) << "Computing tr(mag^n) for n even up to " << par().maxPow
                 << "..." << std::endl;

    std::vector<Result> result;
    ResultWriter        writer(RESULT_FILE_NAME(par().output));
    auto                &phi = envGet(Field, par().field);

    auto m2 = sum(phi), mn = m2;

    m2 = -m2*m2;
    mn = 1.;
    for (unsigned int n = 2; n <= par().maxPow; n += 2)
    {
        Result r;

        mn = mn*m2;
        r.op    = "tr(mag^" + std::to_string(n) + ")";
        r.value = TensorRemove(trace(mn)).real();
        result.push_back(r);
    }
    write(writer, "trmag", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_TrMag_hpp_
