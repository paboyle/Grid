#ifndef Hadrons_MScalarSUN_TrPhi_hpp_
#define Hadrons_MScalarSUN_TrPhi_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         TrPhi                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class TrPhiPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TrPhiPar,
                                    std::string,  field,
                                    unsigned int, maxPow,
                                    std::string,  output);
};

template <typename SImpl>
class TTrPhi: public Module<TrPhiPar>
{
public:
    typedef typename SImpl::Field        Field;
    typedef typename SImpl::ComplexField ComplexField;
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::string, op,
                                        Complex,     value);
    };
public:
    // constructor
    TTrPhi(const std::string name);
    // destructor
    virtual ~TTrPhi(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    // output name generator
    std::string outName(const unsigned int n);
};

MODULE_REGISTER_NS(TrPhiSU2, TTrPhi<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_NS(TrPhiSU3, TTrPhi<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_NS(TrPhiSU4, TTrPhi<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_NS(TrPhiSU5, TTrPhi<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_NS(TrPhiSU6, TTrPhi<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                 TTrPhi implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TTrPhi<SImpl>::TTrPhi(const std::string name)
: Module<TrPhiPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TTrPhi<SImpl>::getInput(void)
{
    std::vector<std::string> in = {par().field};
    
    return in;
}

template <typename SImpl>
std::vector<std::string> TTrPhi<SImpl>::getOutput(void)
{
    std::vector<std::string> out;

    for (unsigned int n = 2; n <= par().maxPow; n += 2)
    {
        out.push_back(outName(n));
    }
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTrPhi<SImpl>::setup(void)
{
    if (par().maxPow < 2)
    {
        HADRON_ERROR(Size, "'maxPow' should be at least equal to 2");
    }
    envTmpLat(Field, "phi2");
    envTmpLat(Field, "buf");
    for (unsigned int n = 2; n <= par().maxPow; n += 2)
    {
        envCreateLat(ComplexField, outName(n));
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTrPhi<SImpl>::execute(void)
{
    LOG(Message) << "Computing tr(phi^n) for n even up to " << par().maxPow
                 << "..." << std::endl; 

    std::vector<Result> result;
    auto                &phi = envGet(Field, par().field);

    envGetTmp(Field, phi2);
    envGetTmp(Field, buf);
    buf  = 1.;
    phi2 = -phi*phi; 
    for (unsigned int n = 2; n <= par().maxPow; n += 2)
    {
        auto &phin = envGet(ComplexField, outName(n));

        buf  = buf*phi2;
        phin = trace(buf);
        if (!par().output.empty())
        {
            Result r;

            r.op    = "phi" + std::to_string(n);
            r.value = TensorRemove(sum(phin));
            result.push_back(r);
        }
    }
    if (result.size() > 0)
    {
        ResultWriter writer(RESULT_FILE_NAME(par().output));

        write(writer, "trphi", result);
    }
}

// output name generator ///////////////////////////////////////////////////////
template <typename SImpl>
std::string TTrPhi<SImpl>::outName(const unsigned int n)
{
    return getName() + "_" + std::to_string(n);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_TrPhi_hpp_
