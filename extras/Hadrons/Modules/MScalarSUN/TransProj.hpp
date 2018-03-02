#ifndef Hadrons_MScalarSUN_TransProj_hpp_
#define Hadrons_MScalarSUN_TransProj_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/Modules/MScalarSUN/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Transverse projection                              *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class TransProjPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TransProjPar,
                                    std::string,  op,
                                    DiffType,     type,
                                    std::string,  output);
};

template <typename SImpl>
class TTransProj: public Module<TransProjPar>
{
public:
    typedef typename SImpl::Field        Field;
    typedef typename SImpl::ComplexField ComplexField;
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::string, op,
                                        Complex    , value);
    };
public:
    // constructor
    TTransProj(const std::string name);
    // destructor
    virtual ~TTransProj(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution

    virtual void execute(void);
};

MODULE_REGISTER_NS(TransProjSU2, TTransProj<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_NS(TransProjSU3, TTransProj<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_NS(TransProjSU4, TTransProj<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_NS(TransProjSU5, TTransProj<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_NS(TransProjSU6, TTransProj<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                        TTransProj implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TTransProj<SImpl>::TTransProj(const std::string name)
: Module<TransProjPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TTransProj<SImpl>::getInput(void)
{
    std::vector<std::string> in = {par().op};
    
    return in;
}

template <typename SImpl>
std::vector<std::string> TTransProj<SImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    for (unsigned int nu = mu; nu < env().getNd(); ++nu)
    {
        out.push_back(varName(getName(), mu, nu));
    }
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTransProj<SImpl>::setup(void)
{
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    for (unsigned int nu = mu; nu < env().getNd(); ++nu)
    {
        envCreateLat(ComplexField, varName(getName(), mu, nu));
    }
    envTmpLat(ComplexField, "buf1");
    envTmpLat(ComplexField, "buf2");
    envTmpLat(ComplexField, "lap");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTransProj<SImpl>::execute(void)
{
    LOG(Message) << "Computing (delta_mu,nu d^2 - d_mu*d_nu)*op using " 
                 << par().type << " derivatives and op= '" << par().op 
                 << "'" << std::endl; 

    std::vector<Result> result;
    auto                &op = envGet(ComplexField, par().op);

    envGetTmp(ComplexField, buf1);
    envGetTmp(ComplexField, buf2);
    envGetTmp(ComplexField, lap);
    lap = zero;
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        dmu(buf1, op, mu, par().type);
        dmu(buf2, buf1, mu, par().type);
        lap += buf2;
    }
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    for (unsigned int nu = mu; nu < env().getNd(); ++nu)
    {
        auto &out = envGet(ComplexField, varName(getName(), mu, nu));
        dmu(buf1, op, mu, par().type);
        dmu(buf2, buf1, nu, par().type);
        out = -buf2;
        if (mu == nu)
        {
            out += lap;
        }
        if (!par().output.empty())
        {
            Result r;

            r.op    = "(delta_" + std::to_string(mu) + "," + std::to_string(nu)
                      + " d^2 - d_" + std::to_string(mu) + "*d_" 
                      + std::to_string(nu) + ")*op";
            r.value = TensorRemove(sum(out));
            result.push_back(r);
        }
    }
    if (result.size() > 0)
    {
        saveResult(par().output, "transproj", result);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_TransProj_hpp_
