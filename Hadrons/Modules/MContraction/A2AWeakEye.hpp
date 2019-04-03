#ifndef Hadrons_MContraction_A2AWeakEye_hpp_
#define Hadrons_MContraction_A2AWeakEye_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         A2AWeakEye                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2AWeakEyePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AWeakEyePar,
                                    std::string,  propT0,
                                    std::string,  propLoop,
                                    unsigned int, dtmin,
                                    unsigned int, dtmax,
                                    std::string,  output);
};

template <typename FImpl>
class TA2AWeakEye: public Module<A2AWeakEyePar>
{
  public:
    FERM_TYPE_ALIASES(FImpl, );
    class Metadata : Serializable
    {
      public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        Gamma::Algebra, op,
                                        unsigned int, trace,
                                        unsigned int, dtmin,
                                        unsigned int, dtmax);
    };
    typedef Correlator<Metadata> Result;

  public:
    // constructor
    TA2AWeakEye(const std::string name);
    // destructor
    virtual ~TA2AWeakEye(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(A2AWeakEye, TA2AWeakEye<FIMPL>, MContraction);

/******************************************************************************
 *                 TA2AWeakEye implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AWeakEye<FImpl>::TA2AWeakEye(const std::string name)
: Module<A2AWeakEyePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AWeakEye<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().propT0, par().propLoop};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TA2AWeakEye<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AWeakEye<FImpl>::setup(void)
{
    int nt = env().getDim(Tp);

    envTmp(std::vector<TComplex>, "corrTmp",    1, nt);
    envTmp(std::vector<Complex>,  "corrSaucer", 1, nt);
    envTmp(std::vector<Complex>,  "corrEye",    1, nt);

    envTmp(std::vector<PropagatorField>, "propT0G5", 1, nt, envGetGrid(PropagatorField));

    envTmp(ComplexField, "saucerField", 1, envGetGrid(ComplexField));
    envTmp(ComplexField, "eyeField",    1, envGetGrid(ComplexField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AWeakEye<FImpl>::execute(void)
{
    auto &propT0   = envGet(std::vector<PropagatorField>, par().propT0);
    auto &propLoop = envGet(PropagatorField, par().propLoop);
    int dtmin = par().dtmin;
    int dtmax = par().dtmax;

    GridBase *grid = propT0[0]._grid;
    int nt = env().getDim(Tp);

    // Implicit gamma-5
    auto G5 = Gamma(Gamma::Algebra::Gamma5);
    envGetTmp(std::vector<PropagatorField>, propT0G5);
    for (int t = 0; t < nt; t++)
    {
        propT0G5[t] = propT0[t] * G5;
    }

    std::vector<Result> result;
    Result res;

    envGetTmp(std::vector<TComplex>, corrTmp);
    envGetTmp(std::vector<Complex>,  corrSaucer);
    envGetTmp(std::vector<Complex>,  corrEye);

    envGetTmp(ComplexField, saucerField);
    envGetTmp(ComplexField, eyeField);

    LOG(Message) << "Computing A2A weak eye diagrams using propagators: "
                 << par().propT0 << " and " << par().propLoop << "." << std::endl;
    LOG(Message) << " dt " << dtmin << "..." << dtmax << std::endl;

    for (auto &G : Gamma::gall)
    {
        std::vector<Gamma> GG({G});
        res.info.op = G.g;
        res.info.dtmin = dtmin;
        res.info.dtmax = dtmax;

        for (int t = 0; t < nt; t++)
        {
            corrSaucer[t] = 0.0;
            corrEye[t] = 0.0;
        }
        for (int t0 = 0; t0 < nt; t0++)
        {
            A2Autils<FImpl>::ContractFourQuarkColourDiagonal(propT0G5[t0], propLoop, GG, GG, eyeField, saucerField);

            sliceSum(saucerField, corrTmp, Tp);

            for (int t = 0; t < nt; t++)
            {
                corrSaucer[t] += 2.0 * corrTmp[(t + t0) % nt]()()();
            }

            sliceSum(eyeField, corrTmp, Tp);

            for (int t = 0; t < nt; t++)
            {
                corrEye[t] += 2.0 * corrTmp[(t + t0) % nt]()()();
            }
        }

        res.corr.clear();
        for (int t = 0; t < nt; t++)
        {
            res.corr.push_back(corrSaucer[t]);
        }
        res.info.trace = 1;
        result.push_back(res);
        res.corr.clear();
        for (int t = 0; t < nt; t++)
        {
            res.corr.push_back(corrEye[t]);
        }
        res.info.trace = 2;
        result.push_back(res);
    }
    LOG(Message) << "Writing results to " << par().output << "." << std::endl;
    saveResult(par().output, "A2AWeakEye", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AWeakEye_hpp_
