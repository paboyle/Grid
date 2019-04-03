#ifndef Hadrons_MContraction_A2AWeakNonEye_hpp_
#define Hadrons_MContraction_A2AWeakNonEye_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         A2AWeakNonEye                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2AWeakNonEyePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AWeakNonEyePar,
                                    std::string,  propT0,
                                    std::string,  propT1,
                                    unsigned int, dtmin,
                                    unsigned int, dtmax,
                                    std::string,  output);
};

template <typename FImpl>
class TA2AWeakNonEye: public Module<A2AWeakNonEyePar>
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
    TA2AWeakNonEye(const std::string name);
    // destructor
    virtual ~TA2AWeakNonEye(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(A2AWeakNonEye, TA2AWeakNonEye<FIMPL>, MContraction);

/******************************************************************************
 *                 TA2AWeakNonEye implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AWeakNonEye<FImpl>::TA2AWeakNonEye(const std::string name)
: Module<A2AWeakNonEyePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AWeakNonEye<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().propT0, par().propT1};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TA2AWeakNonEye<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AWeakNonEye<FImpl>::setup(void)
{
    int nt = env().getDim(Tp);
    
    envTmp(std::vector<TComplex>, "corrTmp",       1, nt);
    envTmp(std::vector<Complex>,  "corrWing",      1, nt);
    envTmp(std::vector<Complex>,  "corrConnected", 1, nt);

    envTmp(std::vector<PropagatorField>, "propT0G5", 1, nt, envGetGrid(PropagatorField));
    envTmp(std::vector<PropagatorField>, "propT1G5", 1, nt, envGetGrid(PropagatorField));

    envTmp(ComplexField, "wingField",      1, envGetGrid(ComplexField));
    envTmp(ComplexField, "connectedField", 1, envGetGrid(ComplexField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AWeakNonEye<FImpl>::execute(void)
{
    auto &propT0 = envGet(std::vector<PropagatorField>, par().propT0);
    auto &propT1 = envGet(std::vector<PropagatorField>, par().propT1);
    int dtmin = par().dtmin;
    int dtmax = par().dtmax;
    int nt = env().getDim(Tp);

    Gamma G5 = Gamma(Gamma::Algebra::Gamma5);
    envGetTmp(std::vector<PropagatorField>, propT0G5);
    envGetTmp(std::vector<PropagatorField>, propT1G5);
    for (int t = 0; t < nt; t++)
    {
        propT0G5[t] = propT0[t] * G5;
        propT1G5[t] = propT1[t] * G5;
    }

    std::vector<Result> result;
    Result              res;

    envGetTmp(std::vector<TComplex>, corrTmp);
    envGetTmp(std::vector<Complex>,  corrWing);
    envGetTmp(std::vector<Complex>,  corrConnected);

    envGetTmp(ComplexField, connectedField);
    envGetTmp(ComplexField, wingField);

    LOG(Message) << "Computing A2A weak non-eye diagrams using propagators: "
                 << par().propT0 << " and " << par().propT1 << "." << std::endl;
    LOG(Message) << " dt " << dtmin << "..." << dtmax << std::endl;
    for (auto &G : Gamma::gall)
    {
        std::vector<Gamma> GG({G});
        res.info.op    = G.g;
        res.info.dtmin = dtmin;
        res.info.dtmax = dtmax;
        for(int t = 0; t < nt; t++)
        {
            corrConnected[t] = 0.0;
            corrWing[t]      = 0.0;
        }
        for (int t0 = 0; t0 < nt; t0++)
        {
            for (int dt = dtmin; dt < dtmax; dt++)
            {
                int t1 = (t0 + dt) % nt;
                A2Autils<FImpl>::ContractFourQuarkColourDiagonal(propT0G5[t0], propT1G5[t1], GG, GG, wingField, connectedField);

                sliceSum(connectedField, corrTmp, Tp);
                for (int t = 0; t < nt; t++)
                {
                    corrConnected[t] += 2.0 * corrTmp[(t + t0) % nt]()()();
                }

                sliceSum(wingField, corrTmp, Tp);
                for (int t = 0; t < nt; t++)
                {
                    corrWing[t] += 2.0 * corrTmp[(t + t0) % nt]()()();
                }
            }
        }
        res.corr.clear();
        for (int t = 0; t < nt; t++)
        {
            res.corr.push_back(corrConnected[t]);
        }
        res.info.trace = 1;
        result.push_back(res);
        res.corr.clear();
        for (int t = 0; t < nt; t++)
        {
            res.corr.push_back(corrWing[t]);
        }
        res.info.trace = 2;
        result.push_back(res);
    }

    LOG(Message) << "Writing results to " << par().output << "." << std::endl;
    saveResult(par().output, "A2AWeakNonEye", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AWeakNonEye_hpp_
