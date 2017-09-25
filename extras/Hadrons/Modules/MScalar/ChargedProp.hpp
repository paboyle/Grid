#ifndef Hadrons_MScalar_ChargedProp_hpp_
#define Hadrons_MScalar_ChargedProp_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Charged scalar propagator                            *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalar)

class ChargedPropPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ChargedPropPar,
                                    std::string, emField,
                                    std::string, source,
                                    double,      mass,
                                    double,      charge,
                                    std::string, output,
                                    std::vector<std::string>, outputMom);
};

class TChargedProp: public Module<ChargedPropPar>
{
public:
    SCALAR_TYPE_ALIASES(SIMPL,);
    typedef PhotonR::GaugeField     EmField;
    typedef PhotonR::GaugeLinkField EmComp;
public:
    // constructor
    TChargedProp(const std::string name);
    // destructor
    virtual ~TChargedProp(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void momD1(ScalarField &s, FFT &fft);
    void momD2(ScalarField &s, FFT &fft);
private:
    std::string                freeMomPropName_, GFSrcName_, prop0Name_,
                               propQName_, propSunName_, propTadName_;
    std::vector<std::string>   phaseName_;
    ScalarField                *freeMomProp_, *GFSrc_, *prop0_;
    std::vector<ScalarField *> phase_;
    EmField                    *A;
};

MODULE_REGISTER_NS(ChargedProp, TChargedProp, MScalar);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalar_ChargedProp_hpp_
