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
                                    std::string, output);
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
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void makeCaches(void);
    void momD1(ScalarField &s, FFT &fft);
    void momD2(ScalarField &s, FFT &fft);
private:
    bool                       freeMomPropDone_, GFSrcDone_, phasesDone_;
    std::string                freeMomPropName_, GFSrcName_, fftName_;
    std::vector<std::string>   phaseName_;
    std::vector<ScalarField *> phase_;
};

MODULE_REGISTER_NS(ChargedProp, TChargedProp, MScalar);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalar_ChargedProp_hpp_
