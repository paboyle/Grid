#ifndef Hadrons_ScalarVP_hpp_
#define Hadrons_ScalarVP_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         ScalarVP                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalar)

class ScalarVPPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ScalarVPPar,
                                    std::string, emField,
                                    std::string, source,
                                    double,      mass,
                                    double,      charge,
                                    std::string, output);
};

class TScalarVP: public Module<ScalarVPPar>
{
public:
    SCALAR_TYPE_ALIASES(SIMPL,);
    typedef PhotonR::GaugeField     EmField;
    typedef PhotonR::GaugeLinkField EmComp;
public:
    // constructor
    TScalarVP(const std::string name);
    // destructor
    virtual ~TScalarVP(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void chargedProp(ScalarField &prop_q, ScalarField &prop_sun,
                            ScalarField &prop_tad, ScalarField &GFSrc,
                            FFT &fft);
    void momD1(ScalarField &s, FFT &fft);
    void momD2(ScalarField &s, FFT &fft);
private:
    std::string                                 freeMomPropName_, GFSrcName_,
                                                prop0Name_, propQName_,
                                                propSunName_, propTadName_;
    std::vector<std::string>                    phaseName_, muPropQName_;
    std::vector<std::vector<std::string> >      vpTensorName_;
    ScalarField                                 *freeMomProp_, *GFSrc_,
                                                *prop0_;
    std::vector<ScalarField *>                  phase_;
    EmField                                     *A;
};

MODULE_REGISTER_NS(ScalarVP, TScalarVP, MScalar);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_ScalarVP_hpp_
