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
                                    std::string, scalarProp,
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
    std::string  prop0Name_, propD1Name_, propD1D1Name_, propD2Name_;
};

MODULE_REGISTER_NS(ScalarVP, TScalarVP, MScalar);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_ScalarVP_hpp_
