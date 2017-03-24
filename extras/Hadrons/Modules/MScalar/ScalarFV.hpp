#ifndef Hadrons_ScalarFV_hpp_
#define Hadrons_ScalarFV_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         ScalarFV                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalar)

class ScalarFVPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ScalarFVPar,
                                    std::string, emField,
                                    std::string, source,
                                    double,      mass,
                                    double,      charge,
                                    std::string, output,
                                    unsigned int, i);
};

class TScalarFV: public Module<ScalarFVPar>
{
public:
    // constructor
    TScalarFV(const std::string name);
    // destructor
    virtual ~TScalarFV(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(ScalarFV, TScalarFV, MScalar);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_ScalarFV_hpp_
