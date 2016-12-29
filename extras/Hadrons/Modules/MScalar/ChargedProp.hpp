#ifndef Hadrons_ChargedProp_hpp_
#define Hadrons_ChargedProp_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         ChargedProp                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalar)

class ChargedPropPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ChargedPropPar,
                                    unsigned int, i);
};

class TChargedProp: public Module<ChargedPropPar>
{
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
};

MODULE_REGISTER_NS(ChargedProp, TChargedProp, MScalar);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_ChargedProp_hpp_
