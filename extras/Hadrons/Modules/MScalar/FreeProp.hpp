#ifndef Hadrons_FreeProp_hpp_
#define Hadrons_FreeProp_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         FreeProp                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalar)

class FreePropPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FreePropPar,
                                    std::string, source,
                                    double,      mass);
};

class TFreeProp: public Module<FreePropPar>
{
public:
    SCALAR_TYPE_ALIASES(SIMPL,);
public:
    // constructor
    TFreeProp(const std::string name);
    // destructor
    virtual ~TFreeProp(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(FreeProp, TFreeProp, MScalar);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_FreeProp_hpp_
