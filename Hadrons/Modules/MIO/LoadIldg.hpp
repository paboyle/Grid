#ifndef Hadrons_MIO_LoadIldg_hpp_
#define Hadrons_MIO_LoadIldg_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Load an ILDG configuration                           *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadIldgPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadIldgPar,
                                    std::string, file);
};

template <typename GImpl>
class TLoadIldg: public Module<LoadIldgPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TLoadIldg(const std::string name);
    // destructor
    virtual ~TLoadIldg(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadIldg,  TLoadIldg<GIMPL>,  MIO);

/******************************************************************************
*                       TLoadIldg implementation                              * 
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TLoadIldg<GImpl>::TLoadIldg(const std::string name)
: Module<LoadIldgPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TLoadIldg<GImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename GImpl>
std::vector<std::string> TLoadIldg<GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TLoadIldg<GImpl>::setup(void)
{
    envCreateLat(GaugeField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TLoadIldg<GImpl>::execute(void)
{
    FieldMetaData header;
    std::string   fileName = par().file + "."
                             + std::to_string(vm().getTrajectory());
    LOG(Message) << "Loading ILDG configuration from file '" << fileName
                 << "'" << std::endl;
    // Is this right? Or should it be U, with a peek to get Umu?
    auto &Umu = envGet(GaugeField, getName());
    
    IldgReader IR;
    IR.open(fileName);
    IR.readConfiguration(Umu,header);
    IR.close();
}



END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadIldg_hpp_
