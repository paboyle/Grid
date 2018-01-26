#ifndef Hadrons_MScalar_VPCounterTerms_hpp_
#define Hadrons_MScalar_VPCounterTerms_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         VPCounterTerms                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalar)

class VPCounterTermsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(VPCounterTermsPar,
                                    std::string, source,
                                    double,      mass,
                                    std::string, output,
                                    std::vector<std::string>, outputMom);
};

class TVPCounterTerms: public Module<VPCounterTermsPar>
{
public:
    SCALAR_TYPE_ALIASES(SIMPL,);
public:
    // constructor
    TVPCounterTerms(const std::string name);
    // destructor
    virtual ~TVPCounterTerms(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void writeVP(const ScalarField &vp, std::string dsetName);
private:
    std::string                freeMomPropName_, GFSrcName_, phatsqName_, prop0Name_,
                               twoscalarName_, twoscalarVertexName_,
                               psquaredName_, psquaredVertexName_;
    std::vector<std::string>   phaseName_, momPhaseName_;
    std::vector<ScalarField *> phase_, momPhase_;
};

MODULE_REGISTER_NS(VPCounterTerms, TVPCounterTerms, MScalar);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalar_VPCounterTerms_hpp_
