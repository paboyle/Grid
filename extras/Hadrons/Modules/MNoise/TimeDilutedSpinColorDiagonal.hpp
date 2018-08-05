#ifndef Hadrons_MNoise_TimeDilutedSpinColorDiagonal_hpp_
#define Hadrons_MNoise_TimeDilutedSpinColorDiagonal_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *             Generate time diluted spin-color diagonal noise                *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNoise)

template <typename FImpl>
class TTimeDilutedSpinColorDiagonal: public Module<NoPar>
{
public:
    // constructor
    TTimeDilutedSpinColorDiagonal(const std::string name);
    // destructor
    virtual ~TTimeDilutedSpinColorDiagonal(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(TimeDilutedSpinColorDiagonal, TTimeDilutedSpinColorDiagonal<FIMPL>, MNoise);

/******************************************************************************
 *              TTimeDilutedSpinColorDiagonal implementation                  *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TTimeDilutedSpinColorDiagonal<FImpl>::TTimeDilutedSpinColorDiagonal(const std::string name)
: Module<NoPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TTimeDilutedSpinColorDiagonal<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TTimeDilutedSpinColorDiagonal<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTimeDilutedSpinColorDiagonal<FImpl>::setup(void)
{
    envCreateDerived(DilutedNoise<FImpl>, 
                     TimeDilutedSpinColorDiagonalNoise<FImpl>,
                     getName(), 1, env().getGrid());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTimeDilutedSpinColorDiagonal<FImpl>::execute(void)
{
    auto &noise = envGet(DilutedNoise<FImpl>, getName());

    noise.generateNoise(*env().get4dRng());
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MNoise_TimeDilutedSpinColorDiagonal_hpp_
