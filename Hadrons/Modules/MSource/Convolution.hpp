#ifndef Hadrons_MSource_Convolution_hpp_
#define Hadrons_MSource_Convolution_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Convolution                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class ConvolutionPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ConvolutionPar,
                                    std::string, field,
                                    std::string, filter,
                                    std::string, mom);
};

template <typename FImpl>
class TConvolution: public Module<ConvolutionPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TConvolution(const std::string name);
    // destructor
    virtual ~TConvolution(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::vector<int> mom_;
};

MODULE_REGISTER_TMP(Convolution, TConvolution<FIMPL>, MSource);

/******************************************************************************
 *                 TConvolution implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TConvolution<FImpl>::TConvolution(const std::string name)
: Module<ConvolutionPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TConvolution<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().field, par().filter};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TConvolution<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TConvolution<FImpl>::setup(void)
{
     mom_ = strToVec<int>(par().mom);
     if(mom_.size() != env().getNd()) {
         HADRONS_ERROR(Size, std::string("momentum has ")
                 + std::to_string(mom_.size()) + " instead of "
                 + std::to_string(env().getNd()) + " components");
     }

    envCreateLat(PropagatorField, getName());
    envTmpLat(ComplexField, "momfield");
    envTmp(FFT, "fft", 1, env().getGrid());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TConvolution<FImpl>::execute(void)
{
    auto &filter = envGet(ComplexField, par().filter);
    auto &field  = envGet(PropagatorField, par().field);
    auto &out    = envGet(PropagatorField, getName());
    envGetTmp(ComplexField, momfield);
    envGetTmp(FFT, fft);

    std::vector<int> mask(env().getNd(), 1);
    mask.back()=0; //transform only the spatial dimensions

    startTimer("Fourier transform");
    fft.FFT_dim_mask(momfield, filter, mask, FFT::forward);
    fft.FFT_dim_mask(out,      field, mask, FFT::forward);
    stopTimer("Fourier transform");

    startTimer("momentum-space multiplication");
    out=momfield*out;
    stopTimer("momentum-space multiplication");

    startTimer("inserting momentum");
    for(int mu=0; mu<env().getNd(); mu++)
    {
       if(mom_[mu]!=0)
       {
          out=Cshift(out, mu, -mom_[mu]);
       }
    }
    stopTimer("inserting momentum");

    startTimer("Fourier transform");
    fft.FFT_dim_mask(out, out, mask, FFT::backward);
    stopTimer("Fourier transform");
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_Convolution_hpp_
