#ifndef Hadrons_MScalarSUN_TimeMomProbe_hpp_
#define Hadrons_MScalarSUN_TimeMomProbe_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/Modules/MScalarSUN/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *          n-point functions O(t,p)*tr(phi(t_1,p_1)*...*phi(t_n,p_n))        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class TimeMomProbePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TimeMomProbePar,
                                    std::string,              field,
                                    std::vector<std::string>, op,
                                    std::vector<std::vector<std::string>>, timeMom,
                                    std::string,              output);
};

class TimeMomProbeResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TimeMomProbeResult,
                                    std::string,                   op,
                                    std::vector<std::vector<int>>, timeMom,
                                    std::vector<Complex>,          data);
};

template <typename SImpl>
class TTimeMomProbe: public Module<TimeMomProbePar>
{
public:
    typedef typename SImpl::Field                    Field;
    typedef typename SImpl::SiteField::scalar_object Site;
    typedef typename SImpl::ComplexField             ComplexField;
    typedef          std::vector<Complex>            SlicedOp;
public:
    // constructor
    TTimeMomProbe(const std::string name);
    // destructor
    virtual ~TTimeMomProbe(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void vectorModulo(std::vector<int> &v);
};

MODULE_REGISTER_TMP(TimeMomProbeSU2, TTimeMomProbe<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(TimeMomProbeSU3, TTimeMomProbe<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(TimeMomProbeSU4, TTimeMomProbe<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(TimeMomProbeSU5, TTimeMomProbe<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(TimeMomProbeSU6, TTimeMomProbe<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                        TTimeMomProbe implementation                        *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TTimeMomProbe<SImpl>::TTimeMomProbe(const std::string name)
: Module<TimeMomProbePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TTimeMomProbe<SImpl>::getInput(void)
{
    std::vector<std::string> in = par().op;
    
    in.push_back(par().field);

    return in;
}

template <typename SImpl>
std::vector<std::string> TTimeMomProbe<SImpl>::getOutput(void)
{
    std::vector<std::string> out;
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTimeMomProbe<SImpl>::setup(void)
{
    envTmpLat(ComplexField, "ftBuf");
    envTmpLat(Field, "ftMatBuf");
}

// execution ///////////////////////////////////////////////////////////////////
// NB: time is direction 0
template <typename SImpl>
void TTimeMomProbe<SImpl>::vectorModulo(std::vector<int> &v)
{
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        auto d = env().getDim(mu);
        v[mu] = ((v[mu] % d) + d) % d;
    }
}

template <typename SImpl>
void TTimeMomProbe<SImpl>::execute(void)
{
    const unsigned int                           nd = env().getNd();
    const unsigned int                           nt = env().getDim(0);
    double                                       partVol = 1.;
    std::set<std::vector<int>>                   timeMomSet;
    std::vector<std::vector<std::vector<int>>>   timeMom;
    std::vector<std::vector<int>>                transferMom;
    FFT                                          fft(env().getGrid());
    std::vector<int>                             dMask(nd, 1);
    std::vector<TimeMomProbeResult>              result;
    std::map<std::string, std::vector<SlicedOp>> slicedOp;
    std::vector<SlicedOp>                        slicedProbe;
    auto                                         &phi = envGet(Field, par().field);

    envGetTmp(ComplexField, ftBuf);
    envGetTmp(Field, ftMatBuf);
    dMask[0] = 0;
    for (unsigned int mu = 1; mu < nd; ++mu)
    {
        partVol *= env().getDim(mu);
    }
    timeMom.resize(par().timeMom.size());
    for (unsigned int p = 0; p < timeMom.size(); ++p)
    {
        for (auto &tms: par().timeMom[p])
        {
            std::vector<int> tm = strToVec<int>(tms);
            
            timeMom[p].push_back(tm);
            timeMomSet.insert(tm);
        }
        transferMom.push_back(std::vector<int>(nd - 1, 0));
        for (auto &tm: timeMom[p])
        {
            for (unsigned int j = 1; j < nd; ++j)
            {
                transferMom[p][j - 1] -= tm[j];
            }
        }
        LOG(Message) << "Probe " << p << " (" << timeMom[p].size() << " points) : " << std::endl;
        LOG(Message) << "  phi(t_i, p_i) for (t_i, p_i) in " << timeMom[p] << std::endl;
        LOG(Message) << "  operator with momentum " << transferMom[p] << std::endl;
    }
    LOG(Message) << "FFT: field '" << par().field << "'" << std::endl;
    fft.FFT_dim_mask(ftMatBuf, phi, dMask, FFT::forward);
    slicedProbe.resize(timeMom.size());
    for (unsigned int p = 0; p < timeMom.size(); ++p)
    {
        std::vector<int> qt;

        LOG(Message) << "Making probe " << p << std::endl;
        slicedProbe[p].resize(nt);
        for (unsigned int t = 0; t < nt; ++t)
        {
            Site acc;
            
            for (unsigned int i = 0; i < timeMom[p].size(); ++i)
            {
                Site buf;

                qt     = timeMom[p][i];
                qt[0] += t;
                vectorModulo(qt);
                peekSite(buf, ftMatBuf, qt);
                if (i == 0)
                {
                    acc = buf;
                }
                else
                {
                    acc *= buf;
                }
            }
            slicedProbe[p][t] = TensorRemove(trace(acc));
        }
        //std::cout << slicedProbe[p]<< std::endl;
    }
    for (auto &o: par().op)
    {
        auto &op = envGet(ComplexField, o);

        slicedOp[o].resize(transferMom.size());
        LOG(Message) << "FFT: operator '" << o << "'" << std::endl;
        fft.FFT_dim_mask(ftBuf, op, dMask, FFT::forward);
        //std::cout << ftBuf << std::endl;
        for (unsigned int p = 0; p < transferMom.size(); ++p)
        {
            std::vector<int> qt(nd, 0);

            for (unsigned int j = 1; j < nd; ++j)
            {
                qt[j] = transferMom[p][j - 1];
            }
            slicedOp[o][p].resize(nt);
            for (unsigned int t = 0; t < nt; ++t)
            {
                TComplex buf;

                qt[0] = t;
                vectorModulo(qt);
                peekSite(buf, ftBuf, qt);
                slicedOp[o][p][t] = TensorRemove(buf);
            }
            //std::cout << ftBuf << std::endl;
            //std::cout << slicedOp[o][p] << std::endl;
        }
    }
    LOG(Message) << "Making correlators" << std::endl;
    for (auto &o: par().op)
    for (unsigned int p = 0; p < timeMom.size(); ++p)
    {
        TimeMomProbeResult r;

        LOG(Message) << "  <" << o << " probe_" << p << ">" << std::endl;
        r.op      = o;
        r.timeMom = timeMom[p];
        r.data    = makeTwoPoint(slicedOp[o][p], slicedProbe[p], 1./partVol);
        result.push_back(r);
    }
    saveResult(par().output, "timemomprobe", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_TimeMomProbe_hpp_
