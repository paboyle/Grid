#ifndef Hadrons_MSolver_Guesser_hpp_
#define Hadrons_MSolver_Guesser_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MSolver)

template <typename FImpl, int nBasis>
std::shared_ptr<LinearFunction<typename FImpl::FermionField>> 
makeGuesser(const std::string epackName)
{
    typedef typename FImpl::FermionField                  FermionField;
    typedef FermionEigenPack<FImpl>                       EPack;
    typedef CoarseFermionEigenPack<FImpl, nBasis>         CoarseEPack;
    typedef DeflatedGuesser<FermionField>                 FineGuesser;
    typedef LocalCoherenceDeflatedGuesser<
        FermionField, typename CoarseEPack::CoarseField>  CoarseGuesser;

    std::shared_ptr<LinearFunction<typename FImpl::FermionField>> guesserPt;

    DEFINE_ENV_LAMBDA;

    if (epackName.empty())
    {
        guesserPt.reset(new ZeroGuesser<FermionField>());
    }
    else
    {
        try
        {
            auto &epack = envGetDerived(EPack, CoarseEPack, epackName);
            
            LOG(Message) << "using low-mode deflation with coarse eigenpack '"
                         << epackName << "' (" 
                         << epack.evecCoarse.size() << " modes)" << std::endl;
            guesserPt.reset(new CoarseGuesser(epack.evec, epack.evecCoarse,
                                              epack.evalCoarse));
        }
        catch (Exceptions::ObjectType &e)
        {
            auto &epack = envGet(EPack, epackName);

            LOG(Message) << "using low-mode deflation with eigenpack '"
                         << epackName << "' (" 
                         << epack.evec.size() << " modes)" << std::endl;
            guesserPt.reset(new FineGuesser(epack.evec, epack.eval));
        }
    }

    return guesserPt;
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif