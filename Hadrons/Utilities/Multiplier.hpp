#ifndef Hadrons_Multiplier_hpp_
#define Hadrons_Multiplier_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

namespace Multiplier
{
    class TrajRange: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(TrajRange,
                                        unsigned int, start,
                                        unsigned int, end,
                                        unsigned int, step,
                                        std::string,  exclude);

        inline std::vector<unsigned int> getTrajectoryList(void)
        {
            std::vector<unsigned int> excVec = strToVec<unsigned int>(exclude);
            std::vector<unsigned int> list;

            for (unsigned int t = start; t < end; t += step)
            {
                if (std::find(excVec.begin(), excVec.end(), t) == excVec.end())
                {
                    list.push_back(t);
                }
            }

            return list;
        }
    };
    
    class GlobalPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(GlobalPar,
                                        TrajRange, trajCounter,
                                        unsigned int, nt,
                                        std::string, diskVectorDir,
                                        unsigned int, blockSize,
                                        unsigned int, cacheSize,
                                        std::string, output);
    };

    class A2AMatrixPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMatrixPar,
                                        std::string, file,
                                        std::string, dataset,
                                        std::string, name);
    };

    class ProductPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(ProductPar,
                                        std::string, terms,
                                        std::vector<std::string>, times);
    };

    class MultiplierMetadata : Serializable
    {
      public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(MultiplierMetadata,
                                        std::vector<Real>, momentum,
                                        std::vector<Gamma::Algebra>, gammas);
    };
}

END_HADRONS_NAMESPACE

#endif // Hadrons_Multiplier_hpp_
