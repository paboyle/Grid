#ifndef Hadrons_LanczosUtils_hpp_
#define Hadrons_LanczosUtils_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/algorithms/iterative/LocalCoherenceLanczos.h>

BEGIN_HADRONS_NAMESPACE

// Lanczos type
#ifndef HADRONS_DEFAULT_LANCZOS_NBASIS
#define HADRONS_DEFAULT_LANCZOS_NBASIS 60
#endif

template <typename T>
struct EigenPack
{
    typedef T VectorType;
    std::vector<RealD> eval;
    std::vector<T>     evec;
    
    EigenPack(void) = default;

    EigenPack(const size_t size, GridBase *grid)
    {
        resize(size, grid);
    }

    void resize(const size_t size, GridBase *grid)
    {
        eval.resize(size);
        evec.resize(size, grid);
    }

    void read(const std::string fileStem)
    {
        std::string     evecFilename = fileStem + "_evec.bin";
        std::string     evalFilename = fileStem + "_eval.xml";
        emptyUserRecord record;
        ScidacReader    binReader;
        XmlReader       xmlReader(evalFilename);

        LOG(Message) << "Reading " << evec.size() << " eigenvectors from '" 
                     << evecFilename << "'" << std::endl;
        binReader.open(evecFilename);
        for(int k = 0; k < evec.size(); ++k) 
        {
            binReader.readScidacFieldRecord(evec[k], record);
        }
        binReader.close();
        LOG(Message) << "Reading " << eval.size() << " eigenvalues from '" 
                     << evalFilename << "'" << std::endl;
        Grid::read(xmlReader, "evals", eval);
    }

    void write(const std::string fileStem)
    {
        std::string     evecFilename = fileStem + "_evec.bin";
        std::string     evalFilename = fileStem + "_eval.xml";
        emptyUserRecord record;
        ScidacWriter    binWriter;
        XmlWriter       xmlWriter(evalFilename);

        LOG(Message) << "Writing " << evec.size() << " eigenvectors to '" 
                     << evecFilename << "'" << std::endl;
        binWriter.open(fileStem + "_evec.bin");
        for(int k = 0; k < evec.size(); ++k) 
        {
            binWriter.writeScidacFieldRecord(evec[k], record);
        }
        binWriter.close();
        LOG(Message) << "Writing " << eval.size() << " eigenvalues to '" 
                     << evalFilename << "'" << std::endl;
        Grid::write(xmlWriter, "evals", eval);
    }
};

template <typename FImpl>
using FineEigenPack = EigenPack<typename FImpl::FermionField>;

template <typename FImpl, int nBasis>
using CoarseEigenPack = EigenPack<
    typename LocalCoherenceLanczos<typename FImpl::SiteSpinor, 
                                   typename FImpl::SiteComplex, 
                                   nBasis>::CoarseField>;

END_HADRONS_NAMESPACE

#endif // Hadrons_LanczosUtils_hpp_