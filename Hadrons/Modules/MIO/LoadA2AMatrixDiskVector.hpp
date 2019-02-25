#ifndef Hadrons_MIO_LoadA2AMatrixDiskVector_hpp_
#define Hadrons_MIO_LoadA2AMatrixDiskVector_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DiskVector.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         LoadA2AMatrixDiskVector                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadA2AMatrixDiskVectorPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadA2AMatrixDiskVectorPar,
                                    std::string,  dirstem,
                                    std::string,  dataset,
                                    unsigned int, ni,
                                    unsigned int, nj);
};

template <typename FImpl>
class TLoadA2AMatrixDiskVector: public Module<LoadA2AMatrixDiskVectorPar>
{
  public:
    FERM_TYPE_ALIASES(FImpl, );
    // constructor
    TLoadA2AMatrixDiskVector(const std::string name);
    // destructor
    virtual ~TLoadA2AMatrixDiskVector(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadA2AMatrixDiskVector, TLoadA2AMatrixDiskVector<FIMPL>, MIO);

/******************************************************************************
 *                 TLoadA2AMatrixDiskVector implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadA2AMatrixDiskVector<FImpl>::TLoadA2AMatrixDiskVector(const std::string name)
: Module<LoadA2AMatrixDiskVectorPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadA2AMatrixDiskVector<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TLoadA2AMatrixDiskVector<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadA2AMatrixDiskVector<FImpl>::setup(void)
{
    int nt = env().getDim(Tp);
    std::string dataset = par().dataset;
    std::string dvFile = getName() + "_" + dataset + "_DV_" + std::to_string(vm().getTrajectory());
    GridBase *grid = env().getGrid();

    envCreate(EigenDiskVector<ComplexD>, getName(), 1, dvFile, nt, 1, true, grid);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadA2AMatrixDiskVector<FImpl>::execute(void)
{
    int nt = env().getDim(Tp);
    int ni = par().ni;
    int nj = par().nj;
    std::string dirstem  = par().dirstem;
    std::string dataset  = par().dataset;

    GridBase *grid = env().getGrid();

    LOG(Message) << "Meson Field size: " << nt << "*" << ni << "*" << nj
                 << " (filesize " << sizeString(nt * ni * nj * sizeof(Complex))
                 << "/momentum/bilinear)" << std::endl;

    auto &mesonFieldDV = envGet(EigenDiskVector<ComplexD>, getName());
    std::string file = dirstem + "." + std::to_string(vm().getTrajectory()) + "/" + dataset + ".h5";

    double t;
    LOG(Message) << "-- Loading '" << file << "'..." << std::endl;
    A2AMatrixIo<HADRONS_A2AM_IO_TYPE> mfIO(file, dataset, nt, ni, nj);
    mfIO.load(mesonFieldDV, &t, grid);
    LOG(Message) << "Read " << mfIO.getSize() << " bytes in " << t << " usec, " << mfIO.getSize() / t * 1.0e6 / 1024 / 1024 << " MB/s" << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadA2AMatrixDiskVector_hpp_
