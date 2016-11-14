#include <Global.hpp>
#include <WilsonLoops.h>

using namespace Grid;
using namespace QCD;
using namespace QedFVol;

template <class S> 
class QedGimpl 
{
public:
  typedef S Simd;

  template <typename vtype>
  using iImplGaugeLink  = iScalar<iScalar<iScalar<vtype>>>;
  template <typename vtype>
  using iImplGaugeField = iVector<iScalar<iScalar<vtype>>, Nd>;

  typedef iImplGaugeLink<Simd> SiteGaugeLink;
  typedef iImplGaugeField<Simd> SiteGaugeField;

  typedef Lattice<SiteGaugeLink> GaugeLinkField; // bit ugly naming; polarised
                                                 // gauge field, lorentz... all
                                                 // ugly
  typedef Lattice<SiteGaugeField> GaugeField;
};

typedef QedGimpl<vComplex>              QedGimplR;
typedef PeriodicGaugeImpl<QedGimplR>    QedPeriodicGimplR;
typedef Photon<QedGimplR>               PhotonR;
typedef PhotonR::GaugeField             EmField;
typedef PhotonR::GaugeLinkField         EmComp;

int main(int argc, char *argv[])
{
    // parse command line
    std::string parameterFileName;
    
    if (argc < 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file> [Grid options]";
        std::cerr << std::endl;
        std::exit(EXIT_FAILURE);
    }
    parameterFileName = argv[1];
    
    // initialization
    Grid_init(&argc, &argv);
    QedFVolLogError.Active(GridLogError.isActive());
    QedFVolLogWarning.Active(GridLogWarning.isActive());
    QedFVolLogMessage.Active(GridLogMessage.isActive());
    QedFVolLogIterative.Active(GridLogIterative.isActive());
    QedFVolLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;
    
    // QED stuff
    std::vector<int> latt_size   = GridDefaultLatt();
    std::vector<int> simd_layout = GridDefaultSimd(4, vComplex::Nsimd());
    std::vector<int> mpi_layout  = GridDefaultMpi();
    GridCartesian    grid(latt_size,simd_layout,mpi_layout);
    GridParallelRNG  pRNG(&grid);
    PhotonR          photon(PhotonR::Gauge::Feynman,
                            PhotonR::ZmScheme::QedL);
    EmField          a(&grid);
    EmField          expA(&grid);

    Real wlA, logWlA;

    pRNG.SeedRandomDevice();
    photon.StochasticField(a, pRNG);

    // Exponentiate photon field
    Complex imag_unit(0, 1);
    expA = exp(imag_unit*a);

    // Calculate Wilson loops
    for(int i=1; i<=10; i++){
        LOG(Message) << i << 'x' << i << " Wilson loop" << std::endl;
        wlA = NewWilsonLoops<QedPeriodicGimplR>::avgWilsonLoop(expA, i, i) * 3;
        logWlA = -2*log(wlA);
        LOG(Message) << "-2log(W) average: " << logWlA << std::endl;
        wlA = NewWilsonLoops<QedPeriodicGimplR>::avgTimelikeWilsonLoop(expA, i, i) * 3;
        logWlA = -2*log(wlA);
        LOG(Message) << "-2log(W) timelike: " << logWlA << std::endl;
        wlA = NewWilsonLoops<QedPeriodicGimplR>::avgSpatialWilsonLoop(expA, i, i) * 3;
        logWlA = -2*log(wlA);
        LOG(Message) << "-2log(W) spatial: " << logWlA << std::endl;
    }

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
