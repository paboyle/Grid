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

const int NCONFIGS = 10;
const int NWILSON = 10;

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

    Complex imag_unit(0, 1);

    Real wlA;
    std::vector<Real> logWlAvg(NWILSON, 0.0), logWlTime(NWILSON, 0.0), logWlSpace(NWILSON, 0.0);

    pRNG.SeedRandomDevice();

    LOG(Message) << "Wilson loop calculation beginning" << std::endl;
    for(int ic = 0; ic < NCONFIGS; ic++){
        LOG(Message) << "Configuration " << ic <<std::endl;
        photon.StochasticField(a, pRNG);

        // Exponentiate photon field
        expA = exp(imag_unit*a);

        // Calculate Wilson loops
        for(int iw=1; iw<=NWILSON; iw++){
            wlA = NewWilsonLoops<QedPeriodicGimplR>::avgWilsonLoop(expA, iw, iw) * 3;
            logWlAvg[iw-1] -= 2*log(wlA);
            wlA = NewWilsonLoops<QedPeriodicGimplR>::avgTimelikeWilsonLoop(expA, iw, iw) * 3;
            logWlTime[iw-1] -= 2*log(wlA);
            wlA = NewWilsonLoops<QedPeriodicGimplR>::avgSpatialWilsonLoop(expA, iw, iw) * 3;
            logWlSpace[iw-1] -= 2*log(wlA);
        }
    }
    LOG(Message) << "Wilson loop calculation completed" << std::endl;
    
    // Calculate Wilson loops
    for(int iw=1; iw<=10; iw++){
        LOG(Message) << iw << 'x' << iw << " Wilson loop" << std::endl;
        LOG(Message) << "-2log(W) average: " << logWlAvg[iw-1]/NCONFIGS << std::endl;
        LOG(Message) << "-2log(W) timelike: " << logWlTime[iw-1]/NCONFIGS << std::endl;
        LOG(Message) << "-2log(W) spatial: " << logWlSpace[iw-1]/NCONFIGS << std::endl;
    }

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
