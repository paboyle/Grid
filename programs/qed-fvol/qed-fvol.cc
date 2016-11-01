#include <Global.hpp>

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

typedef QedGimpl<vComplex>      QedGimplR;
typedef Photon<QedGimplR>       PhotonR;
typedef PhotonR::GaugeField     EmField;
typedef PhotonR::GaugeLinkField EmComp;

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

    pRNG.SeedRandomDevice();
    photon.StochasticField(a, pRNG);

    // Calculate log of plaquette
    EmComp              plaqA(&grid);
    EmComp              wlA(&grid);
    EmComp              tmp(&grid);
    std::vector<EmComp> a_comp(4, &grid);

    for (int dir = 0; dir < Nd; dir++) {
      a_comp[dir] = PeekIndex<LorentzIndex>(a, dir);
    }

    plaqA = zero;
    wlA = zero;

    for(int mu = 1; mu < Nd; mu++) {
        for(int nu = 0; nu < mu; nu++) {
            tmp = a_comp[mu] + Cshift(a_comp[nu], mu, 1) - Cshift(a_comp[mu], nu, 1) - a_comp[nu];
            plaqA = plaqA + cos(tmp);

            tmp = a_comp[mu] + Cshift(a_comp[mu], mu, 1)
                  + Cshift(a_comp[nu], mu, 2) + Cshift(Cshift(a_comp[nu], mu, 2), nu, 1)
                  - Cshift(Cshift(a_comp[mu], nu, 2), mu, 1) - Cshift(a_comp[mu], nu, 2)
                  - Cshift(a_comp[nu], nu, 1) - a_comp[nu];
            wlA = wlA + cos(tmp);
        }
    }

    Real vol = grid.gSites();
    Real faces = (1.0 * Nd * (Nd - 1)) / 2.0;

    Complex avgPlaqA = sum(trace(plaqA));
    avgPlaqA = avgPlaqA / vol / faces;

    Complex avgWlA = sum(trace(wlA));
    avgWlA = avgWlA / vol / faces;

    TComplex tplaqsite;
    LatticeComplex plaqtrace = trace(plaqA);
    std::vector<int> site0 = {0,0,0,0};
    peekSite(tplaqsite, plaqtrace, site0);
    Complex plaqsite = TensorRemove(tplaqsite);

    LOG(Message) << "Plaquette average: " << avgPlaqA << std::endl;
    LOG(Message) << "2x2 Wilson Loop average: " << avgWlA << std::endl;
    LOG(Message) << "Plaquette (one site): " << plaqsite / faces << std::endl;

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
