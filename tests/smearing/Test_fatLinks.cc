/* 
 * Test_fatLinks.cc                                                               
 * 
 * D. Clarke 
 * 
 * Test the various constructs used to make fat links. 
 * 
 */


#include <Grid/Grid.h>
using namespace Grid;

template <class Gimpl> class : public Gimpl {
public:

    INHERIT_GIMPL_TYPES(Gimpl);

    typedef typename Gimpl::GaugeLinkField GaugeMat;
    typedef typename Gimpl::GaugeField GaugeLorentz;

    static void staple(GaugeMat &plaq, const std::vector<GaugeMat> &U, const int mu, const int nu) {
    }

}

int main (int argc, char **argv)
{
    Grid_init(&argc, &argv);

    Coordinate simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();
    Coordinate latt_size   = GridDefaultLatt();

    GridCartesian spacetime(latt_size,simd_layout,mpi_layout);

    PeriodicGimplD::Field U(&spacetime);

    Grid_finalize();
}