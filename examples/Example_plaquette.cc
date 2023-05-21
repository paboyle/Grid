/* 
 * Example_plaquette.cc                                                               
 * 
 * D. Clarke 
 * 
 * Here I just want to create an incredibly simple main to get started with GRID and get used
 * to its syntax. If the reader is like me, they vaguely understand something about lattice coding,
 * they don't know a ton of C++, don't know much of the fine details, and certainly know nothing about GRID.
 *
 * Once you've made a new executable, like this one, you can bootstrap.sh again. At this point,
 * the code should be able to find your new executable. You can tell that bootstrap.sh worked by
 * having a look at Make.inc. You should see your executable inside there.
 *
 * Warning: This code illustrative only, not well tested, and not meant for production use. The best
 * way to read this code is to start at the main.
 * 
 */


// All your mains should have this
#include <Grid/Grid.h>
using namespace Grid;


// This copies what already exists in WilsonLoops.h. The point here is to be pedagogical and explain in
// detail what everything does so we can see how GRID works.
template <class Gimpl> class WLoops : public Gimpl {
public:
    // Gimpl seems to be an arbitrary class. Within this class, it is expected that certain types are
    // already defined, things like Scalar and Field. This macro includes a bunch of #typedefs that
    // implement this equivalence at compile time.
    // WARNING: The first time you include this or take it out, the compile time will increase a lot.
    INHERIT_GIMPL_TYPES(Gimpl);

    // Some example Gimpls can be found in GaugeImplementations.h, at the bottom. These are in turn built
    // out of GaugeImplTypes, which can be found in GaugeImplTypes.h. The GaugeImplTypes contain the base
    // field/vector/link/whatever types. These inherit from iScalar, iVector, and iMatrix objects, which
    // are sort of the building blocks for gerenal math objects. The "i" at the beginning of these names
    // indicates that they should be for internal use only. It seems like these base types have the
    // acceleration, e.g. SIMD or GPU or what-have-you, abstracted away. How you accelerate these things
    // appears to be controlled through a template parameter called vtype.

    // The general math/physics objects, such as a color matrix, are built up by nesting these objects.
    // For instance a general color matrix has two color indices, so it's built up like
    //     iScalar<iScalar<iMatrix<vtype ...
    // where the levels going from the inside out are color, spin, then Lorentz indices. Scalars have
    // no indices, so it's what we use when such an index isn't needed. Lattice objects are made by one
    // higher level of indexing using iVector.

    // These types will be used for U and U_mu objects, respectively.
    typedef typename Gimpl::GaugeLinkField GaugeMat;
    typedef typename Gimpl::GaugeField GaugeLorentz;

    // U_mu_nu(x)
    static void dirPlaquette(GaugeMat &plaq, const std::vector<GaugeMat> &U, const int mu, const int nu) {
        // These CovShift calls seem to carry out the multiplication already. A positive shift moves the lattice 
        // site x_mu = 1 in the RHS to x_mu = 0 in the result.
        plaq = Gimpl::CovShiftForward(U[mu],mu,
                    Gimpl::CovShiftForward(U[nu],nu,
                        Gimpl::CovShiftBackward(U[mu],mu,
                            Gimpl::CovShiftIdentityBackward(U[nu], nu))));
    }

    // tr U_mu_nu(x)
    static void traceDirPlaquette(ComplexField &plaq, const std::vector<GaugeMat> &U, const int mu, const int nu) {
        // This .Grid() syntax seems to get the pointer to the GridBase. Apparently this is needed as argument
        // to instantiate a Lattice object.
        GaugeMat sp(U[0].Grid());
        dirPlaquette(sp, U, mu, nu);
        plaq = trace(sp);
    }

    // sum_mu_nu tr U_mu_nu(x)
    static void sitePlaquette(ComplexField &Plaq, const std::vector<GaugeMat> &U) {
        ComplexField sitePlaq(U[0].Grid());
        Plaq = Zero();
        // Nd=4 and Nc=3 are set as global constants in QCD.h
        for (int mu = 1; mu < Nd; mu++) {
            for (int nu = 0; nu < mu; nu++) {
                traceDirPlaquette(sitePlaq, U, mu, nu);
                Plaq = Plaq + sitePlaq;
            }
        }
    }

    // sum_mu_nu_x Re tr U_mu_nu(x)
    static RealD sumPlaquette(const GaugeLorentz &Umu) {
        std::vector<GaugeMat> U(Nd, Umu.Grid());
        for (int mu = 0; mu < Nd; mu++) {
            // Umu is a GaugeLorentz object, and as such has a non-trivial Lorentz index. We can
            // access the element in the mu Lorentz index with this PeekIndex syntax.
            U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
        }
        ComplexField Plaq(Umu.Grid());
        sitePlaquette(Plaq, U);
        // I guess this should be the line that sums over all space-time sites.
        auto Tp = sum(Plaq);
        // Until now, we have been working with objects inside the tensor nest. This TensorRemove gets
        // rid of the tensor nest to return whatever is inside.
        auto p  = TensorRemove(Tp);
        return p.real();
    }

    // < Re tr U_mu_nu(x) >
    static RealD avgPlaquette(const GaugeLorentz &Umu) {
        // Real double type
        RealD sumplaq = sumPlaquette(Umu);
        // gSites() is the number of global sites. there is also lSites() for local sites.
        double vol = Umu.Grid()->gSites();
        // The number of orientations. 4*3/2=6 for Nd=4, as known.
        double faces = (1.0 * Nd * (Nd - 1)) / 2.0;
        return sumplaq / vol / faces / Nc;
    }
};


// Next we show an example of how to construct an input parameter class. We first inherit
// from Serializable. Then all class data members have to be defined using the
// GRID_SERIALIZABLE_CLASS_MEMBERS macro. This variadic macro allows for arbitrarily many
// class data members. In the below case, we make a parameter file holding the configuration
// name. Here, it expects the name to be labeled with "conf_name" in the configuration file. 
struct ConfParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(
        ConfParameters,
        std::string, conf_name);

    template <class ReaderClass>
    ConfParameters(Reader<ReaderClass>& Reader){
        // If we are reading an XML file, it should be structured like:
        // <grid>
        //   <parameters>
        //     <conf_name>l20t20b06498a_nersc.302500</conf_name>
        //   </parameters>
        // </grid>
        read(Reader, "parameters", *this);
    }
};



// This syntax lets you pass command line arguments to main. An asterisk means that what follows is
// a pointer. Two asterisks means what follows is a pointer to an array. 
int main (int argc, char **argv)
{
    // This initializes Grid. Some command line options include
    //   --mpi n.n.n.n
    //   --threads n
    //   --grid n.n.n.n
    Grid_init(&argc, &argv);

    // This is where you would specify a custom lattice size, if not from the command line. Here
    // Nd is a global quantity that is currently set to 4.
    Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();
    Coordinate latt_size   = GridDefaultLatt();

    // Instantiate the spacetime Grid on which everything will be built.
    GridCartesian GRID(latt_size,simd_layout,mpi_layout);

    // The PeriodicGimplD type is what you want for gauge matrices. There is also a LatticeGaugeFieldD
    // type that you can use, which will work perfectly with what follows. 
    PeriodicGimplD::Field U(&GRID);

    // Here we read in the parameter file params.json to get conf_name. The last argument is what the
    // top organizational level is called in the param file. 
    XmlReader Reader("Example_plaquette.xml",false, "grid");
    ConfParameters param(Reader);  

    // Load a lattice from SIMULATeQCD into U. SIMULATeQCD finds plaquette = 0.6381995717
    FieldMetaData header;
    NerscIO::readConfiguration(U, header, param.conf_name);

    // Let's see what we find.
    RealD plaq = WLoops<PeriodicGimplD>::avgPlaquette(U);

    // This is how you make log messages.
    std::cout << GridLogMessage << std::setprecision(std::numeric_limits<Real>::digits10 + 1) << "Plaquette = " << plaq << std::endl;

    // To wrap things up.
    Grid_finalize();
}