/* 
 * Test_fatLinks.cc                                                               
 * 
 * D. Clarke 
 * 
 * Test the various constructs used to make fat links. 
 * 
 */

#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>
#include <Grid/qcd/smearing/HISQSmearing.h>
using namespace Grid;


// Make the logger work like Python print()
template<typename ... Args>
inline std::string sjoin(Args&&... args) noexcept {
    std::ostringstream msg;
    (msg << ... << args);
    return msg.str();
}
template <typename ... Args>
inline void Grid_log(Args&&... args) {
    std::string msg = sjoin(std::forward<Args>(args)...);
    std::cout << GridLogMessage << msg << std::endl;
}

struct fatParams: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(
        fatParams,
        std::string, conf_in,
        std::string, conf_out);

    template <class ReaderClass>
    fatParams(Reader<ReaderClass>& Reader){
        read(Reader, "parameters", *this);
    }
};

//
// one method: input --> fat
// another   : input --> long (naik)
// another   : input --> unitarize
//

int main (int argc, char **argv)
{
    // Initialize the Grid
    Grid_init(&argc,&argv);
    Coordinate latt_size   = GridDefaultLatt();
    Coordinate simd_layout = GridDefaultSimd(Nd,vComplexD::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();
    Grid_log("mpi = ",mpi_layout);
    Grid_log("simd = ",simd_layout);
    Grid_log("latt = ",latt_size);
    GridCartesian GRID(latt_size,simd_layout,mpi_layout);

    // Instantiate the LatticeGaugeField objects holding thin (Umu) and fat (U_smr) links
    LatticeGaugeField Umu(&GRID);
    LatticeGaugeField U_smr(&GRID);

    // Read in the parameter file
    XmlReader Reader("fatParams.xml",false, "grid");
    fatParams param(Reader);  
    FieldMetaData header;

    // Read the configuration into Umu
    NerscIO::readConfiguration(Umu, header, param.conf_in);

    // Smear Umu and store result in U_smr
    Smear_HISQ_fat<LatticeGaugeField> hisq_fat(&GRID,1/8.,0.,1/16.,1/64.,1/384.,0.);
    hisq_fat.smear(U_smr,Umu);

    NerscIO::writeConfiguration(U_smr,param.conf_out,"HISQ");

    // Test a C-style instantiation 
    double path_coeff[6] = {1, 2, 3, 4, 5, 6};
    Smear_HISQ_fat<LatticeGaugeField> hisq_fat_Cstyle(&GRID,path_coeff);

    // Make sure result doesn't change w.r.t. a trusted lattice
    NerscIO::readConfiguration(Umu, header, "nersc.l8t4b3360.3link.control");
    LatticeGaugeField diff(&GRID);
    diff = Umu-U_smr;
    auto absDiff = norm2(diff)/norm2(Umu);
    Grid_log(" |Umu-U|/|Umu| =",absDiff);

    Grid_finalize();
}