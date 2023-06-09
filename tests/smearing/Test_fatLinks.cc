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
    Grid_init(&argc,&argv);

    Coordinate latt_size   = GridDefaultLatt();
    Coordinate simd_layout = GridDefaultSimd(Nd,vComplexD::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();

    Grid_log("mpi = ",mpi_layout);
    Grid_log("simd = ",simd_layout);
    Grid_log("latt = ",latt_size);

    GridCartesian GRID(latt_size,simd_layout,mpi_layout);

    XmlReader Reader("fatParams.xml",false, "grid");
    fatParams param(Reader);  

    LatticeGaugeField Umu(&GRID);
    LatticeGaugeField U_smr(&GRID);
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, param.conf_in);

    Smear_HISQ_fat hisq_fat(&GRID);

    hisq_fat.smear(U_smr,Umu);

    NerscIO::writeConfiguration(U_smr,param.conf_out,"HISQ");

    Grid_finalize();
}