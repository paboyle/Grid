/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: HMC/ComputeWilsonFlow.cc

Copyright (C) 2017

Author: Guido Cossu <guido.cossu@ed.ac.uk>
Author: Shuhei Yamamoto <syamamoto@bnl.gov>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>
#include <string>

namespace Grid{
  struct WFParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(WFParameters,
            int, steps,
            double, step_size,
            int, meas_interval,
	    double, maxTau, // for the adaptive algorithm
	    int, meas_interval_density,
	    std::string, path); 
       

    template <class ReaderClass >
    WFParameters(Reader<ReaderClass>& Reader){
      read(Reader, "WilsonFlow", *this);
    }

  };

  struct ConfParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(ConfParameters,
	   std::string, conf_path,
           std::string, conf_prefix,
	   std::string, conf_smr_prefix,
           std::string, rng_prefix,
	   int, StartConfiguration,
	   int, EndConfiguration,
           int, Skip);
  
    template <class ReaderClass >
    ConfParameters(Reader<ReaderClass>& Reader){
      read(Reader, "Configurations", *this);
    }

  };
}

template <class T> void writeFile(T& in, std::string const fname){  
#ifdef HAVE_LIME
  // Ref: https://github.com/paboyle/Grid/blob/feature/scidac-wp1/tests/debug/Test_general_coarse_hdcg_phys48.cc#L111
  std::cout << Grid::GridLogMessage << "Writes to: " << fname << std::endl;
  Grid::emptyUserRecord record;
  Grid::ScidacWriter WR(in.Grid()->IsBoss());
  WR.open(fname);
  WR.writeScidacFieldRecord(in,record,0);
  WR.close();
#endif
  // What is the appropriate way to throw error?
}


int main(int argc, char **argv) {
  using namespace Grid;
  
  Grid_init(&argc, &argv);
  GridLogLayout();

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size, simd_layout, mpi_layout);
  
  std::vector<int> seeds({1, 2, 3, 4, 5});
  GridSerialRNG sRNG;
  GridParallelRNG pRNG(&Grid);
  pRNG.SeedFixedIntegers(seeds);

  LatticeGaugeField Umu(&Grid), Uflow(&Grid);
  
  typedef Grid::XmlReader       Serialiser;
  Serialiser Reader("input.xml", false, "root");
  WFParameters WFPar(Reader);
  ConfParameters CPar(Reader);
  CheckpointerParameters CPPar(CPar.conf_path+CPar.conf_prefix, CPar.conf_path+CPar.conf_smr_prefix, CPar.conf_path+CPar.rng_prefix);
  NerscHmcCheckpointer<PeriodicGimplR> CPNersc(CPPar);

  for (int conf = CPar.StartConfiguration; conf <= CPar.EndConfiguration; conf+= CPar.Skip){

  CPNersc.CheckpointRestore(conf, Umu, sRNG, pRNG);

  std::cout << std::setprecision(15);
  std::cout << GridLogMessage << "Initial plaquette: "<< WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu) << std::endl;
  
  std::string file_pre  = WFPar.path;
  std::string file_post = CPar.conf_prefix + "." + std::to_string(conf);

  WilsonFlow<PeriodicGimplR> WF(WFPar.step_size,WFPar.steps,WFPar.meas_interval);
  WF.addMeasurement(WFPar.meas_interval_density, [&file_pre,&file_post,&conf](int step, RealD t, const typename PeriodicGimplR::GaugeField &U){
    
    typedef typename PeriodicGimplR::GaugeLinkField GaugeMat;
    typedef typename PeriodicGimplR::ComplexField ComplexField;
    
    assert(Nd == 4);

    // NOTE:
    // Ideally, turn the folloing into methods of the appropriate class
    /////////////   Compute Energy Density via Clover Leaf    /////////////////////////////////////////////////
    ///// Taken from qcd/smearing/WilsonFlow.h
    //         For plq, use static sitePlaquette from class WilsonLoops in Grid/qcd/utils/WilsonLoops.h and divide it by #faces=(1.0 * Nd * (Nd - 1)) / 2.0, ncol=3
    //E = 1/2 tr( F_munu F_munu )
    //However as  F_numu = -F_munu, only need to sum the trace of the squares of the following 6 field strengths:
    //F_01 F_02 F_03   F_12 F_13  F_23
    GaugeMat F(U.Grid());
    //LatticeComplexD R(U.Grid());
    ComplexField R(U.Grid());
    R = Zero();
  
    for(int mu=0;mu<3;mu++){
      for(int nu=mu+1;nu<4;nu++){
	WilsonLoops<PeriodicGimplR>::FieldStrength(F, U, mu, nu);
	R = R + trace(F*F);
      }
    }
    R = (-1.0) * R;
    
    //// Taken from qcd/utils/WilsonLoops.h
    
    // Bx = -iF(y,z), By = -iF(z,y), Bz = -iF(x,y)
    GaugeMat Bx(U.Grid()), By(U.Grid()), Bz(U.Grid());
    WilsonLoops<PeriodicGimplR>::FieldStrength(Bx, U, Ydir, Zdir);
    WilsonLoops<PeriodicGimplR>::FieldStrength(By, U, Zdir, Xdir);
    WilsonLoops<PeriodicGimplR>::FieldStrength(Bz, U, Xdir, Ydir);

    // Ex = -iF(t,x), Ey = -iF(t,y), Ez = -iF(t,z)
    GaugeMat Ex(U.Grid()), Ey(U.Grid()), Ez(U.Grid());
    WilsonLoops<PeriodicGimplR>::FieldStrength(Ex, U, Tdir, Xdir);
    WilsonLoops<PeriodicGimplR>::FieldStrength(Ey, U, Tdir, Ydir);
    WilsonLoops<PeriodicGimplR>::FieldStrength(Ez, U, Tdir, Zdir);

    double coeff = 8.0/(32.0*M_PI*M_PI);
    ComplexField qfield = coeff*trace(Bx*Ex + By*Ey + Bz*Ez);
    //ComplexField qfield Plq(U.Grid());
    //WilsonLoops<PeriodicGimplR>::sitePlaquette(Plq, U);
    //double coeff = 2.0 / (1.0 * Nd * (Nd - 1)) / 3.0;
    //Plq = coeff * Plq;

    int tau = std::round(t);
    std::string efile = file_pre + "E_dnsty_" + std::to_string(tau) + "_" + file_post;
    writeFile(R,efile);
    std::string tfile = file_pre + "Top_dnsty_" + std::to_string(tau) + "_" + file_post;
    writeFile(qfield,tfile);

    RealD E = real(sum(R))/ RealD(U.Grid()->gSites());
    RealD T = real( sum(qfield) );
    Coordinate scoor; for (int mu=0; mu < Nd; mu++) scoor[mu] = 0;
    RealD E0 = real(peekSite(R,scoor));
    RealD T0 = real(peekSite(qfield,scoor));
    std::cout << GridLogMessage << "[WilsonFlow] Saved energy density (clover) & topo. charge density: "  << conf << " " << step << "  " << tau << "  "
	      << "(E_avg,T_sum) " << E << " " << T << " (E, T at origin) " << E0 << " " << T0 << std::endl;
    
  });
  
  int t=WFPar.maxTau;
  WF.smear(Uflow, Umu);

  RealD WFlow_plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(Uflow);
  RealD WFlow_TC   = WilsonLoops<PeriodicGimplR>::TopologicalCharge(Uflow);
  RealD WFlow_T0   = WF.energyDensityPlaquette(t,Uflow); // t
  RealD WFlow_EC   = WF.energyDensityCloverleaf(t,Uflow);
  std::cout << GridLogMessage << "Plaquette          "<< conf << "   " << WFlow_plaq << std::endl;
  std::cout << GridLogMessage << "T0                 "<< conf << "   " << WFlow_T0 << std::endl;
  std::cout << GridLogMessage << "TC0                 "<< conf << "   " << WFlow_EC << std::endl;
  std::cout << GridLogMessage << "TopologicalCharge  "<< conf << "   " << WFlow_TC   << std::endl;

  std::cout<< GridLogMessage << " Admissibility check:\n";
  const double sp_adm = 0.067;                // admissible threshold
  const double pl_adm = 1.0-sp_adm/Nc;
  std::cout << GridLogMessage << "   (pl_adm =" << pl_adm << ")\n";

  // Need min and reduce min for this function
  //double sp_max = NC_*(1.0-stpl.plaq_min(U,pl_adm));
  double sp_ave = Nc*(1.0-WFlow_plaq);

  //std::cout<< GridLogMessage << "   sp_max = "        << sp_max <<"\n";
  std::cout<< GridLogMessage << "   sp_ave = "        << sp_ave <<"\n";
  std::cout<< GridLogMessage << "   (sp_admissible = "<< sp_adm <<")\n";
  //std::cout<< GridLogMessage << "   sp_admissible - sp_max = "<<sp_adm-sp_max <<"\n";
  std::cout<< GridLogMessage << "   sp_admissible - sp_ave = "<<sp_adm-sp_ave <<"\n";
  }
  Grid_finalize();
}  // main


/*
Input file example


JSON

{
    "WilsonFlow":{
	"steps": 200,
	"step_size": 0.01,
	"meas_interval": 50,
  "maxTau": 2.0
    },
    "Configurations":{
	"conf_prefix": "ckpoint_lat",
	"rng_prefix": "ckpoint_rng",
	"StartConfiguration": 3000,
	"EndConfiguration": 3000,
	"Skip": 5
    }
}


*/
