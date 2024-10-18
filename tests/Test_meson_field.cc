/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: tests/core/Test_meson_field.cc

Copyright (C) 2015-2018

Author: Felix Erben <felix.erben@ed.ac.uk>

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

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/

#include <Grid/Grid.h>
#include <Grid/qcd/utils/A2Autils.h>

using namespace Grid;

const int TSRC = 0;  //timeslice where rho is nonzero
const int VDIM = 5; //length of each vector

typedef typename DomainWallFermionD::ComplexField ComplexField;
typedef typename DomainWallFermionD::FermionField FermionField;

int main(int argc, char *argv[])
{
  // initialization
  Grid_init(&argc, &argv);
  std::cout << GridLogMessage << "Grid initialized" << std::endl;

  // Lattice and rng setup 
  Coordinate latt_size   = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(4, vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();
  GridCartesian    grid(latt_size,simd_layout,mpi_layout);
  int Nt = GridDefaultLatt()[Tp];
  Lattice<iScalar<vInteger>> t(&grid);
  LatticeCoordinate(t, Tp);
  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG(&grid);
  pRNG.SeedFixedIntegers(seeds);

  // MesonField lhs and rhs vectors
  std::vector<FermionField> phi(VDIM,&grid);
  std::cout << GridLogMessage << "Initialising random meson fields" << std::endl;
  for (unsigned int i = 0; i < VDIM; ++i){
    random(pRNG,phi[i]);
  }
  std::cout << GridLogMessage << "Meson fields initialised, rho non-zero only for t = " << TSRC << std::endl;

  // Gamma matrices used in the contraction
  std::vector<Gamma::Algebra> Gmu = {
    Gamma::Algebra::GammaX,
    Gamma::Algebra::GammaY,
    Gamma::Algebra::GammaZ,
    Gamma::Algebra::GammaT
  };

  // momentum phases e^{ipx}
  std::vector<std::vector<double>> momenta = {
	  {0.,0.,0.},
	  {1.,0.,0.},
	  {1.,1.,0.},
	  {1.,1.,1.},
	  {2.,0.,0.}
  };
  // 5 momenta x VDIMxVDIM = 125 calls (x 16 spins) 1.4s => 1400/125 ~10ms per call
  std::cout << GridLogMessage << "Meson fields will be created for " << Gmu.size() << " Gamma matrices and " << momenta.size() << " momenta." << std::endl;

  std::cout << GridLogMessage << "Computing complex phases" << std::endl;
  std::vector<ComplexField> phases(momenta.size(),&grid);
  ComplexField coor(&grid);
  Complex Ci(0.0,1.0);
  for (unsigned int j = 0; j < momenta.size(); ++j)
  {
    phases[j] = Zero();
    for(unsigned int mu = 0; mu < momenta[j].size(); mu++)
    {
      LatticeCoordinate(coor, mu);
      phases[j] = phases[j] + momenta[j][mu]/GridDefaultLatt()[mu]*coor;
    }
    phases[j] = exp((Real)(2*M_PI)*Ci*phases[j]);
  }
  std::cout << GridLogMessage << "Computing complex phases done." << std::endl;

  Eigen::Tensor<ComplexD,5, Eigen::RowMajor> Mpp(momenta.size(),Gmu.size(),Nt,VDIM,VDIM);
  Eigen::Tensor<ComplexD,5, Eigen::RowMajor> Mpp_gpu(momenta.size(),Gmu.size(),Nt,VDIM,VDIM);

  // timer
  double start,stop;

  //execute meson field routine
  std::cout << GridLogMessage << "Meson Field Warmup Begin" << std::endl;
  A2Autils<WilsonImplR>::MesonField(Mpp,&phi[0],&phi[0],Gmu,phases,Tp);
  std::cout << GridLogMessage << "Meson Field Timing Begin" << std::endl;
  start = usecond();
  A2Autils<WilsonImplR>::MesonField(Mpp,&phi[0],&phi[0],Gmu,phases,Tp);
  stop = usecond();
  std::cout << GridLogMessage << "M(phi,phi) created, execution time " << stop-start << " us" << std::endl;

  std::cout << GridLogMessage << "Meson Field GPU Warmup Begin" << std::endl;
  A2Autils<WilsonImplR>::MesonFieldGPU(Mpp_gpu,&phi[0],&phi[0],Gmu,phases,Tp);
  std::cout << GridLogMessage << "Meson Field GPU Timing Begin" << std::endl;
  start = usecond();
  A2Autils<WilsonImplR>::MesonFieldGPU(Mpp_gpu,&phi[0],&phi[0],Gmu,phases,Tp);
  stop = usecond();
  std::cout << GridLogMessage << "M_gpu(phi,phi) created, execution time " << stop-start << " us" << std::endl;

  for(int mom=0;mom<momenta.size();mom++){
    for(int mu=0;mu<Gmu.size();mu++){
      for(int t=0;t<Nt;t++){
	for(int v=0;v<VDIM;v++){
	  for(int w=0;w<VDIM;w++){
	    std::cout << GridLogMessage
		      << " " << mom
		      << " " << mu
		      << " " << t
		      << " " << v
		      << " " << w
		      << " " << Mpp_gpu(mom,mu,t,v,w)
		      << " " << Mpp(mom,mu,t,v,w) << std::endl;
	  }
	}
      }
    }
  }
  
  std::string FileName = "Meson_Fields";
#ifdef HAVE_HDF5
  using Default_Reader = Grid::Hdf5Reader;
  using Default_Writer = Grid::Hdf5Writer;
  FileName.append(".h5");
#else
  using Default_Reader = Grid::BinaryReader;
  using Default_Writer = Grid::BinaryWriter;
  FileName.append(".bin");
#endif
  {
    Default_Writer w(FileName);
    write(w,"phi_phi",Mpp);
    write(w,"phi_phi_gpu",Mpp_gpu);
  }
  // epilogue
  std::cout << GridLogMessage << "Grid is finalizing now" << std::endl;
  Grid_finalize();
  
  return EXIT_SUCCESS;
}
