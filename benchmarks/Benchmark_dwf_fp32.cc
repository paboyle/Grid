 /*************************************************************************************
    Grid physics library, www.github.com/paboyle/Grid
    Source file: ./benchmarks/Benchmark_dwf.cc
    Copyright (C) 2015

    Author: Peter Boyle <paboyle@ph.ed.ac.uk>
    Author: paboyle <paboyle@ph.ed.ac.uk>

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
    /*  END LEGAL */
#include <Grid/Grid.h>
#ifdef GRID_CUDA
#define CUDA_PROFILE
#endif

#ifdef CUDA_PROFILE
#include <cuda_profiler_api.h>
#endif

using namespace std;
using namespace Grid;

////////////////////////
/// Move to domains ////
////////////////////////

struct DomainDecomposition
{
  Coordinate Block;

  DomainDecomposition(const Coordinate &_Block): Block(_Block){ assert(Block.size()==Nd);};
  
  template<class Field>
  void ProjectDomain(Field &f,Integer domain)
  {
    GridBase *grid = f.Grid();
    int dims = grid->Nd();
    int isDWF= (dims==Nd+1);
    assert((dims==Nd)||(dims==Nd+1));

    Field   zz(grid);  zz = Zero();
    LatticeInteger coor(grid);
    LatticeInteger domaincoor(grid);
    LatticeInteger mask(grid); mask = Integer(1);
    LatticeInteger zi(grid);     zi = Integer(0);
    for(int d=0;d<Nd;d++){
      Integer B= Block[d];
      if ( B ) {
	LatticeCoordinate(coor,d+isDWF);
	domaincoor = mod(coor,B);
	mask = where(domaincoor==Integer(0),zi,mask);
	mask = where(domaincoor==Integer(B-1),zi,mask);
      }
    }
    if ( !domain )
      f = where(mask==Integer(1),f,zz);
    else 
      f = where(mask==Integer(0),f,zz);
  };
};

template<typename MomentaField>
struct DirichletFilter: public MomentumFilterBase<MomentaField>
{
  Coordinate Block;
  
  DirichletFilter(const Coordinate &_Block): Block(_Block) {}

  // Edge detect using domain projectors
  void applyFilter (MomentaField &U) const override
  {
    DomainDecomposition Domains(Block);
    GridBase *grid = U.Grid();
    LatticeInteger  coor(grid);
    LatticeInteger  face(grid);
    LatticeInteger  one(grid);   one = 1;
    LatticeInteger  zero(grid); zero = 0;
    LatticeInteger  omega(grid);
    LatticeInteger  omegabar(grid);
    LatticeInteger  tmp(grid);

    omega=one;    Domains.ProjectDomain(omega,0);
    omegabar=one; Domains.ProjectDomain(omegabar,1);
    
    LatticeInteger nface(grid); nface=Zero();
    
    MomentaField projected(grid); projected=Zero();
    typedef decltype(PeekIndex<LorentzIndex>(U,0)) MomentaLinkField;
    MomentaLinkField  Umu(grid);
    MomentaLinkField   zz(grid); zz=Zero();

    int dims = grid->Nd();
    Coordinate Global=grid->GlobalDimensions();
    assert(dims==Nd);

    for(int mu=0;mu<Nd;mu++){

      if ( Block[mu]!=0 ) {

	Umu = PeekIndex<LorentzIndex>(U,mu);

	// Upper face 
 	tmp = Cshift(omegabar,mu,1);
	tmp = tmp + omega;
	face = where(tmp == Integer(2),one,zero );

 	tmp = Cshift(omega,mu,1);
	tmp = tmp + omegabar;
	face = where(tmp == Integer(2),one,face );

	Umu = where(face,zz,Umu);

	PokeIndex<LorentzIndex>(U, Umu, mu);
      }
    }
  }
};



Gamma::Algebra Gmu [] = {
			 Gamma::Algebra::GammaX,
			 Gamma::Algebra::GammaY,
			 Gamma::Algebra::GammaZ,
			 Gamma::Algebra::GammaT
};

void Benchmark(int Ls, std::vector<int> Dirichlet);

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);


  int threads = GridThread::GetThreads();

  int Ls=16;
  for(int i=0;i<argc;i++) {
    if(std::string(argv[i]) == "-Ls"){
      std::stringstream ss(argv[i+1]); ss >> Ls;
    }
  }
  std::vector<int> Dirichlet(5,0);
  Benchmark(Ls,Dirichlet);
  Coordinate latt4  = GridDefaultLatt();
  Coordinate mpi    = GridDefaultMpi();
  Coordinate shm;
  GlobalSharedMemory::GetShmDims(mpi,shm);
  /*
  Dirichlet = std::vector<int>({0,
				latt4[0]/mpi[0] * shm[0],
				latt4[1]/mpi[1] * shm[1],
				latt4[2]/mpi[2] * shm[2],
				latt4[3]/mpi[3] * shm[3]});
  */
  Dirichlet = std::vector<int>({0,
				latt4[0]/mpi[0] ,
				latt4[1]/mpi[1] ,
				latt4[2]/mpi[2] ,
				latt4[3]/mpi[3] });
  
  std::cout << " Dirichlet block "<< Dirichlet<< std::endl;
  Benchmark(Ls,Dirichlet);
  Grid_finalize();
  exit(0);
}
void Benchmark(int Ls, std::vector<int> Dirichlet)
{
  Coordinate latt4 = GridDefaultLatt();
  GridLogLayout();

  long unsigned int single_site_flops = 8*Nc*(7+16*Nc);

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  GridCartesian         * sUGrid   = SpaceTimeGrid::makeFourDimDWFGrid(GridDefaultLatt(),GridDefaultMpi());
  GridRedBlackCartesian * sUrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(sUGrid);
  GridCartesian         * sFGrid   = SpaceTimeGrid::makeFiveDimDWFGrid(Ls,UGrid);
  GridRedBlackCartesian * sFrbGrid = SpaceTimeGrid::makeFiveDimDWFRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});

  std::cout << GridLogMessage << "Initialising 4d RNG" << std::endl;
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedUniqueString(std::string("The 4D RNG"));

  std::cout << GridLogMessage << "Initialising 5d RNG" << std::endl;
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedUniqueString(std::string("The 5D RNG"));

  LatticeFermionF src   (FGrid); random(RNG5,src);
  RealD N2 = 1.0/::sqrt(norm2(src));
  src = src*N2;

  LatticeFermionF result(FGrid); result=Zero();
  LatticeFermionF    ref(FGrid);    ref=Zero();
  LatticeFermionF    tmp(FGrid);
  LatticeFermionF    err(FGrid);

  std::cout << GridLogMessage << "Drawing gauge field" << std::endl;
  LatticeGaugeFieldF Umu(UGrid);
  SU<Nc>::HotConfiguration(RNG4,Umu);
  std::cout << GridLogMessage << "Random gauge initialised " << std::endl;

  ////////////////////////////////////
  // Apply BCs
  ////////////////////////////////////
  std::cout << GridLogMessage << "Applying BCs " << std::endl;
  Coordinate Block(4);
  for(int d=0;d<4;d++)  Block[d]= Dirichlet[d+1];

  std::cout << GridLogMessage << "Dirichlet Block " << Block<< std::endl;
  DirichletFilter<LatticeGaugeFieldF> Filter(Block);
  Filter.applyFilter(Umu);
  
  ////////////////////////////////////
  // Naive wilson implementation
  ////////////////////////////////////
  // replicate across fifth dimension
  //  LatticeGaugeFieldF Umu5d(FGrid);
  std::vector<LatticeColourMatrixF> U(4,UGrid);
  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }
  std::cout << GridLogMessage << "Setting up Cshift based reference " << std::endl;

  if (1)
  {
    ref = Zero();
    for(int mu=0;mu<Nd;mu++){

      tmp = Cshift(src,mu+1,1);
      {
	autoView( tmp_v  , tmp  , CpuWrite);
	autoView( U_v  , U[mu]  , CpuRead);
	for(int ss=0;ss<U[mu].Grid()->oSites();ss++){
	  for(int s=0;s<Ls;s++){
	    tmp_v[Ls*ss+s] = U_v[ss]*tmp_v[Ls*ss+s];
	  }
	}
      }
      ref=ref + tmp - Gamma(Gmu[mu])*tmp;

      {
	autoView( tmp_v  , tmp  , CpuWrite);
	autoView( U_v  , U[mu]  , CpuRead);
	autoView( src_v, src    , CpuRead);
	for(int ss=0;ss<U[mu].Grid()->oSites();ss++){
	  for(int s=0;s<Ls;s++){
	    tmp_v[Ls*ss+s] = adj(U_v[ss])*src_v[Ls*ss+s];
	  }
	}
      }
      tmp =Cshift(tmp,mu+1,-1);
      ref=ref + tmp + Gamma(Gmu[mu])*tmp;
    }
    ref = -0.5*ref;
  }

  RealD mass=0.1;
  RealD M5  =1.8;

  RealD NP = UGrid->_Nprocessors;
  RealD NN = UGrid->NodeCount();

  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Kernel options --dslash-generic, --dslash-unroll, --dslash-asm" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Benchmarking DomainWallFermionR::Dhop                  "<<std::endl;
  std::cout << GridLogMessage<< "* Vectorising space-time by "<<vComplexF::Nsimd()<<std::endl;
  std::cout << GridLogMessage<< "* VComplexF size is "<<sizeof(vComplexF)<< " B"<<std::endl;
  if ( sizeof(RealF)==4 )   std::cout << GridLogMessage<< "* SINGLE precision "<<std::endl;
  if ( sizeof(RealF)==8 )   std::cout << GridLogMessage<< "* DOUBLE precision "<<std::endl;
#ifdef GRID_OMP
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsAndCompute ) std::cout << GridLogMessage<< "* Using Overlapped Comms/Compute" <<std::endl;
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsThenCompute) std::cout << GridLogMessage<< "* Using sequential comms compute" <<std::endl;
#endif
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptHandUnroll) std::cout << GridLogMessage<< "* Using Nc=3       WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptInlineAsm ) std::cout << GridLogMessage<< "* Using Asm Nc=3   WilsonKernels" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;

  DomainWallFermionF Dw(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  Dw.DirichletBlock(Dirichlet);
  int ncall =300;

  if (1) {
    FGrid->Barrier();
    Dw.Dhop(src,result,0);
    std::cout<<GridLogMessage<<"Called warmup"<<std::endl;
    double t0=usecond();
    for(int i=0;i<ncall;i++){
      __SSC_START;
      Dw.Dhop(src,result,0);
      __SSC_STOP;
    }
    double t1=usecond();
    FGrid->Barrier();

    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=single_site_flops*volume*ncall;

    auto nsimd = vComplex::Nsimd();
    auto simdwidth = sizeof(vComplex);

    // RF: Nd Wilson * Ls, Nd gauge * Ls, Nc colors
    double data_rf = volume * ((2*Nd+1)*Nd*Nc + 2*Nd*Nc*Nc) * simdwidth / nsimd * ncall / (1024.*1024.*1024.);

    // mem: Nd Wilson * Ls, Nd gauge, Nc colors
    double data_mem = (volume * (2*Nd+1)*Nd*Nc + (volume/Ls) *2*Nd*Nc*Nc) * simdwidth / nsimd * ncall / (1024.*1024.*1024.);

    std::cout<<GridLogMessage << "Called Dw "<<ncall<<" times in "<<t1-t0<<" us"<<std::endl;
    std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t1-t0)<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per rank =  "<< flops/(t1-t0)/NP<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per node =  "<< flops/(t1-t0)/NN<<std::endl;
    std::cout<<GridLogMessage << "RF  GiB/s (base 2) =   "<< 1000000. * data_rf/((t1-t0))<<std::endl;
    std::cout<<GridLogMessage << "mem GiB/s (base 2) =   "<< 1000000. * data_mem/((t1-t0))<<std::endl;
    err = ref-result;
    std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<<std::endl;

    if(( norm2(err)>1.0e-4) ) {
      std::cout<<GridLogMessage << "WRONG RESULT" << std::endl;
      FGrid->Barrier();
      exit(-1);
    }
    assert (norm2(err)< 1.0e-4 );
  }

  if (1)
  { // Naive wilson dag implementation
    ref = Zero();
    for(int mu=0;mu<Nd;mu++){

      //    ref =  src - Gamma(Gamma::Algebra::GammaX)* src ; // 1+gamma_x
      tmp = Cshift(src,mu+1,1);
      {
	autoView( ref_v, ref, CpuWrite);
	autoView( tmp_v, tmp, CpuRead);
	autoView( U_v  , U[mu]  , CpuRead);
	for(int ss=0;ss<U[mu].Grid()->oSites();ss++){
	  for(int s=0;s<Ls;s++){
	    int i=s+Ls*ss;
	    ref_v[i]+= U_v[ss]*(tmp_v[i] + Gamma(Gmu[mu])*tmp_v[i]); ;
	  }
	}
      }
      
      {
	autoView( tmp_v  , tmp  , CpuWrite);
	autoView( U_v  , U[mu]  , CpuRead);
	autoView( src_v, src    , CpuRead);
	for(int ss=0;ss<U[mu].Grid()->oSites();ss++){
	  for(int s=0;s<Ls;s++){
	    tmp_v[Ls*ss+s] = adj(U_v[ss])*src_v[Ls*ss+s];
	  }
	}
      }
      //      tmp =adj(U[mu])*src;
      tmp =Cshift(tmp,mu+1,-1);
      {
	autoView( ref_v, ref, CpuWrite);
	autoView( tmp_v, tmp, CpuRead);
	for(int i=0;i<ref_v.size();i++){
	  ref_v[i]+= tmp_v[i] - Gamma(Gmu[mu])*tmp_v[i]; ;
	}
      }
    }
    ref = -0.5*ref;
  }
  //  dump=1;
  Dw.Dhop(src,result,1);
  std::cout << GridLogMessage << "Compare to naive wilson implementation Dag to verify correctness" << std::endl;
  std::cout<<GridLogMessage << "Called DwDag"<<std::endl;
  std::cout<<GridLogMessage << "norm dag result "<< norm2(result)<<std::endl;
  std::cout<<GridLogMessage << "norm dag ref    "<< norm2(ref)<<std::endl;
  err = ref-result;
  std::cout<<GridLogMessage << "norm dag diff   "<< norm2(err)<<std::endl;

  if (  norm2(err) > 1.0e-4 ) {
    std::cout << "Error vector is\n" <<err << std::endl;
    std::cout << "Ref   vector is\n" <<ref << std::endl;
    std::cout << "Result  vector is\n" <<result << std::endl;
  }
  assert((norm2(err)<1.0e-4));
  
  LatticeFermionF src_e (FrbGrid);
  LatticeFermionF src_o (FrbGrid);
  LatticeFermionF r_e   (FrbGrid);
  LatticeFermionF r_o   (FrbGrid);
  LatticeFermionF r_eo  (FGrid);

  std::cout<<GridLogMessage << "Calling Deo and Doe and //assert Deo+Doe == Dunprec"<<std::endl;
  pickCheckerboard(Even,src_e,src);
  pickCheckerboard(Odd,src_o,src);

  std::cout<<GridLogMessage << "src_e"<<norm2(src_e)<<std::endl;
  std::cout<<GridLogMessage << "src_o"<<norm2(src_o)<<std::endl;


  // S-direction is INNERMOST and takes no part in the parity.
  std::cout << GridLogMessage<< "*********************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Benchmarking DomainWallFermionF::DhopEO                "<<std::endl;
  std::cout << GridLogMessage<< "* Vectorising space-time by "<<vComplexF::Nsimd()<<std::endl;
  if ( sizeof(RealF)==4 )   std::cout << GridLogMessage<< "* SINGLE precision "<<std::endl;
  if ( sizeof(RealF)==8 )   std::cout << GridLogMessage<< "* DOUBLE precision "<<std::endl;
#ifdef GRID_OMP
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsAndCompute ) std::cout << GridLogMessage<< "* Using Overlapped Comms/Compute" <<std::endl;
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsThenCompute) std::cout << GridLogMessage<< "* Using sequential comms compute" <<std::endl;
#endif
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptHandUnroll) std::cout << GridLogMessage<< "* Using Nc=3       WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptInlineAsm ) std::cout << GridLogMessage<< "* Using Asm Nc=3   WilsonKernels" <<std::endl;
  std::cout << GridLogMessage<< "*********************************************************" <<std::endl;
  {
    FGrid->Barrier();
    Dw.DhopEO(src_o,r_e,DaggerNo);
    double t0=usecond();
    for(int i=0;i<ncall;i++){
#ifdef CUDA_PROFILE
      if(i==10) cudaProfilerStart();
#endif
      Dw.DhopEO(src_o,r_e,DaggerNo);
#ifdef CUDA_PROFILE
      if(i==20) cudaProfilerStop();
#endif
    }
    double t1=usecond();
    FGrid->Barrier();

    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=(single_site_flops*volume*ncall)/2.0;

    std::cout<<GridLogMessage << "Deo mflop/s =   "<< flops/(t1-t0)<<std::endl;
    std::cout<<GridLogMessage << "Deo mflop/s per rank   "<< flops/(t1-t0)/NP<<std::endl;
    std::cout<<GridLogMessage << "Deo mflop/s per node   "<< flops/(t1-t0)/NN<<std::endl;
  }
  Dw.DhopEO(src_o,r_e,DaggerNo);
  Dw.DhopOE(src_e,r_o,DaggerNo);
  Dw.Dhop  (src  ,result,DaggerNo);

  std::cout<<GridLogMessage << "r_e"<<norm2(r_e)<<std::endl;
  std::cout<<GridLogMessage << "r_o"<<norm2(r_o)<<std::endl;
  std::cout<<GridLogMessage << "res"<<norm2(result)<<std::endl;

  setCheckerboard(r_eo,r_o);
  setCheckerboard(r_eo,r_e);

  err = r_eo-result;
  std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<<std::endl;
  assert(norm2(err)<1.0e-4);

  pickCheckerboard(Even,src_e,err);
  pickCheckerboard(Odd,src_o,err);
  std::cout<<GridLogMessage << "norm diff even  "<< norm2(src_e)<<std::endl;
  std::cout<<GridLogMessage << "norm diff odd   "<< norm2(src_o)<<std::endl;

  assert(norm2(src_e)<1.0e-4);
  assert(norm2(src_o)<1.0e-4);
}
