    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./benchmarks/Benchmark_wilson.cc

    Copyright (C) 2018

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

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


#include "Grid/util/Profiling.h"

template<class vobj>
void sliceInnerProductMesonField(std::vector< std::vector<ComplexD> > &mat, 
				 const std::vector<Lattice<vobj> > &lhs,
				 const std::vector<Lattice<vobj> > &rhs,
				 int orthogdim) 
{
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  
  int Lblock = lhs.size();
  int Rblock = rhs.size();

  GridBase *grid = lhs[0]._grid;
  
  const int    Nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();
  int Nt     = grid->GlobalDimensions()[orthogdim];

  assert(mat.size()==Lblock*Rblock);
  for(int t=0;t<mat.size();t++){
    assert(mat[t].size()==Nt);
  }

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  // will locally sum vectors first
  // sum across these down to scalars
  // splitting the SIMD
  std::vector<vector_type,alignedAllocator<vector_type> > lvSum(rd*Lblock*Rblock);
  parallel_for (int r = 0; r < rd * Lblock * Rblock; r++){
    lvSum[r] = zero;
  }

  std::vector<scalar_type > lsSum(ld*Lblock*Rblock,scalar_type(0.0));             

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];
  
  std::cout << GridLogMessage << " Entering first parallel loop "<<std::endl;
  // Parallelise over t-direction doesn't expose as much parallelism as needed for KNL
  parallel_for(int r=0;r<rd;r++){

    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int ss= so+n*stride+b;
	for(int i=0;i<Lblock;i++){
	  auto left = conjugate(lhs[i]._odata[ss]);
	  for(int j=0;j<Rblock;j++){
	    int idx = i+Lblock*j+Lblock*Rblock*r;
	    auto right = rhs[j]._odata[ss];
	    vector_type vv = left()(0)(0) * right()(0)(0)
	      +              left()(0)(1) * right()(0)(1)
	      +              left()(0)(2) * right()(0)(2)
              +              left()(1)(0) * right()(1)(0)
	      +              left()(1)(1) * right()(1)(1)
	      +              left()(1)(2) * right()(1)(2)
              +              left()(2)(0) * right()(2)(0)
	      +              left()(2)(1) * right()(2)(1)
	      +              left()(2)(2) * right()(2)(2)
              +              left()(3)(0) * right()(3)(0)
	      +              left()(3)(1) * right()(3)(1)
	      +              left()(3)(2) * right()(3)(2);
	    lvSum[idx]=lvSum[idx]+vv;
	  }
	}
      }
    }
  }

  std::cout << GridLogMessage << " Entering second parallel loop "<<std::endl;
  // Sum across simd lanes in the plane, breaking out orthog dir.
  parallel_for(int rt=0;rt<rd;rt++){

    std::vector<int> icoor(Nd);

    for(int i=0;i<Lblock;i++){
    for(int j=0;j<Rblock;j++){

      iScalar<vector_type> temp; 
      std::vector<iScalar<scalar_type> > extracted(Nsimd);               

      temp._internal = lvSum[i+Lblock*j+Lblock*Rblock*rt];

      extract(temp,extracted);

      for(int idx=0;idx<Nsimd;idx++){

	grid->iCoorFromIindex(icoor,idx);

	int ldx =rt+icoor[orthogdim]*rd;
      
	int ij_dx = i+Lblock*j+Lblock*Rblock*ldx;
	lsSum[ij_dx]=lsSum[ij_dx]+extracted[idx]._internal;

      }
    }}
  }

  std::cout << GridLogMessage << " Entering non parallel loop "<<std::endl;
  for(int t=0;t<fd;t++)
  {
    int pt = t / ld; // processor plane
    int lt = t % ld;
    for(int i=0;i<Lblock;i++){
    for(int j=0;j<Rblock;j++){
      if (pt == grid->_processor_coor[orthogdim]){
        int ij_dx = i + Lblock * j + Lblock * Rblock * lt;
        mat[i+j*Lblock][t] = lsSum[ij_dx];
      }
      else{
        mat[i+j*Lblock][t] = scalar_type(0.0);
      }
    }}
  }
  std::cout << GridLogMessage << " Done "<<std::endl;
  // defer sum over nodes.
  return;
}

template<class vobj>
void sliceInnerProductMesonFieldGamma(std::vector< std::vector<ComplexD> > &mat, 
				      const std::vector<Lattice<vobj> > &lhs,
				      const std::vector<Lattice<vobj> > &rhs,
				      int orthogdim,
				      std::vector<Gamma::Algebra> gammas) 
{
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  
  int Lblock = lhs.size();
  int Rblock = rhs.size();

  GridBase *grid = lhs[0]._grid;
  
  const int    Nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();
  int Nt     = grid->GlobalDimensions()[orthogdim];
  int Ngamma = gammas.size();

  assert(mat.size()==Lblock*Rblock*Ngamma);
  for(int t=0;t<mat.size();t++){
    assert(mat[t].size()==Nt);
  }

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  // will locally sum vectors first
  // sum across these down to scalars
  // splitting the SIMD
  int MFrvol = rd*Lblock*Rblock*Ngamma;
  int MFlvol = ld*Lblock*Rblock*Ngamma;

  std::vector<vector_type,alignedAllocator<vector_type> > lvSum(MFrvol);
  parallel_for (int r = 0; r < MFrvol; r++){
    lvSum[r] = zero;
  }

  std::vector<scalar_type > lsSum(MFlvol);             
  parallel_for (int r = 0; r < MFlvol; r++){
    lsSum[r]=scalar_type(0.0);
  }

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];
  
  std::cout << GridLogMessage << " Entering first parallel loop "<<std::endl;

  // Parallelise over t-direction doesn't expose as much parallelism as needed for KNL
  parallel_for(int r=0;r<rd;r++){

    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int ss= so+n*stride+b;
	for(int i=0;i<Lblock;i++){
	    auto left = conjugate(lhs[i]._odata[ss]);
	    for(int j=0;j<Rblock;j++){
	  for(int mu=0;mu<Ngamma;mu++){

        auto right = Gamma(gammas[mu])*rhs[j]._odata[ss];

	      vector_type vv = left()(0)(0) * right()(0)(0)
		+              left()(0)(1) * right()(0)(1)
		+              left()(0)(2) * right()(0)(2)
		+              left()(1)(0) * right()(1)(0)
		+              left()(1)(1) * right()(1)(1)
		+              left()(1)(2) * right()(1)(2)
		+              left()(2)(0) * right()(2)(0)
		+              left()(2)(1) * right()(2)(1)
		+              left()(2)(2) * right()(2)(2)
		+              left()(3)(0) * right()(3)(0)
		+              left()(3)(1) * right()(3)(1)
		+              left()(3)(2) * right()(3)(2);

	      int idx = mu+i*Ngamma+Lblock*Ngamma*j+Ngamma*Lblock*Rblock*r;

	      lvSum[idx]=lvSum[idx]+vv;
	    }
	  }
	}
      }
    }
  }

  std::cout << GridLogMessage << " Entering second parallel loop "<<std::endl;
  // Sum across simd lanes in the plane, breaking out orthog dir.
  parallel_for(int rt=0;rt<rd;rt++){

    iScalar<vector_type> temp; 
    std::vector<int> icoor(Nd);
    std::vector<iScalar<scalar_type> > extracted(Nsimd);               

    for(int i=0;i<Lblock;i++){
    for(int j=0;j<Rblock;j++){
    for(int mu=0;mu<Ngamma;mu++){

      int ij_rdx = mu+i*Ngamma+Ngamma*Lblock*j+Ngamma*Lblock*Rblock*rt;
      temp._internal = lvSum[ij_rdx];

      extract(temp,extracted);

      for(int idx=0;idx<Nsimd;idx++){

	grid->iCoorFromIindex(icoor,idx);

	int ldx =rt+icoor[orthogdim]*rd;
      
	int ij_ldx = mu+i*Ngamma+Ngamma*Lblock*j+Ngamma*Lblock*Rblock*ldx;
	lsSum[ij_ldx]=lsSum[ij_ldx]+extracted[idx]._internal;

      }
    }}}
  }

  std::cout << GridLogMessage << " Entering non parallel loop "<<std::endl;
  for(int t=0;t<fd;t++)
  {
    int pt = t / ld; // processor plane
    int lt = t % ld;
    for(int i=0;i<Lblock;i++){
    for(int j=0;j<Rblock;j++){
    for(int mu=0;mu<Ngamma;mu++){
      if (pt == grid->_processor_coor[orthogdim]){
        int ij_dx = mu+i*Ngamma+Ngamma*Lblock*j+Ngamma*Lblock*Rblock* lt;
        mat[mu+i*Ngamma+j*Lblock*Ngamma][t] = lsSum[ij_dx];
      }
      else{
        mat[mu+i*Ngamma+j*Lblock*Ngamma][t] = scalar_type(0.0);
      }
    }}}
  }
  std::cout << GridLogMessage << " Done "<<std::endl;
  // defer sum over nodes.
  return;
}


template<class vobj>
void sliceInnerProductMesonFieldGamma1(std::vector< std::vector<ComplexD> > &mat, 
				      const std::vector<Lattice<vobj> > &lhs,
				      const std::vector<Lattice<vobj> > &rhs,
				      int orthogdim,
				      std::vector<Gamma::Algebra> gammas) 
{

  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iSpinMatrix<scalar_type> SpinMatrix_s;
  
  int Lblock = lhs.size();
  int Rblock = rhs.size();

  GridBase *grid = lhs[0]._grid;
  
  const int    Nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();
  int Nt     = grid->GlobalDimensions()[orthogdim];
  int Ngamma = gammas.size();

  assert(mat.size()==Lblock*Rblock*Ngamma);
  for(int t=0;t<mat.size();t++){
    assert(mat[t].size()==Nt);
  }

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  // will locally sum vectors first
  // sum across these down to scalars
  // splitting the SIMD
  int MFrvol = rd*Lblock*Rblock;
  int MFlvol = ld*Lblock*Rblock;

  Vector<SpinMatrix_v > lvSum(MFrvol);
  parallel_for (int r = 0; r < MFrvol; r++){
    lvSum[r] = zero;
  }

  Vector<SpinMatrix_s > lsSum(MFlvol);             
  parallel_for (int r = 0; r < MFlvol; r++){
    lsSum[r]=scalar_type(0.0);
  }

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];
  
  std::cout << GridLogMessage << " Entering first parallel loop "<<std::endl;

  // Parallelise over t-direction doesn't expose as much parallelism as needed for KNL
  parallel_for(int r=0;r<rd;r++){

    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int ss= so+n*stride+b;
	for(int i=0;i<Lblock;i++){

	  auto left = conjugate(lhs[i]._odata[ss]);
	  for(int j=0;j<Rblock;j++){

	    SpinMatrix_v vv;
	    auto right = rhs[j]._odata[ss];
	    for(int s1=0;s1<Ns;s1++){
	    for(int s2=0;s2<Ns;s2++){
	     vv()(s2,s1)() = left()(s1)(0) * right()(s2)(0)
		+             left()(s1)(1) * right()(s2)(1)
		+             left()(s1)(2) * right()(s2)(2);
	    }}

	    int idx = i+Lblock*j+Lblock*Rblock*r;

	    lvSum[idx]=lvSum[idx]+vv;
	  
	  }
	  }
	}
      }
    }

  std::cout << GridLogMessage << " Entering second parallel loop "<<std::endl;
  // Sum across simd lanes in the plane, breaking out orthog dir.
  parallel_for(int rt=0;rt<rd;rt++){

    std::vector<int> icoor(Nd);
    std::vector<SpinMatrix_s> extracted(Nsimd);               

    for(int i=0;i<Lblock;i++){
    for(int j=0;j<Rblock;j++){

      int ij_rdx = i+Lblock*j+Lblock*Rblock*rt;

      extract(lvSum[ij_rdx],extracted);

      for(int idx=0;idx<Nsimd;idx++){

	grid->iCoorFromIindex(icoor,idx);

	int ldx    = rt+icoor[orthogdim]*rd;

	int ij_ldx = i+Lblock*j+Lblock*Rblock*ldx;

	lsSum[ij_ldx]=lsSum[ij_ldx]+extracted[idx];

      }
    }}
  }

  std::cout << GridLogMessage << " Entering third parallel loop "<<std::endl;
  parallel_for(int t=0;t<fd;t++)
  {
    int pt = t / ld; // processor plane
    int lt = t % ld;
    for(int i=0;i<Lblock;i++){
    for(int j=0;j<Rblock;j++){
      if (pt == grid->_processor_coor[orthogdim]){
        int ij_dx = i + Lblock * j + Lblock * Rblock * lt;
    	for(int mu=0;mu<Ngamma;mu++){
	  mat[mu+i*Ngamma+j*Lblock*Ngamma][t] = trace(lsSum[ij_dx]*Gamma(gammas[mu]));
	}
      }
      else{
        for(int mu=0;mu<Ngamma;mu++){
	  mat[mu+i*Ngamma+j*Lblock*Ngamma][t] = scalar_type(0.0);
	}
      }
    }}
  }
  std::cout << GridLogMessage << " Done "<<std::endl;
  // defer sum over nodes.
  return;
}

template<class vobj>
void sliceInnerProductMesonFieldGammaMom(std::vector< std::vector<ComplexD> > &mat, 
					 const std::vector<Lattice<vobj> > &lhs,
					 const std::vector<Lattice<vobj> > &rhs,
					 int orthogdim,
					 std::vector<Gamma::Algebra> gammas,
					 const std::vector<LatticeComplex > &mom) 
{
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iSpinMatrix<scalar_type> SpinMatrix_s;
  
  int Lblock = lhs.size();
  int Rblock = rhs.size();

  GridBase *grid = lhs[0]._grid;
  
  const int    Nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();
  int Nt     = grid->GlobalDimensions()[orthogdim];
  int Ngamma = gammas.size();
  int Nmom   = mom.size();

  assert(mat.size()==Lblock*Rblock*Ngamma*Nmom);
  for(int t=0;t<mat.size();t++){
    assert(mat[t].size()==Nt);
  }

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  // will locally sum vectors first
  // sum across these down to scalars
  // splitting the SIMD
  int MFrvol = rd*Lblock*Rblock*Nmom;
  int MFlvol = ld*Lblock*Rblock*Nmom;

  Vector<SpinMatrix_v > lvSum(MFrvol);
  parallel_for (int r = 0; r < MFrvol; r++){
    lvSum[r] = zero;
  }

  Vector<SpinMatrix_s > lsSum(MFlvol);             
  parallel_for (int r = 0; r < MFlvol; r++){
    lsSum[r]=scalar_type(0.0);
  }

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];
  
  std::cout << GridLogMessage << " Entering first parallel loop "<<std::endl;

  // Parallelise over t-direction doesn't expose as much parallelism as needed for KNL
  parallel_for(int r=0;r<rd;r++){

    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int ss= so+n*stride+b;


	for(int i=0;i<Lblock;i++){

	  auto left = conjugate(lhs[i]._odata[ss]);
	  for(int j=0;j<Rblock;j++){

	    SpinMatrix_v vv;
	    auto right = rhs[j]._odata[ss];
	    for(int s1=0;s1<Ns;s1++){
	    for(int s2=0;s2<Ns;s2++){
	      vv()(s1,s2)() = left()(s1)(0) * right()(s2)(0)
		+             left()(s1)(1) * right()(s2)(1)
		+             left()(s1)(2) * right()(s2)(2);
	    }}
	    
	    // After getting the sitewise product do the mom phase loop
	    int base = Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*r;
	    // Trigger unroll
	    for ( int m=0;m<Nmom;m++){
	      int idx = m+base;
	      auto phase = mom[m]._odata[ss];
	      mac(&lvSum[idx],&vv,&phase);
	    }
	  
	  }
	}
      }
    }
  }

  std::cout << GridLogMessage << " Entering second parallel loop "<<std::endl;
  // Sum across simd lanes in the plane, breaking out orthog dir.
  parallel_for(int rt=0;rt<rd;rt++){

    std::vector<int> icoor(Nd);
    std::vector<SpinMatrix_s> extracted(Nsimd);               


    for(int i=0;i<Lblock;i++){
    for(int j=0;j<Rblock;j++){
    for(int m=0;m<Nmom;m++){

      int ij_rdx = m+Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*rt;

      extract(lvSum[ij_rdx],extracted);

      for(int idx=0;idx<Nsimd;idx++){

	grid->iCoorFromIindex(icoor,idx);

	int ldx    = rt+icoor[orthogdim]*rd;

	int ij_ldx = m+Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*ldx;

	lsSum[ij_ldx]=lsSum[ij_ldx]+extracted[idx];

      }
    }}}
  }

  std::cout << GridLogMessage << " Entering third parallel loop "<<std::endl;
  parallel_for(int t=0;t<fd;t++)
  {
    int pt = t / ld; // processor plane
    int lt = t % ld;
    for(int i=0;i<Lblock;i++){
    for(int j=0;j<Rblock;j++){
      if (pt == grid->_processor_coor[orthogdim]){
	for(int m=0;m<Nmom;m++){
	  int ij_dx = m+Nmom*i + Nmom*Lblock * j + Nmom*Lblock * Rblock * lt;
	  for(int mu=0;mu<Ngamma;mu++){
	    mat[ mu
		+m*Ngamma
		+i*Nmom*Ngamma
		+j*Nmom*Ngamma*Lblock][t] = trace(lsSum[ij_dx]*Gamma(gammas[mu]));
	  }
	}
      }
      else{
	for(int mu=0;mu<Ngamma;mu++){
	for(int m=0;m<Nmom;m++){
	  mat[mu+m*Ngamma+i*Nmom*Ngamma+j*Nmom*Lblock*Ngamma][t] = scalar_type(0.0);
	}}
      }
    }}
  }
  std::cout << GridLogMessage << " Done "<<std::endl;
  // defer sum over nodes.
  return;
}



/*
template void sliceInnerProductMesonField<SpinColourVector>(std::vector< std::vector<ComplexD> > &mat, 
						   const std::vector<Lattice<SpinColourVector> > &lhs,
						   const std::vector<Lattice<SpinColourVector> > &rhs,
						   int orthogdim) ;
*/

std::vector<Gamma::Algebra> Gmu4 ( {
  Gamma::Algebra::GammaX,
  Gamma::Algebra::GammaY,
  Gamma::Algebra::GammaZ,
  Gamma::Algebra::GammaT });

std::vector<Gamma::Algebra> Gmu16 ( {
  Gamma::Algebra::Gamma5,
  Gamma::Algebra::GammaT,
  Gamma::Algebra::GammaTGamma5,
  Gamma::Algebra::GammaX,
  Gamma::Algebra::GammaXGamma5,
  Gamma::Algebra::GammaY,
  Gamma::Algebra::GammaYGamma5,
  Gamma::Algebra::GammaZ,
  Gamma::Algebra::GammaZGamma5,
  Gamma::Algebra::Identity,
  Gamma::Algebra::SigmaXT,
  Gamma::Algebra::SigmaXY,
  Gamma::Algebra::SigmaXZ,
  Gamma::Algebra::SigmaYT,
  Gamma::Algebra::SigmaYZ,
  Gamma::Algebra::SigmaZT
});

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  
  const int Nmom=7;
  int nt = latt_size[Tp];
  uint64_t vol = 1;
  for(int d=0;d<Nd;d++){
    vol = vol*latt_size[d];
  }
  
  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG(&Grid);
  pRNG.SeedFixedIntegers(seeds);


  int Nm = atoi(argv[1]); // number of all modes (high + low)

  std::vector<LatticeFermion> v(Nm,&Grid);
  std::vector<LatticeFermion> w(Nm,&Grid);
  std::vector<LatticeFermion> gammaV(Nm,&Grid);
  std::vector<LatticeComplex> phases(Nmom,&Grid);

  for(int i=0;i<Nm;i++) { 
    random(pRNG,v[i]);
    random(pRNG,w[i]);
  }

  for(int i=0;i<Nmom;i++) { 
    phases[i] = Complex(1.0);
  }

  double flops = vol * (11.0 * 8.0 + 6.0) * Nm*Nm;
  double byte  = vol * (12.0 * sizeof(Complex) ) * Nm*Nm;

  std::vector<ComplexD> ip(nt);
  std::vector<std::vector<ComplexD> > MesonFields   (Nm*Nm);
  std::vector<std::vector<ComplexD> > MesonFields4  (Nm*Nm*4);
  std::vector<std::vector<ComplexD> > MesonFields16 (Nm*Nm*16);
  std::vector<std::vector<ComplexD> > MesonFields161(Nm*Nm*16);
  std::vector<std::vector<ComplexD> > MesonFields16mom (Nm*Nm*16*Nmom);
  std::vector<std::vector<ComplexD> > MesonFieldsRef(Nm*Nm);

  for(int i=0;i<MesonFields.size();i++   )  MesonFields   [i].resize(nt);
  for(int i=0;i<MesonFieldsRef.size();i++)  MesonFieldsRef[i].resize(nt);
  for(int i=0;i<MesonFields4.size();i++  )  MesonFields4  [i].resize(nt);
  for(int i=0;i<MesonFields16.size();i++ )  MesonFields16 [i].resize(nt);
  for(int i=0;i<MesonFields161.size();i++ ) MesonFields161[i].resize(nt);

  for(int i=0;i<MesonFields16mom.size();i++ ) MesonFields16mom [i].resize(nt);

  GridLogMessage.TimingMode(1);

  std::cout<<GridLogMessage << "Running loop with sliceInnerProductVector"<<std::endl;
  double t0 = usecond();
  for(int i=0;i<Nm;i++) { 
  for(int j=0;j<Nm;j++) { 
    sliceInnerProductVector(ip, w[i],v[j],Tp);
    for(int t=0;t<nt;t++){
      MesonFieldsRef[i+j*Nm][t] = ip[t];
    }
  }}
  double t1 = usecond();
  std::cout<<GridLogMessage << "Done "<< (t1-t0) <<" usecond " <<std::endl;
  std::cout<<GridLogMessage << "Done "<< flops/(t1-t0) <<" mflops " <<std::endl;
  std::cout<<GridLogMessage << "Done "<< byte /(t1-t0) <<" MB/s " <<std::endl;

  std::cout<<GridLogMessage << "Running loop with new code for Nt="<<nt<<std::endl;
  t0 = usecond();
  sliceInnerProductMesonField(MesonFields,w,v,Tp);
  t1 = usecond();
  std::cout<<GridLogMessage << "Done "<< (t1-t0) <<" usecond " <<std::endl;
  std::cout<<GridLogMessage << "Done "<< flops/(t1-t0) <<" mflops " <<std::endl;
  std::cout<<GridLogMessage << "Done "<< byte /(t1-t0) <<" MB/s " <<std::endl;


  std::cout<<GridLogMessage << "Running loop with Four gammas code for Nt="<<nt<<std::endl;
  flops = vol * (11.0 * 8.0 + 6.0) * Nm*Nm*4;
  byte  = vol * (12.0 * sizeof(Complex) ) * Nm*Nm
        + vol * ( 2.0 * sizeof(Complex) ) * Nm*Nm* 4;
  t0 = usecond();
  sliceInnerProductMesonFieldGamma(MesonFields4,w,v,Tp,Gmu4);
  t1 = usecond();
  std::cout<<GridLogMessage << "Done "<< (t1-t0) <<" usecond " <<std::endl;
  std::cout<<GridLogMessage << "Done "<< flops/(t1-t0) <<" mflops " <<std::endl;
  std::cout<<GridLogMessage << "Done "<< byte /(t1-t0) <<" MB/s " <<std::endl;

  std::cout<<GridLogMessage << "Running loop with Sixteen gammas code for Nt="<<nt<<std::endl;
  flops = vol * (11.0 * 8.0 + 6.0) * Nm*Nm*16;
  byte  = vol * (12.0 * sizeof(Complex) ) * Nm*Nm
        + vol * ( 2.0 * sizeof(Complex) ) * Nm*Nm* 16;
  t0 = usecond();
  sliceInnerProductMesonFieldGamma(MesonFields16,w,v,Tp,Gmu16);
  t1 = usecond();
  std::cout<<GridLogMessage << "Done "<< (t1-t0) <<" usecond " <<std::endl;
  std::cout<<GridLogMessage << "Done "<< flops/(t1-t0) <<" mflops " <<std::endl;
  std::cout<<GridLogMessage << "Done "<< byte /(t1-t0) <<" MB/s " <<std::endl;


  std::cout<<GridLogMessage << "Running loop with Sixteen gammas code1 for Nt="<<nt<<std::endl;
  flops = vol * ( 2 * 8.0 + 6.0) * Nm*Nm*16;
  byte  = vol * (12.0 * sizeof(Complex) ) * Nm*Nm
        + vol * ( 2.0 * sizeof(Complex) ) * Nm*Nm* 16;
  t0 = usecond();
  sliceInnerProductMesonFieldGamma1(MesonFields161, w, v, Tp, Gmu16);
  t1 = usecond();
  std::cout<<GridLogMessage << "Done "<< (t1-t0) <<" usecond " <<std::endl;
  std::cout<<GridLogMessage << "Done "<< flops/(t1-t0) <<" mflops " <<std::endl;
  std::cout<<GridLogMessage << "Done "<< byte /(t1-t0) <<" MB/s " <<std::endl;

  std::cout<<GridLogMessage << "Running loop with Sixteen gammas "<<Nmom<<" momenta "<<std::endl;
  flops = vol * ( 2 * 8.0 + 6.0 + 8.0*Nmom) * Nm*Nm*16;
  byte  = vol * (12.0 * sizeof(Complex) ) * Nm*Nm
        + vol * ( 2.0 * sizeof(Complex) *Nmom ) * Nm*Nm* 16;
  t0 = usecond();
  sliceInnerProductMesonFieldGammaMom(MesonFields16mom,w,v,Tp,Gmu16,phases);
  t1 = usecond();
  std::cout<<GridLogMessage << "Done "<< (t1-t0) <<" usecond " <<std::endl;
  std::cout<<GridLogMessage << "Done "<< flops/(t1-t0) <<" mflops " <<std::endl;
  std::cout<<GridLogMessage << "Done "<< byte /(t1-t0) <<" MB/s " <<std::endl;



  RealD err = 0;
  RealD err2 = 0;
  ComplexD diff;
  ComplexD diff2;

  for(int i=0;i<Nm;i++) { 
  for(int j=0;j<Nm;j++) { 
    for(int t=0;t<nt;t++){
      diff = MesonFields[i+Nm*j][t] - MesonFieldsRef[i+Nm*j][t];
      err += real(diff*conj(diff));
    }
  }}
  std::cout<<GridLogMessage << "Norm error "<< err <<std::endl;
  
  err = err*0.;
  diff = diff*0.;

  for (int mu = 0; mu < 16; mu++){
    for (int k = 0; k < gammaV.size(); k++){
      gammaV[k] = Gamma(Gmu16[mu]) * v[k];
    }
    for (int i = 0; i < Nm; i++){
      for (int j = 0; j < Nm; j++){
        sliceInnerProductVector(ip, w[i], gammaV[j], Tp);
        for (int t = 0; t < nt; t++){
          MesonFields[i + j * Nm][t] = ip[t];
          diff = MesonFields16[mu+i*16+Nm*16*j][t] - MesonFields161[mu+i*16+Nm*16*j][t];
          diff2 = MesonFields[i+j*Nm][t] - MesonFields161[mu+i*16+Nm*16*j][t];
          err += real(diff*conj(diff));
          err2 += real(diff2*conj(diff2));
        }
      }
    }
  }
  std::cout << GridLogMessage << "Norm error 16 gamma1/16 gamma naive    " << err << std::endl;
  std::cout << GridLogMessage << "Norm error 16 gamma1/sliceInnerProduct " << err2 << std::endl;

  Grid_finalize();
}

