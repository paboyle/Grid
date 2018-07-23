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
  for (int r = 0; r < rd * Lblock * Rblock; r++){
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
#if 1	    
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
#else
	    vector_type vv = TensorRemove(innerProduct(left,right));
#endif
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

/*
template void sliceInnerProductMesonField<SpinColourVector>(std::vector< std::vector<ComplexD> > &mat, 
						   const std::vector<Lattice<SpinColourVector> > &lhs,
						   const std::vector<Lattice<SpinColourVector> > &rhs,
						   int orthogdim) ;
*/

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);

  int nt = latt_size[Tp];
  uint64_t vol = 1;
  for(int d=0;d<Nd;d++){
    vol = vol*latt_size[d];
  }
  
  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG(&Grid);
  pRNG.SeedFixedIntegers(seeds);


  const int Nm = 32; // number of all modes (high + low)

  std::vector<LatticeFermion> v(Nm,&Grid);
  std::vector<LatticeFermion> w(Nm,&Grid);

  for(int i=0;i<Nm;i++) { 
    random(pRNG,v[i]);
    random(pRNG,w[i]);
  }

  double flops = vol * (11.0 * 8.0 + 6.0) * Nm*Nm;
  double byte  = vol * (12.0 * sizeof(Complex) ) * Nm*Nm;

  std::vector<ComplexD> ip(nt);
  std::vector<std::vector<ComplexD> > MesonFields   (Nm*Nm);
  std::vector<std::vector<ComplexD> > MesonFieldsRef(Nm*Nm);

  for(int i=0;i<Nm;i++) { 
  for(int j=0;j<Nm;j++) { 
    MesonFields   [i+j*Nm].resize(nt);
    MesonFieldsRef[i+j*Nm].resize(nt);
  }}

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
  std::cout<<GridLogMessage << "Done "<< flops/(t1-t0) <<" mflops " <<std::endl;
  std::cout<<GridLogMessage << "Done "<< byte /(t1-t0) <<" MB/s " <<std::endl;

  std::cout<<GridLogMessage << "Running loop with new code for Nt="<<nt<<std::endl;
  double t2 = usecond();
  sliceInnerProductMesonField(MesonFields,w,v,Tp);
  double t3 = usecond();
  std::cout<<GridLogMessage << "Done "<< flops/(t3-t2) <<" mflops " <<std::endl;
  std::cout<<GridLogMessage << "Done "<< byte /(t3-t2) <<" MB/s " <<std::endl;

  RealD err = 0;
  ComplexD diff;

  for(int i=0;i<Nm;i++) { 
  for(int j=0;j<Nm;j++) { 
    for(int t=0;t<nt;t++){
      diff = MesonFields[i+Nm*j][t] - MesonFieldsRef[i+Nm*j][t];
      err += real(diff*conj(diff));
    }
  }}
  std::cout<<GridLogMessage << "Norm error "<< err <<std::endl;
  
  Grid_finalize();
}

