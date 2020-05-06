/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/lattice/Lattice_basis.h

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Christoph Lehner <christoph@lhnr.de>

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

#pragma once

NAMESPACE_BEGIN(Grid);

template<class Field>
void basisOrthogonalize(std::vector<Field> &basis,Field &w,int k) 
{
  // If assume basis[j] are already orthonormal,
  // can take all inner products in parallel saving 2x bandwidth
  // Save 3x bandwidth on the second line of loop.
  // perhaps 2.5x speed up.
  // 2x overall in Multigrid Lanczos  
  for(int j=0; j<k; ++j){
    auto ip = innerProduct(basis[j],w);
    w = w - ip*basis[j];
  }
}

template<class VField, class Matrix>
void basisRotate(VField &basis,Matrix& Qt,int j0, int j1, int k0,int k1,int Nm) 
{
  typedef decltype(basis[0]) Field;
  typedef decltype(basis[0].View()) View;
  auto tmp_v = basis[0].AcceleratorView(ViewReadWrite);
  Vector<View> basis_v(basis.size(),tmp_v);
  typedef typename std::remove_reference<decltype(tmp_v[0])>::type vobj;
  GridBase* grid = basis[0].Grid();
      
  for(int k=0;k<basis.size();k++){
    basis_v[k] = basis[k].AcceleratorView(ViewReadWrite);
  }

#ifndef GRID_NVCC
  thread_region
  {
    std::vector < vobj > B(Nm); // Thread private
    thread_for_in_region(ss, grid->oSites(),{
	for(int j=j0; j<j1; ++j) B[j]=0.;
      
	for(int j=j0; j<j1; ++j){
	  for(int k=k0; k<k1; ++k){
	    B[j] +=Qt(j,k) * basis_v[k][ss];
	  }
	}
	for(int j=j0; j<j1; ++j){
	  basis_v[j][ss] = B[j];
	}
      });
  }
#else
  int nrot = j1-j0;
  if (!nrot) // edge case not handled gracefully by Cuda
    return;

  uint64_t oSites   =grid->oSites();
  uint64_t siteBlock=(grid->oSites()+nrot-1)/nrot; // Maximum 1 additional vector overhead

  Vector <vobj> Bt(siteBlock * nrot); 
  auto Bp=&Bt[0];

  // GPU readable copy of matrix
  Vector<double> Qt_jv(Nm*Nm);
  double *Qt_p = & Qt_jv[0];
  thread_for(i,Nm*Nm,{
      int j = i/Nm;
      int k = i%Nm;
      Qt_p[i]=Qt(j,k);
    });

  // Block the loop to keep storage footprint down
  for(uint64_t s=0;s<oSites;s+=siteBlock){

    // remaining work in this block
    int ssites=MIN(siteBlock,oSites-s);

    // zero out the accumulators
    accelerator_for(ss,siteBlock*nrot,vobj::Nsimd(),{
	decltype(coalescedRead(Bp[ss])) z;
	z=Zero();
	coalescedWrite(Bp[ss],z);
      });

    accelerator_for(sj,ssites*nrot,vobj::Nsimd(),{
	
	int j =sj%nrot;
	int jj  =j0+j;
	int ss =sj/nrot;
	int sss=ss+s;

	for(int k=k0; k<k1; ++k){
	  auto tmp = coalescedRead(Bp[ss*nrot+j]);
	  coalescedWrite(Bp[ss*nrot+j],tmp+ Qt_p[jj*Nm+k] * coalescedRead(basis_v[k][sss]));
	}
      });

    accelerator_for(sj,ssites*nrot,vobj::Nsimd(),{
	int j =sj%nrot;
	int jj  =j0+j;
	int ss =sj/nrot;
	int sss=ss+s;
	coalescedWrite(basis_v[jj][sss],coalescedRead(Bp[ss*nrot+j]));
      });
  }
#endif
}

// Extract a single rotated vector
template<class Field>
void basisRotateJ(Field &result,std::vector<Field> &basis,Eigen::MatrixXd& Qt,int j, int k0,int k1,int Nm) 
{
  typedef decltype(basis[0].AcceleratorView()) View;
  typedef typename Field::vector_object vobj;
  GridBase* grid = basis[0].Grid();

  result.Checkerboard() = basis[0].Checkerboard();
  auto result_v=result.AcceleratorView(ViewWrite);
  Vector<View> basis_v(basis.size(),result_v);
  for(int k=0;k<basis.size();k++){
    basis_v[k] = basis[k].AcceleratorView(ViewRead);
  }
  vobj zz=Zero();
  Vector<double> Qt_jv(Nm);
  double * Qt_j = & Qt_jv[0];
  for(int k=0;k<Nm;++k) Qt_j[k]=Qt(j,k);
  accelerator_for(ss, grid->oSites(),vobj::Nsimd(),{
    auto B=coalescedRead(zz);
    for(int k=k0; k<k1; ++k){
      B +=Qt_j[k] * coalescedRead(basis_v[k][ss]);
    }
    coalescedWrite(result_v[ss], B);
  });
}

template<class Field>
void basisReorderInPlace(std::vector<Field> &_v,std::vector<RealD>& sort_vals, std::vector<int>& idx) 
{
  int vlen = idx.size();

  assert(vlen>=1);
  assert(vlen<=sort_vals.size());
  assert(vlen<=_v.size());

  for (size_t i=0;i<vlen;i++) {

    if (idx[i] != i) {

      //////////////////////////////////////
      // idx[i] is a table of desired sources giving a permutation.
      // Swap v[i] with v[idx[i]].
      // Find  j>i for which _vnew[j] = _vold[i],
      // track the move idx[j] => idx[i]
      // track the move idx[i] => i
      //////////////////////////////////////
      size_t j;
      for (j=i;j<idx.size();j++)
	if (idx[j]==i)
	  break;

      assert(idx[i] > i);     assert(j!=idx.size());      assert(idx[j]==i);

      swap(_v[i],_v[idx[i]]); // should use vector move constructor, no data copy
      std::swap(sort_vals[i],sort_vals[idx[i]]);

      idx[j] = idx[i];
      idx[i] = i;
    }
  }
}

inline std::vector<int> basisSortGetIndex(std::vector<RealD>& sort_vals) 
{
  std::vector<int> idx(sort_vals.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(), [&sort_vals](int i1, int i2) {
    return ::fabs(sort_vals[i1]) < ::fabs(sort_vals[i2]);
  });
  return idx;
}

template<class Field>
void basisSortInPlace(std::vector<Field> & _v,std::vector<RealD>& sort_vals, bool reverse) 
{
  std::vector<int> idx = basisSortGetIndex(sort_vals);
  if (reverse)
    std::reverse(idx.begin(), idx.end());
  
  basisReorderInPlace(_v,sort_vals,idx);
}

// PAB: faster to compute the inner products first then fuse loops.
// If performance critical can improve.
template<class Field>
void basisDeflate(const std::vector<Field> &_v,const std::vector<RealD>& eval,const Field& src_orig,Field& result) {
  result = Zero();
  assert(_v.size()==eval.size());
  int N = (int)_v.size();
  for (int i=0;i<N;i++) {
    Field& tmp = _v[i];
    axpy(result,TensorRemove(innerProduct(tmp,src_orig)) / eval[i],tmp,result);
  }
}

NAMESPACE_END(Grid);
