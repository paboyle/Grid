/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/debug/Test_batched_blas.cc

    Copyright (C) 2026

Author: Peter Boyle <pboyle@bnl.gov>

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

//
// Unit tests for the blas-layer dense primitives:
//
//  T1 : GridBLASInverse::inverseBatched (ComplexF and ComplexD) --
//       || A A^{-1} - I ||_max over a random well-conditioned batch.
//       On a CPU build this exercises the Eigen reference (the semantic
//       oracle); the SAME binary exercises rocSOLVER/cuBLAS/oneMKL on a
//       device build.
//  T2 : explicit-leading-dimension gemmBatched == SOFTWARE SPLIT-K in
//       miniature.  Y = op(A,T).X computed (a) in one compact batch=1 call
//       and (b) as NK K-chunks by POINTER OFFSET with lda = full K, partials
//       reduced on the host.  (a)==(b) validates the strided overload AND
//       the production dense-slab apply pattern (arXiv:2409.03904 fig 11).
//
// Hard asserts throughout: regression gate for Grid/algorithms/blas.
//
#include <Grid/Grid.h>
#include <Grid/algorithms/blas/BatchedInverse.h>

#include <random>

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  GridBLAS        blas;
  GridBLASInverse inverse;

  std::mt19937 rng(12345);
  std::uniform_real_distribution<double> dist(-1.0,1.0);

  ////////////////////////////////////////////////////////////////
  // T1a : batched inversion, ComplexF
  ////////////////////////////////////////////////////////////////
  {
    const int64_t N = 64;
    const int  batch = 4;
    const uint64_t elems = (uint64_t)batch*N*N;

    // Random diagonally-dominant batch: A = N*I + R, |R_ij| <= 1
    std::vector<ComplexF> Ahost(elems);
    for(uint64_t i=0;i<elems;i++) Ahost[i] = ComplexF(dist(rng),dist(rng));
    for(int b=0;b<batch;b++)
      for(int64_t d=0;d<N;d++)
	Ahost[(uint64_t)b*N*N + d*N + d] += ComplexF((RealF)N,0.0);

    deviceVector<ComplexF> Adev(elems);   // gets inverted in place
    deviceVector<ComplexF> Aorig(elems);  // untouched copy for the residual
    deviceVector<ComplexF> Cdev(elems);
    acceleratorCopyToDevice(&Ahost[0],&Adev[0], elems*sizeof(ComplexF));
    acceleratorCopyToDevice(&Ahost[0],&Aorig[0],elems*sizeof(ComplexF));

    deviceVector<ComplexF*> Ap(batch);
    deviceVector<ComplexF*> Op(batch);
    deviceVector<ComplexF*> Cp(batch);
    std::vector<ComplexF*> ptr_h(batch);
    for(int b=0;b<batch;b++) ptr_h[b] = &Adev[(uint64_t)b*N*N];
    acceleratorCopyToDevice(&ptr_h[0],&Ap[0],batch*sizeof(ComplexF*));
    for(int b=0;b<batch;b++) ptr_h[b] = &Aorig[(uint64_t)b*N*N];
    acceleratorCopyToDevice(&ptr_h[0],&Op[0],batch*sizeof(ComplexF*));
    for(int b=0;b<batch;b++) ptr_h[b] = &Cdev[(uint64_t)b*N*N];
    acceleratorCopyToDevice(&ptr_h[0],&Cp[0],batch*sizeof(ComplexF*));

    inverse.inverseBatched(N,Ap);                      // A <- A^{-1}

    ComplexF one (1.0,0.0);
    ComplexF zero(0.0,0.0);
    blas.gemmBatched(GridBLAS_OP_N,GridBLAS_OP_N,
		     (int)N,(int)N,(int)N,
		     one, Op, Ap, zero, Cp);           // C = A A^{-1}
    blas.synchronise();

    std::vector<ComplexF> Chost(elems);
    acceleratorCopyFromDevice(&Cdev[0],&Chost[0],elems*sizeof(ComplexF));

    double maxdev = 0.0;
    for(int b=0;b<batch;b++)
      for(int64_t j=0;j<N;j++)
	for(int64_t i=0;i<N;i++){
	  ComplexF expect = (i==j) ? ComplexF(1.0,0.0) : ComplexF(0.0,0.0);
	  ComplexF got = Chost[(uint64_t)b*N*N + j*N + i];
	  maxdev = std::max(maxdev,(double)abs(got-expect));
	}
    std::cout << GridLogMessage << "T1a inverseBatched ComplexF  ||A Ainv - I||_max = "
	      << maxdev << ( maxdev < 1.0e-4 ? "   PASS" : "   FAIL" ) << std::endl;
    GRID_ASSERT(maxdev < 1.0e-4);
  }

  ////////////////////////////////////////////////////////////////
  // T1b : batched inversion, ComplexD
  ////////////////////////////////////////////////////////////////
  {
    const int64_t N = 48;
    const int  batch = 3;
    const uint64_t elems = (uint64_t)batch*N*N;

    std::vector<ComplexD> Ahost(elems);
    for(uint64_t i=0;i<elems;i++) Ahost[i] = ComplexD(dist(rng),dist(rng));
    for(int b=0;b<batch;b++)
      for(int64_t d=0;d<N;d++)
	Ahost[(uint64_t)b*N*N + d*N + d] += ComplexD((RealD)N,0.0);

    deviceVector<ComplexD> Adev(elems);
    deviceVector<ComplexD> Aorig(elems);
    deviceVector<ComplexD> Cdev(elems);
    acceleratorCopyToDevice(&Ahost[0],&Adev[0], elems*sizeof(ComplexD));
    acceleratorCopyToDevice(&Ahost[0],&Aorig[0],elems*sizeof(ComplexD));

    deviceVector<ComplexD*> Ap(batch);
    deviceVector<ComplexD*> Op(batch);
    deviceVector<ComplexD*> Cp(batch);
    std::vector<ComplexD*> ptr_h(batch);
    for(int b=0;b<batch;b++) ptr_h[b] = &Adev[(uint64_t)b*N*N];
    acceleratorCopyToDevice(&ptr_h[0],&Ap[0],batch*sizeof(ComplexD*));
    for(int b=0;b<batch;b++) ptr_h[b] = &Aorig[(uint64_t)b*N*N];
    acceleratorCopyToDevice(&ptr_h[0],&Op[0],batch*sizeof(ComplexD*));
    for(int b=0;b<batch;b++) ptr_h[b] = &Cdev[(uint64_t)b*N*N];
    acceleratorCopyToDevice(&ptr_h[0],&Cp[0],batch*sizeof(ComplexD*));

    inverse.inverseBatched(N,Ap);

    ComplexD one (1.0,0.0);
    ComplexD zero(0.0,0.0);
    blas.gemmBatched(GridBLAS_OP_N,GridBLAS_OP_N,
		     (int)N,(int)N,(int)N,
		     one, Op, Ap, zero, Cp);
    blas.synchronise();

    std::vector<ComplexD> Chost(elems);
    acceleratorCopyFromDevice(&Cdev[0],&Chost[0],elems*sizeof(ComplexD));

    double maxdev = 0.0;
    for(int b=0;b<batch;b++)
      for(int64_t j=0;j<N;j++)
	for(int64_t i=0;i<N;i++){
	  ComplexD expect = (i==j) ? ComplexD(1.0,0.0) : ComplexD(0.0,0.0);
	  ComplexD got = Chost[(uint64_t)b*N*N + j*N + i];
	  maxdev = std::max(maxdev,(double)abs(got-expect));
	}
    std::cout << GridLogMessage << "T1b inverseBatched ComplexD  ||A Ainv - I||_max = "
	      << maxdev << ( maxdev < 1.0e-10 ? "   PASS" : "   FAIL" ) << std::endl;
    GRID_ASSERT(maxdev < 1.0e-10);
  }

  ////////////////////////////////////////////////////////////////
  // T2 : explicit-ld gemmBatched == software split-K miniature.
  // Slab layout as in the dense coarse-coarse apply: A is K x nrows
  // column major (lda=K); Y = op(A,T).X with X K x nrhs (ldb=K).
  ////////////////////////////////////////////////////////////////
  {
    const int nrows = 8;
    const int nrhs  = 4;
    const int K     = 256;
    const int NK    = 8;        // split-K chunks
    const int Kc    = K/NK;

    std::vector<ComplexF> Ahost((uint64_t)K*nrows);
    std::vector<ComplexF> Xhost((uint64_t)K*nrhs);
    for(auto &z : Ahost) z = ComplexF(dist(rng),dist(rng));
    for(auto &z : Xhost) z = ComplexF(dist(rng),dist(rng));

    deviceVector<ComplexF> Adev(Ahost.size());
    deviceVector<ComplexF> Xdev(Xhost.size());
    deviceVector<ComplexF> Yref((uint64_t)nrows*nrhs);
    deviceVector<ComplexF> Ypart((uint64_t)NK*nrows*nrhs);
    acceleratorCopyToDevice(&Ahost[0],&Adev[0],Ahost.size()*sizeof(ComplexF));
    acceleratorCopyToDevice(&Xhost[0],&Xdev[0],Xhost.size()*sizeof(ComplexF));

    ComplexF one (1.0,0.0);
    ComplexF zero(0.0,0.0);

    // (a) reference: one compact batch=1 call (compact lda == K for OP_T)
    {
      deviceVector<ComplexF*> Ap(1), Xp(1), Yp(1);
      std::vector<ComplexF*> h(1);
      h[0]=&Adev[0]; acceleratorCopyToDevice(&h[0],&Ap[0],sizeof(ComplexF*));
      h[0]=&Xdev[0]; acceleratorCopyToDevice(&h[0],&Xp[0],sizeof(ComplexF*));
      h[0]=&Yref[0]; acceleratorCopyToDevice(&h[0],&Yp[0],sizeof(ComplexF*));
      blas.gemmBatched(GridBLAS_OP_T,GridBLAS_OP_N,
		       nrows,nrhs,K,
		       one, Ap, Xp, zero, Yp);
      blas.synchronise();
    }

    // (b) split-K: NK chunk-pointers into the SAME allocations, lda/ldb = K
    {
      deviceVector<ComplexF*> Ap(NK), Xp(NK), Yp(NK);
      std::vector<ComplexF*> h(NK);
      for(int j=0;j<NK;j++) h[j] = &Adev[(uint64_t)j*Kc];        // K-offset slice
      acceleratorCopyToDevice(&h[0],&Ap[0],NK*sizeof(ComplexF*));
      for(int j=0;j<NK;j++) h[j] = &Xdev[(uint64_t)j*Kc];
      acceleratorCopyToDevice(&h[0],&Xp[0],NK*sizeof(ComplexF*));
      for(int j=0;j<NK;j++) h[j] = &Ypart[(uint64_t)j*nrows*nrhs];
      acceleratorCopyToDevice(&h[0],&Yp[0],NK*sizeof(ComplexF*));
      blas.gemmBatched(GridBLAS_OP_T,GridBLAS_OP_N,
		       nrows,nrhs,Kc,
		       one, Ap, /*lda*/ K,
		            Xp, /*ldb*/ K,
		       zero,Yp, /*ldc*/ nrows);
      blas.synchronise();
    }

    std::vector<ComplexF> Yref_h((uint64_t)nrows*nrhs);
    std::vector<ComplexF> Ypart_h((uint64_t)NK*nrows*nrhs);
    acceleratorCopyFromDevice(&Yref[0], &Yref_h[0], Yref_h.size()*sizeof(ComplexF));
    acceleratorCopyFromDevice(&Ypart[0],&Ypart_h[0],Ypart_h.size()*sizeof(ComplexF));

    double maxdev = 0.0;
    double maxval = 0.0;
    for(int i=0;i<nrows*nrhs;i++){
      ComplexF sum(0.0,0.0);
      for(int j=0;j<NK;j++) sum = sum + Ypart_h[(uint64_t)j*nrows*nrhs + i];
      maxdev = std::max(maxdev,(double)abs(sum-Yref_h[i]));
      maxval = std::max(maxval,(double)abs(Yref_h[i]));
    }
    double rel = maxdev/maxval;
    std::cout << GridLogMessage << "T2  split-K strided gemmBatched  max rel dev vs compact = "
	      << rel << ( rel < 1.0e-4 ? "   PASS" : "   FAIL" ) << std::endl;
    GRID_ASSERT(rel < 1.0e-4);
  }

  ////////////////////////////////////////////////////////////////
  // T3 : ComplexD explicit-ld gemmBatched (the RecursiveSchurInverse
  // merge primitive): same split-K miniature as T2, double precision.
  // On device builds this is the FIRST exercise of hipblasZ/cublasZ
  // gemmBatched through the strided overload.
  ////////////////////////////////////////////////////////////////
  {
    const int nrows = 8;
    const int nrhs  = 4;
    const int K     = 256;
    const int NK    = 8;        // split-K chunks
    const int Kc    = K/NK;

    std::vector<ComplexD> Ahost((uint64_t)K*nrows);
    std::vector<ComplexD> Xhost((uint64_t)K*nrhs);
    for(auto &z : Ahost) z = ComplexD(dist(rng),dist(rng));
    for(auto &z : Xhost) z = ComplexD(dist(rng),dist(rng));

    deviceVector<ComplexD> Adev(Ahost.size());
    deviceVector<ComplexD> Xdev(Xhost.size());
    deviceVector<ComplexD> Yref((uint64_t)nrows*nrhs);
    deviceVector<ComplexD> Ypart((uint64_t)NK*nrows*nrhs);
    acceleratorCopyToDevice(&Ahost[0],&Adev[0],Ahost.size()*sizeof(ComplexD));
    acceleratorCopyToDevice(&Xhost[0],&Xdev[0],Xhost.size()*sizeof(ComplexD));

    ComplexD one (1.0,0.0);
    ComplexD zero(0.0,0.0);

    // (a) reference: one compact batch=1 call (compact lda == K for OP_T)
    {
      deviceVector<ComplexD*> Ap(1), Xp(1), Yp(1);
      std::vector<ComplexD*> h(1);
      h[0]=&Adev[0]; acceleratorCopyToDevice(&h[0],&Ap[0],sizeof(ComplexD*));
      h[0]=&Xdev[0]; acceleratorCopyToDevice(&h[0],&Xp[0],sizeof(ComplexD*));
      h[0]=&Yref[0]; acceleratorCopyToDevice(&h[0],&Yp[0],sizeof(ComplexD*));
      blas.gemmBatched(GridBLAS_OP_T,GridBLAS_OP_N,
                       nrows,nrhs,K,
                       one, Ap, Xp, zero, Yp);
      blas.synchronise();
    }

    // (b) split-K: NK chunk-pointers into the SAME allocations, lda/ldb = K
    {
      deviceVector<ComplexD*> Ap(NK), Xp(NK), Yp(NK);
      std::vector<ComplexD*> h(NK);
      for(int j=0;j<NK;j++) h[j] = &Adev[(uint64_t)j*Kc];        // K-offset slice
      acceleratorCopyToDevice(&h[0],&Ap[0],NK*sizeof(ComplexD*));
      for(int j=0;j<NK;j++) h[j] = &Xdev[(uint64_t)j*Kc];
      acceleratorCopyToDevice(&h[0],&Xp[0],NK*sizeof(ComplexD*));
      for(int j=0;j<NK;j++) h[j] = &Ypart[(uint64_t)j*nrows*nrhs];
      acceleratorCopyToDevice(&h[0],&Yp[0],NK*sizeof(ComplexD*));
      blas.gemmBatched(GridBLAS_OP_T,GridBLAS_OP_N,
                       nrows,nrhs,Kc,
                       one, Ap, /*lda*/ K,
                            Xp, /*ldb*/ K,
                       zero,Yp, /*ldc*/ nrows);
      blas.synchronise();
    }

    std::vector<ComplexD> Yref_h((uint64_t)nrows*nrhs);
    std::vector<ComplexD> Ypart_h((uint64_t)NK*nrows*nrhs);
    acceleratorCopyFromDevice(&Yref[0], &Yref_h[0], Yref_h.size()*sizeof(ComplexD));
    acceleratorCopyFromDevice(&Ypart[0],&Ypart_h[0],Ypart_h.size()*sizeof(ComplexD));

    double maxdev = 0.0;
    double maxval = 0.0;
    for(int i=0;i<nrows*nrhs;i++){
      ComplexD sum(0.0,0.0);
      for(int j=0;j<NK;j++) sum = sum + Ypart_h[(uint64_t)j*nrows*nrhs + i];
      maxdev = std::max(maxdev,(double)abs(sum-Yref_h[i]));
      maxval = std::max(maxval,(double)abs(Yref_h[i]));
    }
    double rel = maxdev/maxval;
    std::cout << GridLogMessage << "T3  split-K strided gemmBatched ComplexD  max rel dev vs compact = "
              << rel << ( rel < 1.0e-13 ? "   PASS" : "   FAIL" ) << std::endl;
    GRID_ASSERT(rel < 1.0e-13);
  }

  std::cout << GridLogMessage << "All batched-blas tests PASSED" << std::endl;

  Grid_finalize();
}
