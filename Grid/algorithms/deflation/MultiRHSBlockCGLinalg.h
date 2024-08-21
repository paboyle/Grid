/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: MultiRHSBlockCGLinalg.h

    Copyright (C) 2024

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
#pragma once

NAMESPACE_BEGIN(Grid);


/* Need helper object for BLAS accelerated mrhs blockCG */
template<class Field>
class MultiRHSBlockCGLinalg
{
public:

  typedef typename Field::scalar_type   scalar;
  typedef typename Field::scalar_object scalar_object;

  deviceVector<scalar> BLAS_X;      // nrhs x vol -- the sources
  deviceVector<scalar> BLAS_Y;      // nrhs x vol -- the result
  deviceVector<scalar> BLAS_C;      // nrhs x nrhs -- the coefficients 
  
  MultiRHSBlockCGLinalg() {};
  ~MultiRHSBlockCGLinalg(){ Deallocate(); };
  
  void Deallocate(void)
  {
    BLAS_C.resize(0);
    BLAS_X.resize(0);
    BLAS_Y.resize(0);
  }
  void MaddMatrix(std::vector<Field> &AP, Eigen::MatrixXcd &m , const std::vector<Field> &X,const std::vector<Field> &Y,RealD scale=1.0)
  {
    std::vector<Field> Y_copy(AP.size(),AP[0].Grid());
    for(int r=0;r<AP.size();r++){
      Y_copy[r] = Y[r];
    }
    MulMatrix(AP,m,X);
    for(int r=0;r<AP.size();r++){
      AP[r] = scale*AP[r]+Y_copy[r];
    }
  }
  void MulMatrix(std::vector<Field> &Y, Eigen::MatrixXcd &m , const std::vector<Field> &X)
  {
    typedef typename Field::scalar_type scomplex;
    GridBase *grid;
    uint64_t vol;
    uint64_t words;

    int nrhs = Y.size();
    grid  = X[0].Grid();
    vol   = grid->lSites();
    words = sizeof(scalar_object)/sizeof(scalar);
    int64_t vw = vol * words;

    RealD t0 = usecond();
    BLAS_X.resize(nrhs * vw); // cost free if size doesn't change
    BLAS_Y.resize(nrhs * vw); // cost free if size doesn't change
    BLAS_C.resize(nrhs * nrhs);// cost free if size doesn't change
    RealD t1 = usecond();

    /////////////////////////////////////////////
    // Copy in the multi-rhs sources
    /////////////////////////////////////////////
    for(int r=0;r<nrhs;r++){
      int64_t offset = r*vw;
      autoView(x_v,X[r],AcceleratorRead);
      acceleratorCopyDeviceToDevice(&x_v[0],&BLAS_X[offset],sizeof(scalar_object)*vol);
    }

    // Assumes Eigen storage contiguous
    acceleratorCopyToDevice(&m(0,0),&BLAS_C[0],BLAS_C.size()*sizeof(scalar));
    
  /*
   * in Fortran column major notation (cuBlas order)
   *
   * Xxr = [X1(x)][..][Xn(x)]
   * Yxr = [Y1(x)][..][Ym(x)]
   * Y = X . C
   */
    deviceVector<scalar *> Xd(1);
    deviceVector<scalar *> Yd(1);
    deviceVector<scalar *> Cd(1);

    scalar * Xh = & BLAS_X[0];
    scalar * Yh = & BLAS_Y[0];
    scalar * Ch = & BLAS_C[0];

    acceleratorPut(Xd[0],Xh);
    acceleratorPut(Yd[0],Yh);
    acceleratorPut(Cd[0],Ch);

    RealD t2 = usecond();
    GridBLAS BLAS;
    /////////////////////////////////////////
    // Y = X*C (transpose?)
    /////////////////////////////////////////
    BLAS.gemmBatched(GridBLAS_OP_N,GridBLAS_OP_N, 
    		     vw,nrhs,nrhs,
		     ComplexD(1.0),
		     Xd,
		     Cd,
		     ComplexD(0.0),  // wipe out Y
		     Yd);
    BLAS.synchronise();
    RealD t3 = usecond();

    // Copy back Y = m X 
    for(int r=0;r<nrhs;r++){
      int64_t offset = r*vw;
      autoView(y_v,Y[r],AcceleratorWrite);
      acceleratorCopyDeviceToDevice(&BLAS_Y[offset],&y_v[0],sizeof(scalar_object)*vol);
    }    
    RealD t4 = usecond();
    std::cout << "MulMatrix alloc    took "<< t1-t0<<" us"<<std::endl;
    std::cout << "MulMatrix preamble took "<< t2-t1<<" us"<<std::endl;
    std::cout << "MulMatrix blas     took "<< t3-t2<<" us"<<std::endl;
    std::cout << "MulMatrix copy     took "<< t4-t3<<" us"<<std::endl;
    std::cout << "MulMatrix total "<< t4-t0<<" us"<<std::endl;
  }
  
  void InnerProductMatrix(Eigen::MatrixXcd &m , const std::vector<Field> &X, const std::vector<Field> &Y)
  {
    int nrhs;
    GridBase *grid;
    uint64_t vol;
    uint64_t words;

    nrhs = X.size();
    assert(X.size()==Y.size());
    conformable(X[0],Y[0]);

    grid  = X[0].Grid();
    vol   = grid->lSites();
    words = sizeof(scalar_object)/sizeof(scalar);
    int64_t vw = vol * words;

    RealD t0 = usecond();
    BLAS_X.resize(nrhs * vw); // cost free if size doesn't change
    BLAS_Y.resize(nrhs * vw); // cost free if size doesn't change
    BLAS_C.resize(nrhs * nrhs);// cost free if size doesn't change
    RealD t1 = usecond();

    /////////////////////////////////////////////
    // Copy in the multi-rhs sources
    /////////////////////////////////////////////
    for(int r=0;r<nrhs;r++){
      int64_t offset = r*vw;
      autoView(x_v,X[r],AcceleratorRead);
      acceleratorCopyDeviceToDevice(&x_v[0],&BLAS_X[offset],sizeof(scalar_object)*vol);
      autoView(y_v,Y[r],AcceleratorRead);
      acceleratorCopyDeviceToDevice(&y_v[0],&BLAS_Y[offset],sizeof(scalar_object)*vol);
    }
    RealD t2 = usecond();

  /*
   * in Fortran column major notation (cuBlas order)
   *
   * Xxr = [X1(x)][..][Xn(x)]
   *
   * Yxr = [Y1(x)][..][Ym(x)]
   *
   * C_rs = X^dag Y
   */
    deviceVector<scalar *> Xd(1);
    deviceVector<scalar *> Yd(1);
    deviceVector<scalar *> Cd(1);

    scalar * Xh = & BLAS_X[0];
    scalar * Yh = & BLAS_Y[0];
    scalar * Ch = & BLAS_C[0];

    acceleratorPut(Xd[0],Xh);
    acceleratorPut(Yd[0],Yh);
    acceleratorPut(Cd[0],Ch);

    GridBLAS BLAS;

    RealD t3 = usecond();
    /////////////////////////////////////////
    // C_rs = X^dag Y
    /////////////////////////////////////////
    BLAS.gemmBatched(GridBLAS_OP_C,GridBLAS_OP_N, 
    		     nrhs,nrhs,vw,
		     ComplexD(1.0),
		     Xd,
		     Yd,
		     ComplexD(0.0),  // wipe out C
		     Cd);
    BLAS.synchronise();
    RealD t4 = usecond();

    std::vector<scalar> HOST_C(BLAS_C.size());      // nrhs . nrhs -- the coefficients 
    acceleratorCopyFromDevice(&BLAS_C[0],&HOST_C[0],BLAS_C.size()*sizeof(scalar));
    grid->GlobalSumVector(&HOST_C[0],nrhs*nrhs);

    RealD t5 = usecond();
    for(int rr=0;rr<nrhs;rr++){
      for(int r=0;r<nrhs;r++){
	int off = r+nrhs*rr;
	m(r,rr)=HOST_C[off];
      }
    }
    RealD t6 = usecond();
    uint64_t M=nrhs;
    uint64_t N=nrhs;
    uint64_t K=vw;
    RealD bytes = 1.0*sizeof(ComplexD)*(M*N*2+N*K+M*K);
    RealD flops = 8.0*M*N*K;
    flops = flops/(t4-t3)/1.e3;
    bytes = bytes/(t4-t3)/1.e3;
    std::cout << "InnerProductMatrix m,n,k "<< M<<","<<N<<","<<K<<std::endl;
    std::cout << "InnerProductMatrix alloc t1 "<< t1-t0<<" us"<<std::endl;
    std::cout << "InnerProductMatrix cp    t2 "<< t2-t1<<" us"<<std::endl;
    std::cout << "InnerProductMatrix setup t3 "<< t3-t2<<" us"<<std::endl;
    std::cout << "InnerProductMatrix blas t4 "<< t4-t3<<" us"<<std::endl;
    std::cout << "InnerProductMatrix blas    "<< flops<<" GF/s"<<std::endl;
    std::cout << "InnerProductMatrix blas    "<< bytes<<" GB/s"<<std::endl;
    std::cout << "InnerProductMatrix gsum t5 "<< t5-t4<<" us"<<std::endl;
    std::cout << "InnerProductMatrix cp   t6 "<< t6-t5<<" us"<<std::endl;
    std::cout << "InnerProductMatrix took "<< t6-t0<<" us"<<std::endl;

  }
};

NAMESPACE_END(Grid);
