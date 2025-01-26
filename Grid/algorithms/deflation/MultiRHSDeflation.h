/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: MultiRHSDeflation.h

    Copyright (C) 2023

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


/* Need helper object for BLAS accelerated mrhs projection

   i) MultiRHS Deflation

   Import Evecs -> nev x vol x internal 
   Import vector of Lattice objects -> nrhs x vol x internal
   => Cij (nrhs x Nev) via GEMM.
   => Guess  (nrhs x vol x internal)  = C x evecs (via GEMM)
   Export

   
   ii) MultiRHS block projection

   Import basis -> nblock x nbasis x  (block x internal) 
   Import vector of fine lattice objects -> nblock x nrhs x (block x internal) 

   => coarse_(nrhs x nbasis )^block = via batched GEMM

   iii)   Alternate interface: 
   Import higher dim Lattice object-> vol x nrhs layout
   
*/
template<class Field>
class MultiRHSDeflation
{
public:

  typedef typename Field::scalar_type   scalar;
  typedef typename Field::scalar_object scalar_object;

  int nev;
  std::vector<RealD> eval;
  GridBase *grid;
  uint64_t vol;
  uint64_t words;
  
  deviceVector<scalar> BLAS_E;      //  nev x vol -- the eigenbasis   (up to a 1/sqrt(lambda))
  deviceVector<scalar> BLAS_R;      // nrhs x vol -- the sources
  deviceVector<scalar> BLAS_G;      // nrhs x vol -- the guess
  deviceVector<scalar> BLAS_C;      // nrhs x nev -- the coefficients 
  
  MultiRHSDeflation(){};
  ~MultiRHSDeflation(){ Deallocate(); };
  
  void Deallocate(void)
  {
    nev=0;
    grid=nullptr;
    vol=0;
    words=0;
    BLAS_E.resize(0);
    BLAS_R.resize(0);
    BLAS_C.resize(0);
    BLAS_G.resize(0);
  }
  void Allocate(int _nev,GridBase *_grid)
  {
    nev=_nev;
    grid=_grid;
    vol   = grid->lSites();
    words = sizeof(scalar_object)/sizeof(scalar);
    eval.resize(nev);
    BLAS_E.resize (vol * words * nev );
    std::cout << GridLogMessage << " Allocate for "<<nev<<" eigenvectors and volume "<<vol<<std::endl;
  }
  void ImportEigenVector(Field &evec,RealD &_eval, int ev)
  {
    //    std::cout << " ev " <<ev<<" eval "<<_eval<< std::endl;
    assert(ev<eval.size());
    eval[ev] = _eval;

    int64_t offset = ev*vol*words;
    autoView(v,evec,AcceleratorRead);
    acceleratorCopyDeviceToDevice(&v[0],&BLAS_E[offset],sizeof(scalar_object)*vol);

  }
  void ImportEigenBasis(std::vector<Field> &evec,std::vector<RealD> &_eval)
  {
    ImportEigenBasis(evec,_eval,0,evec.size());
  }
  // Could use to import a batch of eigenvectors
  void ImportEigenBasis(std::vector<Field> &evec,std::vector<RealD> &_eval, int _ev0, int _nev)
  {
    assert(_ev0+_nev<=evec.size());

    Allocate(_nev,evec[0].Grid());
    
    // Imports a sub-batch of eigenvectors, _ev0, ..., _ev0+_nev-1
    for(int e=0;e<nev;e++){
      std::cout << "Importing eigenvector "<<e<<" evalue "<<_eval[_ev0+e]<<std::endl;
      ImportEigenVector(evec[_ev0+e],_eval[_ev0+e],e);
    }
  }
  void DeflateSources(std::vector<Field> &source,std::vector<Field> & guess)
  {
    int nrhs = source.size();
    assert(source.size()==guess.size());
    assert(grid == guess[0].Grid());
    conformable(guess[0],source[0]);

    int64_t vw = vol * words;

    RealD t0 = usecond();
    BLAS_R.resize(nrhs * vw); // cost free if size doesn't change
    BLAS_G.resize(nrhs * vw); // cost free if size doesn't change
    BLAS_C.resize(nev * nrhs);// cost free if size doesn't change

    /////////////////////////////////////////////
    // Copy in the multi-rhs sources
    /////////////////////////////////////////////
    //    for(int r=0;r<nrhs;r++){
    //      std::cout << " source["<<r<<"] = "<<norm2(source[r])<<std::endl;
    //    }
    for(int r=0;r<nrhs;r++){
      int64_t offset = r*vw;
      autoView(v,source[r],AcceleratorRead);
      acceleratorCopyDeviceToDevice(&v[0],&BLAS_R[offset],sizeof(scalar_object)*vol);
    }

  /*
   * in Fortran column major notation (cuBlas order)
   *
   * Exe = [e1(x)][..][en(x)]
   *
   * Rxr = [r1(x)][..][rm(x)]
   *
   * C_er = E^dag R
   * C_er = C_er / lambda_e 
   * G_xr = Exe Cer
   */
    deviceVector<scalar *> Ed(1);
    deviceVector<scalar *> Rd(1);
    deviceVector<scalar *> Cd(1);
    deviceVector<scalar *> Gd(1);

    scalar * Eh = & BLAS_E[0];
    scalar * Rh = & BLAS_R[0];
    scalar * Ch = & BLAS_C[0];
    scalar * Gh = & BLAS_G[0];

    acceleratorPut(Ed[0],Eh);
    acceleratorPut(Rd[0],Rh);
    acceleratorPut(Cd[0],Ch);
    acceleratorPut(Gd[0],Gh);

    GridBLAS BLAS;

    /////////////////////////////////////////
    // C_er = E^dag R
    /////////////////////////////////////////
    BLAS.gemmBatched(GridBLAS_OP_C,GridBLAS_OP_N, 
    		     nev,nrhs,vw,
		     scalar(1.0),
		     Ed,
		     Rd,
		     scalar(0.0),  // wipe out C
		     Cd);
    BLAS.synchronise();

    assert(BLAS_C.size()==nev*nrhs);

    std::vector<scalar> HOST_C(BLAS_C.size());      // nrhs . nev -- the coefficients 
    acceleratorCopyFromDevice(&BLAS_C[0],&HOST_C[0],BLAS_C.size()*sizeof(scalar));
    grid->GlobalSumVector(&HOST_C[0],nev*nrhs);
    for(int e=0;e<nev;e++){
      RealD lam(1.0/eval[e]);
      for(int r=0;r<nrhs;r++){
	int off = e+nev*r;
	HOST_C[off]=HOST_C[off] * lam;
	//	std::cout << "C["<<e<<"]["<<r<<"] ="<<HOST_C[off]<< " eval[e] "<<eval[e] <<std::endl;
      }
    }
    acceleratorCopyToDevice(&HOST_C[0],&BLAS_C[0],BLAS_C.size()*sizeof(scalar));

    
    /////////////////////////////////////////
    // Guess G_xr = Exe Cer
    /////////////////////////////////////////
    BLAS.gemmBatched(GridBLAS_OP_N,GridBLAS_OP_N, 
		     vw,nrhs,nev,
		     scalar(1.0),
		     Ed, // x . nev
		     Cd, // nev . nrhs
		     scalar(0.0),
		     Gd);
    BLAS.synchronise();

    ///////////////////////////////////////
    // Copy out the multirhs
    ///////////////////////////////////////
    for(int r=0;r<nrhs;r++){
      int64_t offset = r*vw;
      autoView(v,guess[r],AcceleratorWrite);
      acceleratorCopyDeviceToDevice(&BLAS_G[offset],&v[0],sizeof(scalar_object)*vol);
    }
    RealD t1 = usecond();
    std::cout << GridLogMessage << "MultiRHSDeflation for "<<nrhs<<" sources with "<<nev<<" eigenvectors took " << (t1-t0)/1e3 <<" ms"<<std::endl;
  }
};

NAMESPACE_END(Grid);
