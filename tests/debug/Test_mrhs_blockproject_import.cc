/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/debug/Test_mrhs_blockproject_import.cc

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
// ImportFineGridMrhsVectors / ExportCoarseGridMrhsVectors against the vector
// of single RHS routines they replace. The BLAS buffers must come out
// identical, so the mrhs path is a pure removal of the unpack and repack.
//
#include <Grid/Grid.h>

using namespace Grid;

const int nbasis = 8;
const int nrhs   = 4;

typedef LatticeFermionD FineField;
typedef iVector<sTComplexD,nbasis> CoarseSiteObj;   // unvectorised coarse space
typedef Lattice<CoarseSiteObj>    CoarseField;

template<class T>
RealD BufDiff(deviceVector<T> &a,deviceVector<T> &b)
{
  GRID_ASSERT(a.size()==b.size());
  std::vector<T> ha(a.size()),hb(b.size());
  acceleratorCopyFromDevice(&a[0],&ha[0],a.size()*sizeof(T));
  acceleratorCopyFromDevice(&b[0],&hb[0],b.size()*sizeof(T));
  RealD num=0.0;
  for(int64_t i=0;i<a.size();i++){
    T d = ha[i]-hb[i];
    num += real(d)*real(d)+imag(d)*imag(d);
  }
  return num;
}
template<class T>
RealD BufNorm(deviceVector<T> &a)
{
  std::vector<T> ha(a.size());
  acceleratorCopyFromDevice(&a[0],&ha[0],a.size()*sizeof(T));
  RealD num=0.0;
  for(int64_t i=0;i<a.size();i++) num += real(ha[i])*real(ha[i])+imag(ha[i])*imag(ha[i]);
  return num;
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate flatt = GridDefaultLatt();
  Coordinate fsimd = GridDefaultSimd(Nd,vComplexD::Nsimd());
  Coordinate fmpi  = GridDefaultMpi();

  Coordinate block({2,2,2,2});
  Coordinate clatt(Nd);
  for(int d=0;d<Nd;d++) clatt[d]=flatt[d]/block[d];

  GridCartesian *FineGrid   = new GridCartesian(flatt,fsimd,fmpi);
  Coordinate csimd(Nd,1);                            // coarse space is unvectorised
  GridCartesian *CoarseGrid = new GridCartesian(clatt,csimd,fmpi);

  // D+1 grids, rhs innermost and unvectorised
  Coordinate fml(1,nrhs),fms(1,1),fmm(1,1);
  Coordinate cml(1,nrhs),cms(1,1),cmm(1,1);
  for(int d=0;d<Nd;d++){
    fml.push_back(flatt[d]); fms.push_back(fsimd[d]); fmm.push_back(fmpi[d]);
    cml.push_back(clatt[d]); cms.push_back(csimd[d]); cmm.push_back(fmpi[d]);
  }
  GridCartesian *FineMulti   = new GridCartesian(fml,fms,fmm);
  GridCartesian *CoarseMulti = new GridCartesian(cml,cms,cmm);

  std::cout << GridLogMessage << "fine "<<flatt<<"  coarse "<<clatt<<"  nrhs "<<nrhs
	    << "  fine Nsimd "<<FineGrid->Nsimd()<<"  coarse Nsimd "<<CoarseGrid->Nsimd()<<std::endl;

  GridParallelRNG pRNG(FineGrid); pRNG.SeedFixedIntegers(std::vector<int>({1,2,3,4}));

  MultiRHSBlockProject<FineField> Projector;
  Projector.Allocate(nbasis,FineGrid,CoarseGrid);

  ///////////////////////////////////////////////////////
  // Same data as nrhs single RHS fields and as one D+1
  ///////////////////////////////////////////////////////
  std::vector<FineField> vecs(nrhs,FineGrid);
  FineField mrhs(FineMulti);
  for(int r=0;r<nrhs;r++){
    random(pRNG,vecs[r]);
    InsertSliceFast(vecs[r],mrhs,r,0);
  }

  int64_t fsz = FineGrid->lSites()*Projector.words*nrhs;
  deviceVector<ComplexD> Fa(fsz), Fb(fsz);

  Projector.ImportFineGridVectors    (vecs,Fa);
  Projector.ImportFineGridMrhsVectors(mrhs,Fb);

  RealD fn = BufNorm(Fa);
  RealD fd = BufDiff(Fa,Fb);
  std::cout << GridLogMessage << "fine import  |vector|^2 = " << fn
	    << "   |vector - mrhs|^2 = " << fd << std::endl;
  GRID_ASSERT( fn > 0.0 );
  GRID_ASSERT( fd == 0.0 );

  ///////////////////////////////////////////////////////
  // Coarse export: fill the BLAS buffer, read it back
  // both ways, and compare the resulting fields
  ///////////////////////////////////////////////////////
  int64_t csz = CoarseGrid->lSites()*nbasis*nrhs;
  deviceVector<ComplexD> Cbuf(csz);
  {
    std::vector<ComplexD> host(csz);
    GridSerialRNG sRNG; sRNG.SeedFixedIntegers(std::vector<int>({9,8,7,6}));
    for(int64_t i=0;i<csz;i++){
      RealD re,im; random(sRNG,re); random(sRNG,im);
      host[i]=ComplexD(re-0.5,im-0.5);
    }
    acceleratorCopyToDevice(&host[0],&Cbuf[0],csz*sizeof(ComplexD));
  }

  std::vector<CoarseField> cvecs(nrhs,CoarseGrid);
  CoarseField cmrhs(CoarseMulti);

  Projector.ExportCoarseGridVectors    (cvecs,Cbuf);
  Projector.ExportCoarseGridMrhsVectors(cmrhs,Cbuf);

  RealD cn=0.0, cd=0.0;
  for(int r=0;r<nrhs;r++){
    CoarseField slice(CoarseGrid);
    ExtractSliceFast(slice,cmrhs,r,0);
    CoarseField e(CoarseGrid);
    e = slice - cvecs[r];
    cn += norm2(cvecs[r]);
    cd += norm2(e);
  }
  std::cout << GridLogMessage << "coarse export |vector|^2 = " << cn
	    << "   |vector - mrhs|^2 = " << cd << std::endl;
  GRID_ASSERT( cn > 0.0 );
  GRID_ASSERT( cd == 0.0 );

  ///////////////////////////////////////////////////////
  // All four blockProject orderings must agree
  ///////////////////////////////////////////////////////
  {
    std::vector<FineField> basis(nbasis,FineGrid);
    for(int b=0;b<nbasis;b++) random(pRNG,basis[b]);
    Projector.ImportBasis(basis);

    std::vector<CoarseField> cvv(nrhs,CoarseGrid);   // vector fine -> vector coarse
    std::vector<CoarseField> cmv(nrhs,CoarseGrid);   // mrhs   fine -> vector coarse
    CoarseField cvm(CoarseMulti);                    // vector fine -> mrhs   coarse
    CoarseField cmm(CoarseMulti);                    // mrhs   fine -> mrhs   coarse

    Projector.blockProject(vecs,cvv);
    Projector.blockProject(mrhs,cmv);
    Projector.blockProject(vecs,cvm);
    Projector.blockProject(mrhs,cmm);

    RealD ref=0.0, d_mv=0.0, d_vm=0.0, d_mm=0.0;
    for(int r=0;r<nrhs;r++){
      CoarseField svm(CoarseGrid),smm(CoarseGrid),e(CoarseGrid);
      ExtractSliceFast(svm,cvm,r,0);
      ExtractSliceFast(smm,cmm,r,0);
      ref  += norm2(cvv[r]);
      e = cmv[r]-cvv[r];  d_mv += norm2(e);
      e = svm    -cvv[r]; d_vm += norm2(e);
      e = smm    -cvv[r]; d_mm += norm2(e);
    }
    std::cout << GridLogMessage << "blockProject |vec->vec|^2 = " << ref << std::endl;
    std::cout << GridLogMessage << "   mrhs->vec diff " << d_mv
	      << "   vec->mrhs diff " << d_vm
	      << "   mrhs->mrhs diff " << d_mm << std::endl;
    GRID_ASSERT( ref > 0.0 );
    GRID_ASSERT( d_mv == 0.0 );
    GRID_ASSERT( d_vm == 0.0 );
    GRID_ASSERT( d_mm == 0.0 );

    ///////////////////////////////////////////////////////
    // Promote must agree across orderings too
    ///////////////////////////////////////////////////////
    std::vector<FineField> fvv(nrhs,FineGrid);
    FineField fmm(FineMulti);
    Projector.blockPromote(fvv,cvv);
    Projector.blockPromote(fmm,cmm);

    RealD pref=0.0, pdiff=0.0;
    for(int r=0;r<nrhs;r++){
      FineField s(FineGrid),e(FineGrid);
      ExtractSliceFast(s,fmm,r,0);
      e = s - fvv[r];
      pref  += norm2(fvv[r]);
      pdiff += norm2(e);
    }
    std::cout << GridLogMessage << "blockPromote |vec|^2 = " << pref
	      << "   mrhs diff " << pdiff << std::endl;
    GRID_ASSERT( pref > 0.0 );
    GRID_ASSERT( pdiff == 0.0 );
  }

  std::cout << GridLogMessage << "Test_mrhs_blockproject_import: ALL PASS" << std::endl;

  Grid_finalize();
}
