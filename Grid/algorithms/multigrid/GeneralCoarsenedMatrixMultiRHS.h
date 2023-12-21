/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/GeneralCoarsenedMatrixMultiRHS.h

    Copyright (C) 2015

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

#include <Grid/algorithms/multigrid/BatchedBlas.h>

NAMESPACE_BEGIN(Grid);


// Move this to accelerator.h
// Also give a copy device.
// Rename acceleratorPut
// Rename acceleratorGet
template<class T> void deviceSet(T& dev,T&host)
{
  acceleratorCopyToDevice(&host,&dev,sizeof(T));
}
template<class T> T deviceGet(T& dev)
{
  T host;
  acceleratorCopyFromDevice(&dev,&host,sizeof(T));
  return host;
}

// Fine Object == (per site) type of fine field
// nbasis      == number of deflation vectors
template<class Fobj,class CComplex,int nbasis>
class MultiGeneralCoarsenedMatrix : public SparseMatrixBase<Lattice<iVector<CComplex,nbasis > > >  {
public:
  typedef typename CComplex::scalar_object SComplex;
  typedef GeneralCoarsenedMatrix<Fobj,CComplex,nbasis> GeneralCoarseOp;
  typedef MultiGeneralCoarsenedMatrix<Fobj,CComplex,nbasis> MultiGeneralCoarseOp;

  typedef iVector<CComplex,nbasis >           siteVector;
  typedef iMatrix<CComplex,nbasis >           siteMatrix;
  typedef iVector<SComplex,nbasis >           calcVector;
  typedef iMatrix<SComplex,nbasis >           calcMatrix;
  typedef Lattice<iScalar<CComplex> >         CoarseComplexField;
  typedef Lattice<siteVector>                 CoarseVector;
  typedef Lattice<iMatrix<CComplex,nbasis > > CoarseMatrix;
  typedef iMatrix<CComplex,nbasis >  Cobj;
  typedef iVector<CComplex,nbasis >  Cvec;
  typedef Lattice< CComplex >   CoarseScalar; // used for inner products on fine field
  typedef Lattice<Fobj >        FineField;
  typedef CoarseVector Field;

  ////////////////////
  // Data members
  ////////////////////
  GridCartesian *       _CoarseGridMulti; 
  GridCartesian *       _CoarseGrid;
  GeneralCoarseOp &     _Op;
  NonLocalStencilGeometry geom;
  PaddedCell Cell;
  GeneralLocalStencil Stencil;

  deviceVector<calcVector> BLAS_B;
  deviceVector<calcVector> BLAS_C;
  std::vector<deviceVector<calcMatrix> > BLAS_A;

  std::vector<deviceVector<ComplexD *> > BLAS_AP;
  std::vector<deviceVector<ComplexD *> > BLAS_BP;
  deviceVector<ComplexD *>               BLAS_CP;

  ///////////////////////
  // Interface
  ///////////////////////
  GridBase      * Grid(void)           { return _CoarseGridMulti; };   // this is all the linalg routines need to know
  GridCartesian * CoarseGrid(void)     { return _CoarseGridMulti; };   // this is all the linalg routines need to know

  MultiGeneralCoarsenedMatrix(GeneralCoarseOp & Op,GridCartesian *CoarseGridMulti) :
    _Op(Op),
    _CoarseGrid(Op.CoarseGrid()),
    _CoarseGridMulti(CoarseGridMulti),
    geom(_CoarseGridMulti,Op.geom.hops,Op.geom.skip+1),
    Cell(Op.geom.Depth(),_CoarseGridMulti),
    Stencil(Cell.grids.back(),geom.shifts) // padded cell stencil
  {
    int32_t padded_sites   = _Op._A[0].Grid()->lSites();
    int32_t unpadded_sites = _CoarseGrid->lSites();
    
    int32_t nrhs  = CoarseGridMulti->FullDimensions()[0];  // # RHS
    int32_t orhs  = nrhs/CComplex::Nsimd();

    /////////////////////////////////////////////////
    // Device data vector storage
    /////////////////////////////////////////////////
    BLAS_A.resize(geom.npoint);
    for(int p=0;p<geom.npoint;p++){
      BLAS_A[p].resize (unpadded_sites); // no ghost zone, npoint elements
    }
    BLAS_B.resize(nrhs *padded_sites);   // includes ghost zone
    BLAS_C.resize(nrhs *unpadded_sites); // no ghost zone

    BLAS_AP.resize(geom.npoint);
    BLAS_BP.resize(geom.npoint);
    for(int p=0;p<geom.npoint;p++){
      BLAS_AP[p].resize(unpadded_sites);
      BLAS_BP[p].resize(unpadded_sites);
    }
    BLAS_CP.resize(unpadded_sites);

    /////////////////////////////////////////////////
    // Pointers to data
    /////////////////////////////////////////////////

    // Site identity mapping for A, C
    for(int p=0;p<geom.npoint;p++){
      for(int ss=0;ss<unpadded_sites;ss++){
	ComplexD *ptr = (ComplexD *)&BLAS_A[p][ss];
	//ComplexD *ptr = (ComplexD *)&BLAS_A[p][0]; std::cout << " A ptr "<<std::hex<<ptr<<std::dec<<" "<<ss<<"/"<<BLAS_A[p].size()<<std::endl;
	deviceSet(BLAS_AP[p][ss],ptr);
      }
    }
    for(int ss=0;ss<unpadded_sites;ss++){
      ComplexD *ptr = (ComplexD *)&BLAS_C[ss*nrhs];
      //ComplexD *ptr = (ComplexD *)&BLAS_C[0];  std::cout << " C ptr "<<std::hex<<ptr<<std::dec<<" "<<ss<<"/"<<BLAS_C.size()<<std::endl;
      deviceSet(BLAS_CP[ss],ptr);
    }

    /////////////////////////////////////////////////
    // Neighbour table is more complicated
    /////////////////////////////////////////////////
    int32_t j=0; // Interior point counter (unpadded)
    for(int32_t s=0;s<padded_sites;s++){ // 4 volume, padded
      int ghost_zone=0;
      for(int32_t point = 0 ; point < geom.npoint; point++){
	int i=s*orhs*geom.npoint+point;
	if( Stencil._entries[i]._wrap ) { // stencil is indexed by the oSite of the CoarseGridMulti, hence orhs factor
	  ghost_zone=1; // If general stencil wrapped in any direction, wrap=1
	}
      }
      //      GeneralStencilEntryReordered tmp;
      if( ghost_zone==0) {
	for(int32_t point = 0 ; point < geom.npoint; point++){
	  int i=s*orhs*geom.npoint+point;
 	  int32_t nbr = Stencil._entries[i]._offset*CComplex::Nsimd(); // oSite -> lSite
	  //	  std::cout << " B ptr "<< nbr<<"/"<<BLAS_B.size()<<std::endl;
	  assert(nbr<BLAS_B.size());
	  ComplexD * ptr = (ComplexD *)&BLAS_B[nbr];
	  //	  ComplexD * ptr = (ComplexD *)&BLAS_B[0];
	  //	  std::cout << " B ptr unpadded "<<std::hex<<ptr<<std::dec<<" "<<s<<"/"<<padded_sites<<std::endl;
	  //	  std::cout << " B ptr   padded "<<std::hex<<ptr<<std::dec<<" "<<j<<"/"<<unpadded_sites<<std::endl;
	  deviceSet(BLAS_BP[point][j],ptr); // neighbour indexing in ghost zone volume
	  //	  auto tmp = deviceGet(*BLAS_BP[point][j]);  // debug trigger SEGV if bad ptr
	}
	j++;
      }
    }
    assert(j==unpadded_sites);
    CopyMatrix();
  }
  template<class vobj> void GridtoBLAS(const Lattice<vobj> &grid,deviceVector<typename vobj::scalar_object> &out)
  {
    std::vector<typename vobj::scalar_object> tmp;
    unvectorizeToLexOrdArray(tmp,grid);
    //    std::cout << "GridtoBLAS volume " <<tmp.size()<<" " << grid.Grid()->lSites()<<" "<<out.size()<<std::endl;
    //    std::cout << "GridtoBLAS site 0 " <<tmp[0]<<std::endl;
    assert(tmp.size()==grid.Grid()->lSites());
    assert(tmp.size()==out.size());
    out.resize(tmp.size());
    acceleratorCopyToDevice(&tmp[0],&out[0],sizeof(typename vobj::scalar_object)*tmp.size());
  }    
  template<class vobj> void BLAStoGrid(Lattice<vobj> &grid,deviceVector<typename vobj::scalar_object> &in)
  {
    std::vector<typename vobj::scalar_object> tmp;
    tmp.resize(in.size());
    //    std::cout << "BLAStoGrid volume " <<tmp.size()<<" "<< grid.Grid()->lSites()<<std::endl;
    assert(in.size()==grid.Grid()->lSites());
    acceleratorCopyFromDevice(&in[0],&tmp[0],sizeof(typename vobj::scalar_object)*in.size());
    vectorizeFromLexOrdArray(tmp,grid);
  }
  void CopyMatrix (void)
  {
    // Clone "A" to be lexicographic in the physics coords
    // Use unvectorisetolexordarray
    // Copy to device
    for(int p=0;p<geom.npoint;p++){
      //Unpadded
      auto Aup = _Op.Cell.Extract(_Op._A[p]);
      //      Coordinate coor({0,0,0,0,0});
      //      auto sval = peekSite(Aup,coor);
      //      std::cout << "CopyMatrix: p "<<p<<" Aup[0] :"<<sval<<std::endl;
      //      sval = peekSite(_Op._A[p],coor);
      //      std::cout << "CopyMatrix: p "<<p<<" _Op._Ap[0] :"<<sval<<std::endl;
      GridtoBLAS(Aup,BLAS_A[p]);
      //      std::cout << "Copy Matrix p "<<p<<" "<< deviceGet(BLAS_A[p][0])<<std::endl;
    }
  }
  void Mdag(const CoarseVector &in, CoarseVector &out)
  {
    this->M(in,out);
  }
  void M (const CoarseVector &in, CoarseVector &out)
  {
    std::cout << GridLogMessage << "New Mrhs coarse"<<std::endl;
    conformable(CoarseGrid(),in.Grid());
    conformable(in.Grid(),out.Grid());
    out.Checkerboard() = in.Checkerboard();

    RealD t_tot;
    RealD t_exch;
    RealD t_GtoB;
    RealD t_BtoG;
    RealD t_mult;

    t_tot=-usecond();
    CoarseVector tin=in;
    t_exch=-usecond();
    CoarseVector pin = Cell.ExchangePeriodic(tin); //padded input
    t_exch+=usecond();

    CoarseVector pout(pin.Grid());

    int npoint = geom.npoint;
    typedef calcMatrix* Aview;
    typedef LatticeView<Cvec> Vview;
      
    const int Nsimd = CComplex::Nsimd();

    RealD flops,bytes;
    int64_t osites=in.Grid()->oSites(); // unpadded
    int64_t unpadded_vol = _CoarseGrid->lSites();
    
    flops = 1.0* npoint * nbasis * nbasis * 8.0 * osites * CComplex::Nsimd();
    bytes = 1.0*osites*sizeof(siteMatrix)*npoint/pin.Grid()->GlobalDimensions()[0]
          + 2.0*osites*sizeof(siteVector)*npoint;
    
    int64_t nrhs  =pin.Grid()->GlobalDimensions()[0];
    assert(nrhs>=1);

    std::cout << GridLogMessage << "New Mrhs GridtoBLAS in sizes "<<in.Grid()->lSites()<<" "<<pin.Grid()->lSites()<<std::endl;
    t_GtoB=-usecond();
    GridtoBLAS(pin,BLAS_B);
    //    out = Zero();
    //    GridtoBLAS(out,BLAS_C);
    t_GtoB+=usecond();

    GridBLAS BLAS;

    t_mult=-usecond();
    for(int p=0;p<geom.npoint;p++){
      RealD c = 1.0;
      if (p==0) c = 0.0;
      ComplexD beta(c);
      //      std::cout << GridLogMessage << "New Mrhs coarse gemmBatched "<<p<<std::endl;
      BLAS.gemmBatched(nbasis,nrhs,nbasis,
		       ComplexD(1.0),
		       BLAS_AP[p], 
		       BLAS_BP[p], 
		       ComplexD(c), 
		       BLAS_CP);
    }
    t_mult+=usecond();
    //    std::cout << GridLogMessage << "New Mrhs coarse BLAStoGrid "<<std::endl;
    t_BtoG=-usecond();
    BLAStoGrid(out,BLAS_C);
    t_BtoG+=usecond();
    t_tot+=usecond();
    //    auto check =deviceGet(BLAS_C[0]);
    //      std::cout << "C[0] "<<check<<std::endl;
    //    Coordinate coor({0,0,0,0,0,0});
    //    peekLocalSite(check,out,coor);
    //    std::cout << "C[0] "<< check<<std::endl;
    std::cout << GridLogMessage << "New Mrhs coarse DONE "<<std::endl;
    std::cout << GridLogMessage<<"Coarse Mult exch "<<t_exch<<" us"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Mult mult "<<t_mult<<" us"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Mult GtoB  "<<t_GtoB<<" us"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Mult BtoG  "<<t_BtoG<<" us"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Mult tot  "<<t_tot<<" us"<<std::endl;
    std::cout << GridLogMessage<<std::endl;
    std::cout << GridLogMessage<<"Coarse Kernel flops "<< flops<<std::endl;
    std::cout << GridLogMessage<<"Coarse Kernel flop/s "<< flops/t_mult<<" mflop/s"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Kernel bytes/s "<< bytes/t_mult/1000<<" GB/s"<<std::endl;
    std::cout << GridLogMessage<<"Coarse overall flops/s "<< flops/t_tot<<" mflop/s"<<std::endl;
    std::cout << GridLogMessage<<"Coarse total bytes   "<< bytes/1e6<<" MB"<<std::endl;
  };
  virtual  void Mdiag    (const Field &in, Field &out){ assert(0);};
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp){assert(0);};
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out){assert(0);};

};
  
NAMESPACE_END(Grid);
