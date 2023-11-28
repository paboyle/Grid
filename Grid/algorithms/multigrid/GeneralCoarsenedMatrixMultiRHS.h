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

NAMESPACE_BEGIN(Grid);

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
  
  std::vector<deviceVector<calcMatrix> > _A;
  std::vector<CoarseVector> MultTemporaries;
  deviceVector<GeneralStencilEntryReordered> StencilMasked;

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
    Stencil(Cell.grids.back(),geom.shifts)
  {
    _A.resize(geom.npoint);
    int32_t padded_sites = _Op._A[0].Grid()->lSites();
    for(int p=0;p<geom.npoint;p++){
      _A[p].resize(padded_sites);
    }
    std::cout << GridLogMessage<<"MultiGeneralCoarsenedMatrix "<<_CoarseGrid->lSites()<<" coarse sites "<<_Op._A[0].Grid()->lSites() <<std::endl;

    StencilMasked.resize(_CoarseGridMulti->oSites()*geom.npoint);
    std::vector<GeneralStencilEntryReordered> StencilTmp;

    int32_t j=0;
    int32_t sites = Stencil._entries.size()/geom.npoint;
    for(int32_t s=0;s<sites;s++){
      int ghost_zone=0;
      for(int32_t point = 0 ; point < geom.npoint; point++){
	int i=s*geom.npoint+point;
	if( Stencil._entries[i]._permute ) {
	  ghost_zone=1;
	}
      }
      GeneralStencilEntryReordered tmp;
      if( ghost_zone==0) {
	for(int32_t point = 0 ; point < geom.npoint; point++){
	  int i=s*geom.npoint+point;
 	  tmp._offset = Stencil._entries[i]._offset;
	  tmp._permute= Stencil._entries[i]._permute; // Should be no premute and j=site
	  tmp._input = s;
	  StencilTmp.push_back(tmp);
	}
	j++;
      }
    }

    std::cout << " oSites " << _CoarseGridMulti->oSites()<<std::endl;
    std::cout << " npoint " << geom.npoint<<std::endl;
    std::cout << " StencilTmp "<<StencilTmp.size();
    assert(_CoarseGridMulti->oSites()*geom.npoint==StencilTmp.size());
    acceleratorCopyToDevice(&StencilTmp[0],&StencilMasked[0],sizeof(GeneralStencilEntryReordered)*StencilTmp.size());
    CopyMatrix();
  }
  void CopyMatrix (void)
  {
    // Clone "A" to be lexicographic in the physics coords
    // Use unvectorisetolexordarray
    // Copy to device
    std::vector<calcMatrix> tmp;
    for(int p=0;p<geom.npoint;p++){
      unvectorizeToLexOrdArray(tmp,_Op._A[p]);
      acceleratorCopyToDevice(&tmp[0],&_A[p][0],sizeof(calcMatrix)*tmp.size());
    }
  }
  void Mdag(const CoarseVector &in, CoarseVector &out)
  {
    this->M(in,out);
  }
  void M (const CoarseVector &in, CoarseVector &out)
  {
    RealD tviews=0;    RealD ttot=0;    RealD tmult=0;   RealD texch=0;    RealD text=0; RealD ttemps=0; RealD tcopy=0;
    RealD tmult2=0;

    ttot=-usecond();
    conformable(CoarseGrid(),in.Grid());
    conformable(in.Grid(),out.Grid());
    out.Checkerboard() = in.Checkerboard();
    CoarseVector tin=in;

    texch-=usecond();
    CoarseVector pin = Cell.ExchangePeriodic(tin);
    texch+=usecond();
    CoarseVector pout(pin.Grid());

    int npoint = geom.npoint;
    typedef calcMatrix* Aview;
    typedef LatticeView<Cvec> Vview;
      
    const int Nsimd = CComplex::Nsimd();

    RealD flops,bytes;
    int64_t nrhs  =pin.Grid()->GlobalDimensions()[0]/Nsimd;
    assert(nrhs>=1);

#if 0
    {
      tviews-=usecond();
      autoView( in_v , pin, AcceleratorRead);
      autoView( out_v , pout, AcceleratorWriteDiscard);
      tviews+=usecond();

      // Static and prereserve to keep UVM region live and not resized across multiple calls
      ttemps-=usecond();
      MultTemporaries.resize(npoint,in.Grid());       
      ttemps+=usecond();

      std::vector<Aview> AcceleratorViewContainer_h;
      std::vector<Vview> AcceleratorVecViewContainer_h; 

      tviews-=usecond();
      for(int p=0;p<npoint;p++) {
	AcceleratorViewContainer_h.push_back( &_A[p][0]);
	AcceleratorVecViewContainer_h.push_back(MultTemporaries[p].View(AcceleratorWrite));
      }
      tviews+=usecond();

      static deviceVector<Aview> AcceleratorViewContainer; AcceleratorViewContainer.resize(npoint);
      static deviceVector<Vview> AcceleratorVecViewContainer; AcceleratorVecViewContainer.resize(npoint); 
      
      auto Aview_p = &AcceleratorViewContainer[0];
      auto Vview_p = &AcceleratorVecViewContainer[0];

      tcopy-=usecond();
      acceleratorCopyToDevice(&AcceleratorViewContainer_h[0],&AcceleratorViewContainer[0],npoint *sizeof(Aview));
      acceleratorCopyToDevice(&AcceleratorVecViewContainer_h[0],&AcceleratorVecViewContainer[0],npoint *sizeof(Vview));
      tcopy+=usecond();

      int32_t bound = _A[0].size();
      int64_t osites=pin.Grid()->oSites();
      flops = 1.0* npoint * nbasis * nbasis * 8.0 * osites * CComplex::Nsimd();
      bytes = 1.0*osites*sizeof(siteMatrix)*npoint/pin.Grid()->GlobalDimensions()[0]
	+ 2.0*osites*sizeof(siteVector)*npoint;

      std::cout << " osites "<<osites <<" bound "<<bound<<std::endl;
      std::cout << " padded local dims   "<<pin.Grid()->LocalDimensions()<<std::endl;
      std::cout << " unpadded local dims "<<in.Grid()->LocalDimensions()<<std::endl;
      tmult-=usecond();
      autoView( Stencil_v  , Stencil, AcceleratorRead);
      accelerator_for(rspb, osites*nbasis*npoint, Nsimd, {
	  typedef decltype(coalescedRead(in_v[0](0))) calcComplex;
	  int32_t ss   = rspb/(nbasis*npoint);
	  int32_t bp   = rspb%(nbasis*npoint);
	  int32_t point= bp/nbasis;
	  int32_t b    = bp%nbasis;
	  assert(ss<bound);
	  auto SE  = Stencil_v.GetEntry(point,ss);
	  if ( SE->_permute == 0 ) {
	    int32_t snbr= SE->_offset;
	    auto nbr = coalescedReadGeneralPermute(in_v[snbr],SE->_permute,Nd);
	    auto res = Aview_p[point][ss](0,b)*nbr(0);
	    for(int bb=1;bb<nbasis;bb++) {
	      res = res + Aview_p[point][ss](bb,b)*nbr(bb);
	    }
	    coalescedWrite(Vview_p[point][ss](b),res);
	  }
      });
      tmult2-=usecond();
      accelerator_for(sb, osites*nbasis, Nsimd, {
	  int ss = sb/nbasis;
	  int b  = sb%nbasis;
	  auto res = coalescedRead(Vview_p[0][ss](b));
	  for(int point=1;point<npoint;point++){
	    res = res + coalescedRead(Vview_p[point][ss](b));
	  }
	  coalescedWrite(out_v[ss](b),res);
      });
      tmult2+=usecond();
      tmult+=usecond();

      for(int p=0;p<npoint;p++) {
	AcceleratorVecViewContainer_h[p].ViewClose();
      }
    }

    text-=usecond();
    out = Cell.Extract(pout);
    text+=usecond();
    ttot+=usecond();
#else
    {
      tviews-=usecond();
      autoView( in_v , pin, AcceleratorRead);
      autoView( out_v , out, AcceleratorWriteDiscard);
      tviews+=usecond();

      // Static and prereserve to keep UVM region live and not resized across multiple calls
      ttemps-=usecond();
      MultTemporaries.resize(npoint,in.Grid());       
      ttemps+=usecond();

      std::vector<Aview> AcceleratorViewContainer_h;
      std::vector<Vview> AcceleratorVecViewContainer_h; 

      tviews-=usecond();
      for(int p=0;p<npoint;p++) {
	AcceleratorViewContainer_h.push_back( &_A[p][0]);
	AcceleratorVecViewContainer_h.push_back(MultTemporaries[p].View(AcceleratorWrite));
      }
      tviews+=usecond();

      static deviceVector<Aview> AcceleratorViewContainer; AcceleratorViewContainer.resize(npoint);
      static deviceVector<Vview> AcceleratorVecViewContainer; AcceleratorVecViewContainer.resize(npoint); 
      
      auto Aview_p = &AcceleratorViewContainer[0];
      auto Vview_p = &AcceleratorVecViewContainer[0];

      tcopy-=usecond();
      acceleratorCopyToDevice(&AcceleratorViewContainer_h[0],&AcceleratorViewContainer[0],npoint *sizeof(Aview));
      acceleratorCopyToDevice(&AcceleratorVecViewContainer_h[0],&AcceleratorVecViewContainer[0],npoint *sizeof(Vview));
      tcopy+=usecond();

      int32_t bound = _A[0].size();
      int64_t osites=in.Grid()->oSites();
      flops = 1.0* npoint * nbasis * nbasis * 8.0 * osites * CComplex::Nsimd();
      bytes = 1.0*osites*sizeof(siteMatrix)*npoint/pin.Grid()->GlobalDimensions()[0]
	+ 2.0*osites*sizeof(siteVector)*npoint;

      std::cout << " osites "<<osites <<" bound "<<bound<< " stencilsize  "<<StencilMasked.size()<<std::endl;
      std::cout << " padded local dims   "<<pin.Grid()->LocalDimensions()<<std::endl;
      std::cout << " unpadded local dims "<<in.Grid()->LocalDimensions()<<std::endl;
      tmult-=usecond();
      auto Stencil_v = &StencilMasked[0];
      accelerator_for(rspb, StencilMasked.size()*nbasis, Nsimd, {
	  typedef decltype(coalescedRead(in_v[0](0))) calcComplex;
	  int32_t ss   = rspb/(nbasis*npoint); // site of unpadded
	  int32_t bp   = rspb%(nbasis*npoint);
	  int32_t point= bp/nbasis;
	  int32_t b    = bp%nbasis;
	  auto SE  = &Stencil_v[ss*npoint+point];
	  int32_t s   = SE->_input;
	  int32_t snbr= SE->_offset;
	  std::cout << " unpadded " << ss<<" padded " << s<< " point "<<point <<" row " <<b<<std::endl;
	  auto nbr = coalescedRead(in_v[snbr]);
	  auto res = Aview_p[point][s](0,b)*nbr(0);
	  for(int bb=1;bb<nbasis;bb++) {
	    res = res + Aview_p[point][s](bb,b)*nbr(bb);
	  }
	  coalescedWrite(Vview_p[point][ss](b),res);
      });
      tmult2-=usecond();
      accelerator_for(sb, osites*nbasis, Nsimd, {
	  int ss = sb/nbasis;
	  int b  = sb%nbasis;
	  auto res = coalescedRead(Vview_p[0][ss](b));
	  for(int point=1;point<npoint;point++){
	    res = res + coalescedRead(Vview_p[point][ss](b));
	  }
	  coalescedWrite(out_v[ss](b),res);
      });
      tmult2+=usecond();
      tmult+=usecond();
      for(int p=0;p<npoint;p++) {
	AcceleratorVecViewContainer_h[p].ViewClose();
      }
    }
    ttot+=usecond();
#endif

    std::cout << GridLogMessage<<"Coarse Mult Aviews "<<tviews<<" us"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Mult exch "<<texch<<" us"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Mult mult "<<tmult<<" us"<<std::endl;
    std::cout << GridLogMessage<<" of which mult2  "<<tmult2<<" us"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Mult ext  "<<text<<" us"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Mult temps "<<ttemps<<" us"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Mult copy  "<<tcopy<<" us"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Mult tot  "<<ttot<<" us"<<std::endl;
    //    std::cout << GridLogMessage<<std::endl;
    std::cout << GridLogMessage<<"Coarse Kernel flop/s "<< flops/tmult<<" mflop/s"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Kernel bytes/s"<< bytes/tmult<<" MB/s"<<std::endl;
    std::cout << GridLogMessage<<"Coarse overall flops/s "<< flops/ttot<<" mflop/s"<<std::endl;
    std::cout << GridLogMessage<<"Coarse total bytes   "<< bytes/1e6<<" MB"<<std::endl;
    
  };
  virtual  void Mdiag    (const Field &in, Field &out){ assert(0);};
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp){assert(0);};
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out){assert(0);};

};
  
NAMESPACE_END(Grid);
