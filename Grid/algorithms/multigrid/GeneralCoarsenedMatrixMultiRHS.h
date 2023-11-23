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
    for(int p=0;p<geom.npoint;p++){
      _A[p].resize(_CoarseGrid->lSites());
    }
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
    conformable(CoarseGrid(),in.Grid());
    conformable(in.Grid(),out.Grid());
    out.Checkerboard() = in.Checkerboard();
    CoarseVector tin=in;

    CoarseVector pin = Cell.ExchangePeriodic(tin);
    CoarseVector pout(pin.Grid());

    int npoint = geom.npoint;
    typedef calcMatrix* Aview;
    typedef LatticeView<Cvec> Vview;
      
    const int Nsimd = CComplex::Nsimd();
    
    int64_t osites=pin.Grid()->oSites();
    int64_t nrhs  =pin.Grid()->GlobalDimensions()[0]/Nsimd;

    {
      autoView( in_v , pin, AcceleratorRead);
      autoView( out_v , pout, AcceleratorWriteDiscard);
      autoView( Stencil_v  , Stencil, AcceleratorRead);

      // Static and prereserve to keep UVM region live and not resized across multiple calls
      MultTemporaries.resize(npoint,pin.Grid());       

      std::vector<Aview> AcceleratorViewContainer_h;
      std::vector<Vview> AcceleratorVecViewContainer_h; 

      for(int p=0;p<npoint;p++) {
	AcceleratorViewContainer_h.push_back( &_A[p][0]);
	AcceleratorVecViewContainer_h.push_back(MultTemporaries[p].View(AcceleratorWrite));
      }

      static deviceVector<Aview> AcceleratorViewContainer; AcceleratorViewContainer.resize(npoint);
      static deviceVector<Vview> AcceleratorVecViewContainer; AcceleratorVecViewContainer.resize(npoint); 
      
      auto Aview_p = &AcceleratorViewContainer[0];
      auto Vview_p = &AcceleratorVecViewContainer[0];

      acceleratorCopyToDevice(&AcceleratorViewContainer_h[0],&AcceleratorViewContainer[0],npoint *sizeof(Aview));
      acceleratorCopyToDevice(&AcceleratorVecViewContainer_h[0],&AcceleratorVecViewContainer[0],npoint *sizeof(Vview));

      accelerator_for(rspb, osites*nbasis*npoint, Nsimd, {
	  typedef decltype(coalescedRead(in_v[0](0))) calcComplex;
	  int32_t ss   = rspb/(nbasis*npoint);
	  int32_t bp   = rspb%(nbasis*npoint);
	  int32_t point= bp/nbasis;
	  int32_t b    = bp%nbasis;
	  auto SE  = Stencil_v.GetEntry(point,ss);
	  int32_t snbr= SE->_offset;
	  auto nbr = coalescedReadGeneralPermute(in_v[snbr],SE->_permute,Nd);
	  auto res = Aview_p[point][ss](0,b)*nbr(0);
	  for(int bb=1;bb<nbasis;bb++) {
	    res = res + Aview_p[point][ss](bb,b)*nbr(bb);
	  }
	  coalescedWrite(Vview_p[point][ss](b),res);
      });
      accelerator_for(sb, osites*nbasis, Nsimd, {
	  int ss = sb/nbasis;
	  int b  = sb%nbasis;
	  auto res = coalescedRead(Vview_p[0][ss](b));
	  for(int point=1;point<npoint;point++){
	    res = res + coalescedRead(Vview_p[point][ss](b));
	  }
	  coalescedWrite(out_v[ss](b),res);
      });
      for(int p=0;p<npoint;p++) {
	AcceleratorVecViewContainer_h[p].ViewClose();
      }
    }

    out = Cell.Extract(pout);

  };
  virtual  void Mdiag    (const Field &in, Field &out){ assert(0);};
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp){assert(0);};
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out){assert(0);};

};
  
NAMESPACE_END(Grid);
