/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/CoarsenedMatrix.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
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
#ifndef  GRID_ALGORITHM_COARSENED_MATRIX_H
#define  GRID_ALGORITHM_COARSENED_MATRIX_H

#include <Grid/qcd/QCD.h> // needed for Dagger(Yes|No), Inverse(Yes|No)

NAMESPACE_BEGIN(Grid);

template<class vobj,class CComplex>
inline void blockMaskedInnerProduct(Lattice<CComplex> &CoarseInner,
				    const Lattice<decltype(innerProduct(vobj(),vobj()))> &FineMask,
				    const Lattice<vobj> &fineX,
				    const Lattice<vobj> &fineY)
{
  typedef decltype(innerProduct(vobj(),vobj())) dotp;

  GridBase *coarse(CoarseInner.Grid());
  GridBase *fine  (fineX.Grid());

  Lattice<dotp> fine_inner(fine); fine_inner.Checkerboard() = fineX.Checkerboard();
  Lattice<dotp> fine_inner_msk(fine);

  // Multiply could be fused with innerProduct
  // Single block sum kernel could do both masks.
  fine_inner = localInnerProduct(fineX,fineY);
  mult(fine_inner_msk, fine_inner,FineMask);
  blockSum(CoarseInner,fine_inner_msk);
}

// Fine Object == (per site) type of fine field
// nbasis      == number of deflation vectors
template<class Fobj,class CComplex,int nbasis>
class CoarsenedMatrix : public CheckerBoardedSparseMatrixBase<Lattice<iVector<CComplex,nbasis > > >  {
public:
    
  typedef iVector<CComplex,nbasis >           siteVector;
  typedef Lattice<CComplex >                  CoarseComplexField;
  typedef Lattice<siteVector>                 CoarseVector;
  typedef Lattice<iMatrix<CComplex,nbasis > > CoarseMatrix;
  typedef iMatrix<CComplex,nbasis >  Cobj;
  typedef Lattice< CComplex >   CoarseScalar; // used for inner products on fine field
  typedef Lattice<Fobj >        FineField;
  typedef CoarseVector FermionField;

  // enrich interface, use default implementation as in FermionOperator ///////
  void Dminus(CoarseVector const& in, CoarseVector& out) { out = in; }
  void DminusDag(CoarseVector const& in, CoarseVector& out) { out = in; }
  void ImportPhysicalFermionSource(CoarseVector const& input, CoarseVector& imported) { imported = input; }
  void ImportUnphysicalFermion(CoarseVector const& input, CoarseVector& imported) { imported = input; }
  void ExportPhysicalFermionSolution(CoarseVector const& solution, CoarseVector& exported) { exported = solution; };
  void ExportPhysicalFermionSource(CoarseVector const& solution, CoarseVector& exported) { exported = solution; };

  ////////////////////
  // Data members
  ////////////////////
  Geometry         geom;
  GridBase *       _grid; 
  GridBase*        _cbgrid;
  int hermitian;

  CartesianStencil<siteVector,siteVector,DefaultImplParams> Stencil; 
  CartesianStencil<siteVector,siteVector,DefaultImplParams> StencilEven;
  CartesianStencil<siteVector,siteVector,DefaultImplParams> StencilOdd;

  std::vector<CoarseMatrix> A;
  std::vector<CoarseMatrix> Aeven;
  std::vector<CoarseMatrix> Aodd;

  CoarseMatrix AselfInv;
  CoarseMatrix AselfInvEven;
  CoarseMatrix AselfInvOdd;

  Vector<RealD> dag_factor;

  ///////////////////////
  // Interface
  ///////////////////////
  GridBase * Grid(void)         { return _grid; };   // this is all the linalg routines need to know
  GridBase * RedBlackGrid()     { return _cbgrid; };

  int ConstEE() { return 0; }

  void M (const CoarseVector &in, CoarseVector &out)
  {
    conformable(_grid,in.Grid());
    conformable(in.Grid(),out.Grid());
    out.Checkerboard() = in.Checkerboard();

    SimpleCompressor<siteVector> compressor;

    Stencil.HaloExchange(in,compressor);
    autoView( in_v , in, AcceleratorRead);
    autoView( out_v , out, AcceleratorWrite);
    autoView( Stencil_v  , Stencil, AcceleratorRead);
    int npoint = geom.npoint;
    typedef LatticeView<Cobj> Aview;
      
    Vector<Aview> AcceleratorViewContainer;
  
    for(int p=0;p<geom.npoint;p++) AcceleratorViewContainer.push_back(A[p].View(AcceleratorRead));
    Aview *Aview_p = & AcceleratorViewContainer[0];

    const int Nsimd = CComplex::Nsimd();
    typedef decltype(coalescedRead(in_v[0])) calcVector;
    typedef decltype(coalescedRead(in_v[0](0))) calcComplex;

    int osites=Grid()->oSites();

    accelerator_for(sss, Grid()->oSites()*nbasis, Nsimd, {
      int ss = sss/nbasis;
      int b  = sss%nbasis;
      calcComplex res = Zero();
      calcVector nbr;
      int ptype;
      StencilEntry *SE;

      for(int point=0;point<npoint;point++){

	SE=Stencil_v.GetEntry(ptype,point,ss);
	  
	if(SE->_is_local) { 
	  nbr = coalescedReadPermute(in_v[SE->_offset],ptype,SE->_permute);
	} else {
	  nbr = coalescedRead(Stencil_v.CommBuf()[SE->_offset]);
	}
	acceleratorSynchronise();

	for(int bb=0;bb<nbasis;bb++) {
	  res = res + coalescedRead(Aview_p[point][ss](b,bb))*nbr(bb);
	}
      }
      coalescedWrite(out_v[ss](b),res);
      });

    for(int p=0;p<geom.npoint;p++) AcceleratorViewContainer[p].ViewClose();
  };

  void Mdag (const CoarseVector &in, CoarseVector &out)
  {
    if(hermitian) {
      // corresponds to Petrov-Galerkin coarsening
      return M(in,out);
    } else {
      // corresponds to Galerkin coarsening
      return MdagNonHermitian(in, out);
    }
  };

  void MdagNonHermitian(const CoarseVector &in, CoarseVector &out)
  {
    conformable(_grid,in.Grid());
    conformable(in.Grid(),out.Grid());
    out.Checkerboard() = in.Checkerboard();

    SimpleCompressor<siteVector> compressor;

    Stencil.HaloExchange(in,compressor);
    autoView( in_v , in, AcceleratorRead);
    autoView( out_v , out, AcceleratorWrite);
    autoView( Stencil_v  , Stencil, AcceleratorRead);
    int npoint = geom.npoint;
    typedef LatticeView<Cobj> Aview;

    Vector<Aview> AcceleratorViewContainer;

    for(int p=0;p<geom.npoint;p++) AcceleratorViewContainer.push_back(A[p].View(AcceleratorRead));
    Aview *Aview_p = & AcceleratorViewContainer[0];

    const int Nsimd = CComplex::Nsimd();
    typedef decltype(coalescedRead(in_v[0])) calcVector;
    typedef decltype(coalescedRead(in_v[0](0))) calcComplex;

    int osites=Grid()->oSites();

    Vector<int> points(geom.npoint, 0);
    for(int p=0; p<geom.npoint; p++)
      points[p] = geom.points_dagger[p];

    auto points_p = &points[0];

    RealD* dag_factor_p = &dag_factor[0];

    accelerator_for(sss, Grid()->oSites()*nbasis, Nsimd, {
      int ss = sss/nbasis;
      int b  = sss%nbasis;
      calcComplex res = Zero();
      calcVector nbr;
      int ptype;
      StencilEntry *SE;

      for(int p=0;p<npoint;p++){
        int point = points_p[p];

	SE=Stencil_v.GetEntry(ptype,point,ss);

	if(SE->_is_local) {
	  nbr = coalescedReadPermute(in_v[SE->_offset],ptype,SE->_permute);
	} else {
	  nbr = coalescedRead(Stencil_v.CommBuf()[SE->_offset]);
	}
	acceleratorSynchronise();

	for(int bb=0;bb<nbasis;bb++) {
	  res = res + dag_factor_p[b*nbasis+bb]*coalescedRead(Aview_p[point][ss](b,bb))*nbr(bb);
	}
      }
      coalescedWrite(out_v[ss](b),res);
      });

    for(int p=0;p<geom.npoint;p++) AcceleratorViewContainer[p].ViewClose();
  }

  void MdirComms(const CoarseVector &in)
  {
    SimpleCompressor<siteVector> compressor;
    Stencil.HaloExchange(in,compressor);
  }
  void MdirCalc(const CoarseVector &in, CoarseVector &out, int point)
  {
    conformable(_grid,in.Grid());
    conformable(_grid,out.Grid());
    out.Checkerboard() = in.Checkerboard();

    typedef LatticeView<Cobj> Aview;
    Vector<Aview> AcceleratorViewContainer;
    for(int p=0;p<geom.npoint;p++) AcceleratorViewContainer.push_back(A[p].View(AcceleratorRead));
    Aview *Aview_p = & AcceleratorViewContainer[0];

    autoView( out_v , out, AcceleratorWrite);
    autoView( in_v  , in, AcceleratorRead);
    autoView( Stencil_v  , Stencil, AcceleratorRead);

    const int Nsimd = CComplex::Nsimd();
    typedef decltype(coalescedRead(in_v[0])) calcVector;
    typedef decltype(coalescedRead(in_v[0](0))) calcComplex;

    accelerator_for(sss, Grid()->oSites()*nbasis, Nsimd, {
      int ss = sss/nbasis;
      int b  = sss%nbasis;
      calcComplex res = Zero();
      calcVector nbr;
      int ptype;
      StencilEntry *SE;

      SE=Stencil_v.GetEntry(ptype,point,ss);
	  
      if(SE->_is_local) { 
	nbr = coalescedReadPermute(in_v[SE->_offset],ptype,SE->_permute);
      } else {
	nbr = coalescedRead(Stencil_v.CommBuf()[SE->_offset]);
      }
      acceleratorSynchronise();

      for(int bb=0;bb<nbasis;bb++) {
	res = res + coalescedRead(Aview_p[point][ss](b,bb))*nbr(bb);
      }
      coalescedWrite(out_v[ss](b),res);
    });
    for(int p=0;p<geom.npoint;p++) AcceleratorViewContainer[p].ViewClose();
  }
  void MdirAll(const CoarseVector &in,std::vector<CoarseVector> &out)
  {
    this->MdirComms(in);
    int ndir=geom.npoint-1;
    if ((out.size()!=ndir)&&(out.size()!=ndir+1)) { 
      std::cout <<"MdirAll out size "<< out.size()<<std::endl;
      std::cout <<"MdirAll ndir "<< ndir<<std::endl;
      assert(0);
    }
    for(int p=0;p<ndir;p++){
      MdirCalc(in,out[p],p);
    }
  };
  void Mdir(const CoarseVector &in, CoarseVector &out, int dir, int disp){

    this->MdirComms(in);

    MdirCalc(in,out,geom.point(dir,disp));
  };

  void Mdiag(const CoarseVector &in, CoarseVector &out)
  {
    int point=geom.npoint-1;
    MdirCalc(in, out, point); // No comms
  };

  void Mooee(const CoarseVector &in, CoarseVector &out) {
    MooeeInternal(in, out, DaggerNo, InverseNo);
  }

  void MooeeInv(const CoarseVector &in, CoarseVector &out) {
    MooeeInternal(in, out, DaggerNo, InverseYes);
  }

  void MooeeDag(const CoarseVector &in, CoarseVector &out) {
    MooeeInternal(in, out, DaggerYes, InverseNo);
  }

  void MooeeInvDag(const CoarseVector &in, CoarseVector &out) {
    MooeeInternal(in, out, DaggerYes, InverseYes);
  }

  void Meooe(const CoarseVector &in, CoarseVector &out) {
    if(in.Checkerboard() == Odd) {
      DhopEO(in, out, DaggerNo);
    } else {
      DhopOE(in, out, DaggerNo);
    }
  }

  void MeooeDag(const CoarseVector &in, CoarseVector &out) {
    if(in.Checkerboard() == Odd) {
      DhopEO(in, out, DaggerYes);
    } else {
      DhopOE(in, out, DaggerYes);
    }
  }

  void Dhop(const CoarseVector &in, CoarseVector &out, int dag) {
    conformable(in.Grid(), _grid); // verifies full grid
    conformable(in.Grid(), out.Grid());

    out.Checkerboard() = in.Checkerboard();

    DhopInternal(Stencil, A, in, out, dag);
  }

  void DhopOE(const CoarseVector &in, CoarseVector &out, int dag) {
    conformable(in.Grid(), _cbgrid);    // verifies half grid
    conformable(in.Grid(), out.Grid()); // drops the cb check

    assert(in.Checkerboard() == Even);
    out.Checkerboard() = Odd;

    DhopInternal(StencilEven, Aodd, in, out, dag);
  }

  void DhopEO(const CoarseVector &in, CoarseVector &out, int dag) {
    conformable(in.Grid(), _cbgrid);    // verifies half grid
    conformable(in.Grid(), out.Grid()); // drops the cb check

    assert(in.Checkerboard() == Odd);
    out.Checkerboard() = Even;

    DhopInternal(StencilOdd, Aeven, in, out, dag);
  }

  void MooeeInternal(const CoarseVector &in, CoarseVector &out, int dag, int inv) {
    out.Checkerboard() = in.Checkerboard();
    assert(in.Checkerboard() == Odd || in.Checkerboard() == Even);

    CoarseMatrix *Aself = nullptr;
    if(in.Grid()->_isCheckerBoarded) {
      if(in.Checkerboard() == Odd) {
        Aself = (inv) ? &AselfInvOdd : &Aodd[geom.npoint-1];
        DselfInternal(StencilOdd, *Aself, in, out, dag);
      } else {
        Aself = (inv) ? &AselfInvEven : &Aeven[geom.npoint-1];
        DselfInternal(StencilEven, *Aself, in, out, dag);
      }
    } else {
      Aself = (inv) ? &AselfInv : &A[geom.npoint-1];
      DselfInternal(Stencil, *Aself, in, out, dag);
    }
    assert(Aself != nullptr);
  }

  void DselfInternal(CartesianStencil<siteVector,siteVector,DefaultImplParams> &st, CoarseMatrix &a,
                       const CoarseVector &in, CoarseVector &out, int dag) {
    int point = geom.npoint-1;
    autoView( out_v, out, AcceleratorWrite);
    autoView( in_v,  in,  AcceleratorRead);
    autoView( st_v,  st,  AcceleratorRead);
    autoView( a_v,   a,   AcceleratorRead);

    const int Nsimd = CComplex::Nsimd();
    typedef decltype(coalescedRead(in_v[0])) calcVector;
    typedef decltype(coalescedRead(in_v[0](0))) calcComplex;

    RealD* dag_factor_p = &dag_factor[0];

    if(dag) {
      accelerator_for(sss, in.Grid()->oSites()*nbasis, Nsimd, {
        int ss = sss/nbasis;
        int b  = sss%nbasis;
        calcComplex res = Zero();
        calcVector nbr;
        int ptype;
        StencilEntry *SE;

        SE=st_v.GetEntry(ptype,point,ss);

        if(SE->_is_local) {
          nbr = coalescedReadPermute(in_v[SE->_offset],ptype,SE->_permute);
        } else {
          nbr = coalescedRead(st_v.CommBuf()[SE->_offset]);
        }
        acceleratorSynchronise();

        for(int bb=0;bb<nbasis;bb++) {
          res = res + dag_factor_p[b*nbasis+bb]*coalescedRead(a_v[ss](b,bb))*nbr(bb);
        }
        coalescedWrite(out_v[ss](b),res);
      });
    } else {
      accelerator_for(sss, in.Grid()->oSites()*nbasis, Nsimd, {
        int ss = sss/nbasis;
        int b  = sss%nbasis;
        calcComplex res = Zero();
        calcVector nbr;
        int ptype;
        StencilEntry *SE;

        SE=st_v.GetEntry(ptype,point,ss);

        if(SE->_is_local) {
          nbr = coalescedReadPermute(in_v[SE->_offset],ptype,SE->_permute);
        } else {
          nbr = coalescedRead(st_v.CommBuf()[SE->_offset]);
        }
        acceleratorSynchronise();

        for(int bb=0;bb<nbasis;bb++) {
          res = res + coalescedRead(a_v[ss](b,bb))*nbr(bb);
        }
        coalescedWrite(out_v[ss](b),res);
      });
    }
  }

  void DhopInternal(CartesianStencil<siteVector,siteVector,DefaultImplParams> &st, std::vector<CoarseMatrix> &a,
                    const CoarseVector &in, CoarseVector &out, int dag) {
    SimpleCompressor<siteVector> compressor;

    st.HaloExchange(in,compressor);
    autoView( in_v,  in,  AcceleratorRead);
    autoView( out_v, out, AcceleratorWrite);
    autoView( st_v , st,  AcceleratorRead);
    typedef LatticeView<Cobj> Aview;

    // determine in what order we need the points
    int npoint = geom.npoint-1;
    Vector<int> points(npoint, 0);
    for(int p=0; p<npoint; p++)
      points[p] = (dag && !hermitian) ? geom.points_dagger[p] : p;

    auto points_p = &points[0];

    Vector<Aview> AcceleratorViewContainer;
    for(int p=0;p<npoint;p++) AcceleratorViewContainer.push_back(a[p].View(AcceleratorRead));
    Aview *Aview_p = & AcceleratorViewContainer[0];

    const int Nsimd = CComplex::Nsimd();
    typedef decltype(coalescedRead(in_v[0])) calcVector;
    typedef decltype(coalescedRead(in_v[0](0))) calcComplex;

    RealD* dag_factor_p = &dag_factor[0];

    if(dag) {
      accelerator_for(sss, in.Grid()->oSites()*nbasis, Nsimd, {
        int ss = sss/nbasis;
        int b  = sss%nbasis;
        calcComplex res = Zero();
        calcVector nbr;
        int ptype;
        StencilEntry *SE;

        for(int p=0;p<npoint;p++){
          int point = points_p[p];
          SE=st_v.GetEntry(ptype,point,ss);

          if(SE->_is_local) {
            nbr = coalescedReadPermute(in_v[SE->_offset],ptype,SE->_permute);
          } else {
            nbr = coalescedRead(st_v.CommBuf()[SE->_offset]);
          }
          acceleratorSynchronise();

          for(int bb=0;bb<nbasis;bb++) {
            res = res + dag_factor_p[b*nbasis+bb]*coalescedRead(Aview_p[point][ss](b,bb))*nbr(bb);
          }
        }
        coalescedWrite(out_v[ss](b),res);
      });
    } else {
      accelerator_for(sss, in.Grid()->oSites()*nbasis, Nsimd, {
        int ss = sss/nbasis;
        int b  = sss%nbasis;
        calcComplex res = Zero();
        calcVector nbr;
        int ptype;
        StencilEntry *SE;

        for(int p=0;p<npoint;p++){
          int point = points_p[p];
          SE=st_v.GetEntry(ptype,point,ss);

          if(SE->_is_local) {
            nbr = coalescedReadPermute(in_v[SE->_offset],ptype,SE->_permute);
          } else {
            nbr = coalescedRead(st_v.CommBuf()[SE->_offset]);
          }
          acceleratorSynchronise();

          for(int bb=0;bb<nbasis;bb++) {
            res = res + coalescedRead(Aview_p[point][ss](b,bb))*nbr(bb);
          }
        }
        coalescedWrite(out_v[ss](b),res);
      });
    }

    for(int p=0;p<npoint;p++) AcceleratorViewContainer[p].ViewClose();
  }
  
  CoarsenedMatrix(GridCartesian &CoarseGrid, int hermitian_=0) 	:
    _grid(&CoarseGrid),
    _cbgrid(new GridRedBlackCartesian(&CoarseGrid)),
    geom(CoarseGrid._ndimension),
    hermitian(hermitian_),
    Stencil(&CoarseGrid,geom.npoint,Even,geom.directions,geom.displacements),
    StencilEven(_cbgrid,geom.npoint,Even,geom.directions,geom.displacements),
    StencilOdd(_cbgrid,geom.npoint,Odd,geom.directions,geom.displacements),
    A(geom.npoint,&CoarseGrid),
    Aeven(geom.npoint,_cbgrid),
    Aodd(geom.npoint,_cbgrid),
    AselfInv(&CoarseGrid),
    AselfInvEven(_cbgrid),
    AselfInvOdd(_cbgrid),
    dag_factor(nbasis*nbasis)
  {
    fillFactor();
  };

  CoarsenedMatrix(GridCartesian &CoarseGrid, GridRedBlackCartesian &CoarseRBGrid, int hermitian_=0) 	:

    _grid(&CoarseGrid),
    _cbgrid(&CoarseRBGrid),
    geom(CoarseGrid._ndimension),
    hermitian(hermitian_),
    Stencil(&CoarseGrid,geom.npoint,Even,geom.directions,geom.displacements),
    StencilEven(&CoarseRBGrid,geom.npoint,Even,geom.directions,geom.displacements),
    StencilOdd(&CoarseRBGrid,geom.npoint,Odd,geom.directions,geom.displacements),
    A(geom.npoint,&CoarseGrid),
    Aeven(geom.npoint,&CoarseRBGrid),
    Aodd(geom.npoint,&CoarseRBGrid),
    AselfInv(&CoarseGrid),
    AselfInvEven(&CoarseRBGrid),
    AselfInvOdd(&CoarseRBGrid),
    dag_factor(nbasis*nbasis)
  {
    fillFactor();
  };

  void fillFactor() {
    Eigen::MatrixXd dag_factor_eigen = Eigen::MatrixXd::Ones(nbasis, nbasis);
    if(!hermitian) {
      const int nb = nbasis/2;
      dag_factor_eigen.block(0,nb,nb,nb) *= -1.0;
      dag_factor_eigen.block(nb,0,nb,nb) *= -1.0;
    }

    // GPU readable prefactor
    thread_for(i, nbasis*nbasis, {
      int j = i/nbasis;
      int k = i%nbasis;
      dag_factor[i] = dag_factor_eigen(j, k);
    });
  }

  void CoarsenOperator(GridBase *FineGrid,LinearOperatorBase<Lattice<Fobj> > &linop,
		       Aggregation<Fobj,CComplex,nbasis> & Subspace)
  {
    typedef Lattice<typename Fobj::tensor_reduced> FineComplexField;
    typedef typename Fobj::scalar_type scalar_type;

    std::cout << GridLogMessage<< "CoarsenMatrix "<< std::endl;

    FineComplexField one(FineGrid); one=scalar_type(1.0,0.0);
    FineComplexField zero(FineGrid); zero=scalar_type(0.0,0.0);

    std::vector<FineComplexField> masks(geom.npoint,FineGrid);
    FineComplexField imask(FineGrid); // contributions from within this block
    FineComplexField omask(FineGrid); // contributions from outwith this block

    FineComplexField evenmask(FineGrid);
    FineComplexField oddmask(FineGrid); 

    FineField     phi(FineGrid);
    FineField     tmp(FineGrid);
    FineField     zz(FineGrid); zz=Zero();
    FineField    Mphi(FineGrid);
    FineField    Mphie(FineGrid);
    FineField    Mphio(FineGrid);
    std::vector<FineField>     Mphi_p(geom.npoint,FineGrid);

    Lattice<iScalar<vInteger> > coor (FineGrid);
    Lattice<iScalar<vInteger> > bcoor(FineGrid);
    Lattice<iScalar<vInteger> > bcb  (FineGrid); bcb = Zero();

    CoarseVector iProj(Grid()); 
    CoarseVector oProj(Grid()); 
    CoarseVector SelfProj(Grid()); 
    CoarseComplexField iZProj(Grid()); 
    CoarseComplexField oZProj(Grid()); 

    CoarseScalar InnerProd(Grid()); 

    std::cout << GridLogMessage<< "CoarsenMatrix Orthog "<< std::endl;
    // Orthogonalise the subblocks over the basis
    blockOrthogonalise(InnerProd,Subspace.subspace);

    // Compute the matrix elements of linop between this orthonormal
    // set of vectors.
    std::cout << GridLogMessage<< "CoarsenMatrix masks "<< std::endl;
    int self_stencil=-1;
    for(int p=0;p<geom.npoint;p++)
    { 
      int dir   = geom.directions[p];
      int disp  = geom.displacements[p];
      A[p]=Zero();
      if( geom.displacements[p]==0){
	self_stencil=p;
      }

      Integer block=(FineGrid->_rdimensions[dir])/(Grid()->_rdimensions[dir]);

      LatticeCoordinate(coor,dir);

      ///////////////////////////////////////////////////////
      // Work out even and odd block checkerboarding for fast diagonal term
      ///////////////////////////////////////////////////////
      if ( disp==1 ) {
	bcb   = bcb + div(coor,block);
      }
	
      if ( disp==0 ) {
	  masks[p]= Zero();
      } else if ( disp==1 ) {
	masks[p] = where(mod(coor,block)==(block-1),one,zero);
      } else if ( disp==-1 ) {
	masks[p] = where(mod(coor,block)==(Integer)0,one,zero);
      }
    }
    evenmask = where(mod(bcb,2)==(Integer)0,one,zero);
    oddmask  = one-evenmask;

    assert(self_stencil!=-1);

    for(int i=0;i<nbasis;i++){

      phi=Subspace.subspace[i];

      std::cout << GridLogMessage<< "CoarsenMatrix vector "<<i << std::endl;
      linop.OpDirAll(phi,Mphi_p);
      linop.OpDiag  (phi,Mphi_p[geom.npoint-1]);

      for(int p=0;p<geom.npoint;p++){ 

	Mphi = Mphi_p[p];

	int dir   = geom.directions[p];
	int disp  = geom.displacements[p];

	if ( (disp==-1) || (!hermitian ) ) {

	  ////////////////////////////////////////////////////////////////////////
	  // Pick out contributions coming from this cell and neighbour cell
	  ////////////////////////////////////////////////////////////////////////
	  omask = masks[p];
	  imask = one-omask;
	
	  for(int j=0;j<nbasis;j++){
	    
	    blockMaskedInnerProduct(oZProj,omask,Subspace.subspace[j],Mphi);
	    
	    autoView( iZProj_v , iZProj, AcceleratorRead) ;
	    autoView( oZProj_v , oZProj, AcceleratorRead) ;
	    autoView( A_p     ,  A[p], AcceleratorWrite);
	    autoView( A_self  , A[self_stencil], AcceleratorWrite);

	    accelerator_for(ss, Grid()->oSites(), Fobj::Nsimd(),{ coalescedWrite(A_p[ss](j,i),oZProj_v(ss)); });
	    if ( hermitian && (disp==-1) ) {
	      for(int pp=0;pp<geom.npoint;pp++){// Find the opposite link and set <j|A|i> = <i|A|j>*
		int dirp   = geom.directions[pp];
		int dispp  = geom.displacements[pp];
		if ( (dirp==dir) && (dispp==1) ){
		  auto sft = conjugate(Cshift(oZProj,dir,1));
		  autoView( sft_v    ,  sft  , AcceleratorWrite);
		  autoView( A_pp     ,  A[pp], AcceleratorWrite);
		  accelerator_for(ss, Grid()->oSites(), Fobj::Nsimd(),{ coalescedWrite(A_pp[ss](i,j),sft_v(ss)); });
		}
	      }
	    }

	  }
	}
      }

      ///////////////////////////////////////////
      // Faster alternate self coupling.. use hermiticity to save 2x
      ///////////////////////////////////////////
      {
	mult(tmp,phi,evenmask);  linop.Op(tmp,Mphie);
	mult(tmp,phi,oddmask );  linop.Op(tmp,Mphio);

	{
	  autoView( tmp_      , tmp, AcceleratorWrite);
	  autoView( evenmask_ , evenmask, AcceleratorRead);
	  autoView( oddmask_  ,  oddmask, AcceleratorRead);
	  autoView( Mphie_    ,  Mphie, AcceleratorRead);
	  autoView( Mphio_    ,  Mphio, AcceleratorRead);
	  accelerator_for(ss, FineGrid->oSites(), Fobj::Nsimd(),{ 
	      coalescedWrite(tmp_[ss],evenmask_(ss)*Mphie_(ss) + oddmask_(ss)*Mphio_(ss));
	    });
	}

	blockProject(SelfProj,tmp,Subspace.subspace);

	autoView( SelfProj_ , SelfProj, AcceleratorRead);
	autoView( A_self  , A[self_stencil], AcceleratorWrite);

	accelerator_for(ss, Grid()->oSites(), Fobj::Nsimd(),{
	  for(int j=0;j<nbasis;j++){
	    coalescedWrite(A_self[ss](j,i), SelfProj_(ss)(j));
	  }
	});

      }
    }
    if(hermitian) {
      std::cout << GridLogMessage << " ForceHermitian, new code "<<std::endl;
    }

    InvertSelfStencilLink(); std::cout << GridLogMessage << "Coarse self link inverted" << std::endl;
    FillHalfCbs(); std::cout << GridLogMessage << "Coarse half checkerboards filled" << std::endl;
  }

  void InvertSelfStencilLink() {
    std::cout << GridLogDebug << "CoarsenedMatrix::InvertSelfStencilLink" << std::endl;
    int localVolume = Grid()->lSites();

    typedef typename Cobj::scalar_object scalar_object;

    autoView(Aself_v,    A[geom.npoint-1], CpuRead);
    autoView(AselfInv_v, AselfInv,         CpuWrite);
    thread_for(site, localVolume, { // NOTE: Not able to bring this to GPU because of Eigen + peek/poke
      Eigen::MatrixXcd selfLinkEigen    = Eigen::MatrixXcd::Zero(nbasis, nbasis);
      Eigen::MatrixXcd selfLinkInvEigen = Eigen::MatrixXcd::Zero(nbasis, nbasis);

      scalar_object selfLink    = Zero();
      scalar_object selfLinkInv = Zero();

      Coordinate lcoor;

      Grid()->LocalIndexToLocalCoor(site, lcoor);
      peekLocalSite(selfLink, Aself_v, lcoor);

      for (int i = 0; i < nbasis; ++i)
        for (int j = 0; j < nbasis; ++j)
          selfLinkEigen(i, j) = static_cast<ComplexD>(TensorRemove(selfLink(i, j)));

      selfLinkInvEigen = selfLinkEigen.inverse();

      for(int i = 0; i < nbasis; ++i)
        for(int j = 0; j < nbasis; ++j)
          selfLinkInv(i, j) = selfLinkInvEigen(i, j);

      pokeLocalSite(selfLinkInv, AselfInv_v, lcoor);
    });
  }

  void FillHalfCbs() {
    std::cout << GridLogDebug << "CoarsenedMatrix::FillHalfCbs" << std::endl;
    for(int p = 0; p < geom.npoint; ++p) {
      pickCheckerboard(Even, Aeven[p], A[p]);
      pickCheckerboard(Odd, Aodd[p], A[p]);
    }
    pickCheckerboard(Even, AselfInvEven, AselfInv);
    pickCheckerboard(Odd, AselfInvOdd, AselfInv);
  }
};

NAMESPACE_END(Grid);
#endif
