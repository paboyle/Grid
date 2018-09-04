#pragma once
#include <Grid/Hadrons/Global.hpp>
#include <unsupported/Eigen/CXX11/Tensor>

BEGIN_HADRONS_NAMESPACE

template <typename FImpl>
class A2Autils 
{
public:
  typedef typename FImpl::ComplexField ComplexField;
  typedef typename FImpl::FermionField FermionField;
  typedef typename FImpl::PropagatorField PropagatorField;

  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iSpinMatrix<scalar_type> SpinMatrix_s;
  typedef iSinglet<vector_type> Scalar_v;
  typedef iSinglet<scalar_type> Scalar_s;

  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;

  static void MesonField(Eigen::Tensor<ComplexD,5> &mat, 
			 const FermionField *lhs_wi,
			 const FermionField *rhs_vj,
			 std::vector<Gamma::Algebra> gammas,
			 const std::vector<ComplexField > &mom,
			 int orthogdim);

  static void PionFieldWVmom(Eigen::Tensor<ComplexD,4> &mat, 
			     const FermionField *wi,
			     const FermionField *vj,
			     const std::vector<ComplexField > &mom,
			     int orthogdim);

  static void PionFieldXX(Eigen::Tensor<ComplexD,3> &mat, 
			  const FermionField *wi,
			  const FermionField *vj,
			  int orthogdim,
			  int g5);
  
  static void PionFieldWV(Eigen::Tensor<ComplexD,3> &mat, 
			  const FermionField *wi,
			  const FermionField *vj,
			  int orthogdim);
  static void PionFieldWW(Eigen::Tensor<ComplexD,3> &mat, 
			  const FermionField *wi,
			  const FermionField *wj,
			  int orthogdim);
  static void PionFieldVV(Eigen::Tensor<ComplexD,3> &mat, 
			  const FermionField *vi,
			  const FermionField *vj,
			  int orthogdim);

  static void ContractWWVV(std::vector<PropagatorField> &WWVV,
			   const Eigen::Tensor<ComplexD,3> &WW_sd,
			   const FermionField *vs,
			   const FermionField *vd);

  static void DeltaFeq2(int dt_min,int dt_max,
			Eigen::Tensor<ComplexD,2> &dF2_fig8,
			Eigen::Tensor<ComplexD,2> &dF2_trtr,
			Eigen::Tensor<ComplexD,1> &den0,
			Eigen::Tensor<ComplexD,1> &den1,
			Eigen::Tensor<ComplexD,3> &WW_sd, 
			const FermionField *vs,
			const FermionField *vd,
			int orthogdim);
};

template<class FImpl>
void A2Autils<FImpl>::MesonField(Eigen::Tensor<ComplexD,5> &mat, 
				 const FermionField *lhs_wi,
				 const FermionField *rhs_vj,
				 std::vector<Gamma::Algebra> gammas,
				 const std::vector<ComplexField > &mom,
				 int orthogdim) 
{
  typedef typename FImpl::SiteSpinor vobj;

  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iSpinMatrix<scalar_type> SpinMatrix_s;
  
  int Lblock = mat.dimension(3); 
  int Rblock = mat.dimension(4);

  GridBase *grid = lhs_wi[0]._grid;
  
  const int    Nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();

  int Nt     = grid->GlobalDimensions()[orthogdim];
  int Ngamma = gammas.size();
  int Nmom   = mom.size();

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  // will locally sum vectors first
  // sum across these down to scalars
  // splitting the SIMD
  int MFrvol = rd*Lblock*Rblock*Nmom;
  int MFlvol = ld*Lblock*Rblock*Nmom;

  Vector<SpinMatrix_v > lvSum(MFrvol);
  parallel_for (int r = 0; r < MFrvol; r++){
    lvSum[r] = zero;
  }

  Vector<SpinMatrix_s > lsSum(MFlvol);             
  parallel_for (int r = 0; r < MFlvol; r++){
    lsSum[r]=scalar_type(0.0);
  }

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];

  // potentially wasting cores here if local time extent too small
  parallel_for(int r=0;r<rd;r++){

    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){

	int ss= so+n*stride+b;

	for(int i=0;i<Lblock;i++){

	  auto left = conjugate(lhs_wi[i]._odata[ss]);

	  for(int j=0;j<Rblock;j++){

	    SpinMatrix_v vv;
	    auto right = rhs_vj[j]._odata[ss];
	    for(int s1=0;s1<Ns;s1++){
	    for(int s2=0;s2<Ns;s2++){
	      vv()(s1,s2)() = left()(s2)(0) * right()(s1)(0)
		+             left()(s2)(1) * right()(s1)(1)
		+             left()(s2)(2) * right()(s1)(2);
	    }}
	    
	    // After getting the sitewise product do the mom phase loop
	    int base = Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*r;
	    for ( int m=0;m<Nmom;m++){
	      int idx = m+base;
	      auto phase = mom[m]._odata[ss];
	      mac(&lvSum[idx],&vv,&phase);
	    }
	  
	  }
	}
      }
    }
  }


  // Sum across simd lanes in the plane, breaking out orthog dir.
  parallel_for(int rt=0;rt<rd;rt++){

    std::vector<int> icoor(Nd);
    std::vector<SpinMatrix_s> extracted(Nsimd);               

    for(int i=0;i<Lblock;i++){
    for(int j=0;j<Rblock;j++){
    for(int m=0;m<Nmom;m++){

      int ij_rdx = m+Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*rt;

      extract(lvSum[ij_rdx],extracted);

      for(int idx=0;idx<Nsimd;idx++){

	grid->iCoorFromIindex(icoor,idx);

	int ldx    = rt+icoor[orthogdim]*rd;

	int ij_ldx = m+Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*ldx;

	lsSum[ij_ldx]=lsSum[ij_ldx]+extracted[idx];

      }
    }}}
  }

  assert(mat.dimension(0) == Nmom);
  assert(mat.dimension(1) == Ngamma);
  assert(mat.dimension(2) == Nt);

  // ld loop and local only??
  int pd = grid->_processors[orthogdim];
  int pc = grid->_processor_coor[orthogdim];
  parallel_for_nest2(int lt=0;lt<ld;lt++)
  {
    for(int pt=0;pt<pd;pt++){
      int t = lt + pt*ld;
      if (pt == pc){
	for(int i=0;i<Lblock;i++){
	  for(int j=0;j<Rblock;j++){
	    for(int m=0;m<Nmom;m++){
	      int ij_dx = m+Nmom*i + Nmom*Lblock * j + Nmom*Lblock * Rblock * lt;
	      for(int mu=0;mu<Ngamma;mu++){
		// this is a bit slow
		mat(m,mu,t,i,j) = trace(lsSum[ij_dx]*Gamma(gammas[mu]));
	      }
	    }
	  }
	}
      } else { 
	const scalar_type zz(0.0);
	for(int i=0;i<Lblock;i++){
	  for(int j=0;j<Rblock;j++){
	    for(int mu=0;mu<Ngamma;mu++){
	      for(int m=0;m<Nmom;m++){
		mat(m,mu,t,i,j) =zz;
	      }
	    }
	  }
	}
      }
    }
  }

  ////////////////////////////////////////////////////////////////////
  // This global sum is taking as much as 50% of time on 16 nodes
  // Vector size is 7 x 16 x 32 x 16 x 16 x sizeof(complex) = 2MB - 60MB depending on volume
  // Healthy size that should suffice
  ////////////////////////////////////////////////////////////////////

  grid->GlobalSumVector(&mat(0,0,0,0,0),Nmom*Ngamma*Nt*Lblock*Rblock);

}


///////////////////////////////////////////////////////////////////
//Meson 
// Interested in
//
//      sum_x,y Trace[ G S(x,tx,y,ty) G S(y,ty,x,tx) ]
//
// Conventional meson field:
//                 
//    = sum_x,y Trace[ sum_j G |v_j(y,ty)> <w_j(x,tx)|  G sum_i |v_i(x,tx) ><w_i(y,ty)| ]
//    = sum_ij sum_x,y < w_j(x,tx)| G |v_i(x,tx) > <w_i(y,ty) (x)|G| v_j(y,ty) >
//    = sum_ij PI_ji(tx) PI_ij(ty)
//
// G5-Hermiticity
//
//      sum_x,y Trace[ G S(x,tx,y,ty) G S(y,ty,x,tx) ]
//    = sum_x,y Trace[ G S(x,tx,y,ty) G g5 S^dag(x,tx,y,ty) g5 ]
//    = sum_x,y Trace[ g5 G sum_j |v_j(y,ty)> <w_j(x,tx)|  G g5 sum_i   (|v_j(y,ty)> <w_i(x,tx)|)^dag ]      --  (*)
//
// NB:  Dag applies to internal indices spin,colour,complex
//
//    = sum_ij sum_x,y Trace[ g5 G |v_j(y,ty)> <w_j(x,tx)|  G g5  |w_i(x,tx)> <v_i(y,ty)| ]
//    = sum_ij sum_x,y <v_i(y,ty)|g5 G |v_j(y,ty)> <w_j(x,tx)|  G g5 |w_i(x,tx)> 
//    = sum_ij  PionVV(ty) PionWW(tx)
//
// (*) is only correct estimator if w_i and w_j come from distinct noise sets to preserve the kronecker
//     expectation value. Otherwise biased.
////////////////////////////////////////////////////////////////////

template<class FImpl>
void A2Autils<FImpl>::PionFieldXX(Eigen::Tensor<ComplexD,3> &mat, 
				  const FermionField *wi,
				  const FermionField *vj,
				  int orthogdim,
				  int g5) 
{
  int Lblock = mat.dimension(1); 
  int Rblock = mat.dimension(2);

  GridBase *grid = wi[0]._grid;
  
  const int    nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();

  int Nt     = grid->GlobalDimensions()[orthogdim];

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  // will locally sum vectors first
  // sum across these down to scalars
  // splitting the SIMD
  int MFrvol = rd*Lblock*Rblock;
  int MFlvol = ld*Lblock*Rblock;

  Vector<vector_type > lvSum(MFrvol);
  parallel_for (int r = 0; r < MFrvol; r++){
    lvSum[r] = zero;
  }

  Vector<scalar_type > lsSum(MFlvol);             
  parallel_for (int r = 0; r < MFlvol; r++){
    lsSum[r]=scalar_type(0.0);
  }

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];

  parallel_for(int r=0;r<rd;r++){

    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){

	int ss= so+n*stride+b;

	for(int i=0;i<Lblock;i++){

	  auto w = conjugate(wi[i]._odata[ss]);
	  if (g5) {
	    w()(2)(0) = - w()(2)(0);
	    w()(2)(1) = - w()(2)(1);
	    w()(2)(2) = - w()(2)(2);
	    w()(3)(0) = - w()(3)(0);
	    w()(3)(1) = - w()(3)(1);
	    w()(3)(2) = - w()(3)(2);
	  }
	  for(int j=0;j<Rblock;j++){

	    auto v = vj[j]._odata[ss];
	    auto vv = v()(0)(0);

	    vv =      w()(0)(0) * v()(0)(0)// Gamma5 Dirac basis explicitly written out
	      +       w()(0)(1) * v()(0)(1)
	      +       w()(0)(2) * v()(0)(2)
	      +       w()(1)(0) * v()(1)(0)
	      +       w()(1)(1) * v()(1)(1)
	      +       w()(1)(2) * v()(1)(2)
	      +       w()(2)(0) * v()(2)(0)
	      +       w()(2)(1) * v()(2)(1)
	      +       w()(2)(2) * v()(2)(2)
	      +       w()(3)(0) * v()(3)(0)
	      +       w()(3)(1) * v()(3)(1)
	      +       w()(3)(2) * v()(3)(2);
	    
	    int idx = i+Lblock*j+Lblock*Rblock*r;
	    lvSum[idx] = lvSum[idx]+vv;
	  }
	}
      }
    }
  }

  // Sum across simd lanes in the plane, breaking out orthog dir.
  parallel_for(int rt=0;rt<rd;rt++){

    std::vector<int> icoor(nd);
    iScalar<vector_type> temp; 
    std::vector<iScalar<scalar_type> > extracted(Nsimd);               

    for(int i=0;i<Lblock;i++){
    for(int j=0;j<Rblock;j++){

      int ij_rdx = i+Lblock*j+Lblock*Rblock*rt;

      temp._internal =lvSum[ij_rdx];
      extract(temp,extracted);

      for(int idx=0;idx<Nsimd;idx++){

	grid->iCoorFromIindex(icoor,idx);

	int ldx    = rt+icoor[orthogdim]*rd;

	int ij_ldx =i+Lblock*j+Lblock*Rblock*ldx;

	lsSum[ij_ldx]=lsSum[ij_ldx]+extracted[idx]._internal;

      }
    }}
  }

  assert(mat.dimension(0) == Nt);
  // ld loop and local only??
  int pd = grid->_processors[orthogdim];
  int pc = grid->_processor_coor[orthogdim];
  parallel_for_nest2(int lt=0;lt<ld;lt++)
  {
    for(int pt=0;pt<pd;pt++){
      int t = lt + pt*ld;
      if (pt == pc){
	for(int i=0;i<Lblock;i++){
	  for(int j=0;j<Rblock;j++){
	    int ij_dx = i + Lblock * j + Lblock * Rblock * lt;
	    mat(t,i,j) = lsSum[ij_dx];
	  }
	}
      } else { 
	const scalar_type zz(0.0);
	for(int i=0;i<Lblock;i++){
	  for(int j=0;j<Rblock;j++){
	    mat(t,i,j) =zz;
	  }
	}
      }
    }
  }

  grid->GlobalSumVector(&mat(0,0,0),Nt*Lblock*Rblock);
}

template<class FImpl>
void A2Autils<FImpl>::PionFieldWVmom(Eigen::Tensor<ComplexD,4> &mat, 
				     const FermionField *wi,
				     const FermionField *vj,
				     const std::vector<ComplexField > &mom,
				     int orthogdim) 
{
  int Lblock = mat.dimension(2); 
  int Rblock = mat.dimension(3);

  GridBase *grid = wi[0]._grid;
  
  const int    nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();

  int Nt     = grid->GlobalDimensions()[orthogdim];
  int Nmom   = mom.size();

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  // will locally sum vectors first
  // sum across these down to scalars
  // splitting the SIMD
  int MFrvol = rd*Lblock*Rblock*Nmom;
  int MFlvol = ld*Lblock*Rblock*Nmom;

  Vector<vector_type > lvSum(MFrvol);
  parallel_for (int r = 0; r < MFrvol; r++){
    lvSum[r] = zero;
  }

  Vector<scalar_type > lsSum(MFlvol);             
  parallel_for (int r = 0; r < MFlvol; r++){
    lsSum[r]=scalar_type(0.0);
  }

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];

  parallel_for(int r=0;r<rd;r++){

    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){

	int ss= so+n*stride+b;

	for(int i=0;i<Lblock;i++){

	  auto w = conjugate(wi[i]._odata[ss]);

	  for(int j=0;j<Rblock;j++){

	    auto v = vj[j]._odata[ss];

	    auto vv = w()(0)(0) * v()(0)(0)// Gamma5 Dirac basis explicitly written out
	      +       w()(0)(1) * v()(0)(1)
	      +       w()(0)(2) * v()(0)(2)
	      +       w()(1)(0) * v()(1)(0)
	      +       w()(1)(1) * v()(1)(1)
	      +       w()(1)(2) * v()(1)(2)
	      -       w()(2)(0) * v()(2)(0)
	      -       w()(2)(1) * v()(2)(1)
	      -       w()(2)(2) * v()(2)(2)
	      -       w()(3)(0) * v()(3)(0)
	      -       w()(3)(1) * v()(3)(1)
	      -       w()(3)(2) * v()(3)(2);

	    
	    // After getting the sitewise product do the mom phase loop
	    int base = Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*r;
	    for ( int m=0;m<Nmom;m++){
	      int idx = m+base;
	      auto phase = mom[m]._odata[ss];
	      mac(&lvSum[idx],&vv,&phase()()());
	    }
	  }
	}
      }
    }
  }


  // Sum across simd lanes in the plane, breaking out orthog dir.
  parallel_for(int rt=0;rt<rd;rt++){

    std::vector<int> icoor(nd);
    iScalar<vector_type> temp; 
    std::vector<iScalar<scalar_type> > extracted(Nsimd);               

      //    std::vector<scalar_type> extracted(Nsimd);               

    for(int i=0;i<Lblock;i++){
    for(int j=0;j<Rblock;j++){
    for(int m=0;m<Nmom;m++){

      int ij_rdx = m+Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*rt;

      temp._internal = lvSum[ij_rdx];
      extract(temp,extracted);

      for(int idx=0;idx<Nsimd;idx++){

	grid->iCoorFromIindex(icoor,idx);

	int ldx    = rt+icoor[orthogdim]*rd;

	int ij_ldx = m+Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*ldx;

	lsSum[ij_ldx]=lsSum[ij_ldx]+extracted[idx]._internal;

      }
    }}}
  }

  assert(mat.dimension(0) == Nmom);
  assert(mat.dimension(1) == Nt);
  // ld loop and local only??
  int pd = grid->_processors[orthogdim];
  int pc = grid->_processor_coor[orthogdim];
  parallel_for_nest2(int lt=0;lt<ld;lt++)
  {
    for(int pt=0;pt<pd;pt++){
      int t = lt + pt*ld;
      if (pt == pc){
	for(int i=0;i<Lblock;i++){
	  for(int j=0;j<Rblock;j++){
	    for(int m=0;m<Nmom;m++){
	      int ij_dx = m+Nmom*i + Nmom*Lblock * j + Nmom*Lblock * Rblock * lt;
	      mat(m,t,i,j) = lsSum[ij_dx];
	    }
	  }
	}
      } else { 
	const scalar_type zz(0.0);
	for(int i=0;i<Lblock;i++){
	  for(int j=0;j<Rblock;j++){
	    for(int m=0;m<Nmom;m++){
	      mat(m,t,i,j) =zz;
	    }
	  }
	}
      }
    }
  }

  grid->GlobalSumVector(&mat(0,0,0,0),Nmom*Nt*Lblock*Rblock);
}

template<class FImpl>
void A2Autils<FImpl>::PionFieldWV(Eigen::Tensor<ComplexD,3> &mat, 
				  const FermionField *wi,
				  const FermionField *vj,
				  int orthogdim) 
{
  const int g5=1;
  PionFieldXX(mat,wi,vj,orthogdim,g5);
}
template<class FImpl>
void A2Autils<FImpl>::PionFieldWW(Eigen::Tensor<ComplexD,3> &mat, 
				  const FermionField *wi,
				  const FermionField *wj,
				  int orthogdim) 
{
  const int nog5=0;
  PionFieldXX(mat,wi,wj,orthogdim,nog5);
}
template<class FImpl>
void A2Autils<FImpl>::PionFieldVV(Eigen::Tensor<ComplexD,3> &mat, 
				  const FermionField *vi,
				  const FermionField *vj,
				  int orthogdim) 
{
  const int nog5=0;
  PionFieldXX(mat,vi,vj,orthogdim,nog5);
}


////////////////////////////////////////////
// Schematic thoughts about more generalised four quark insertion
//
// --  Pupil or fig8 type topology (depending on flavour structure) done below
// --  Have  Bag style   --  WWVV  VVWW
//            _  _      
//           / \/ \   
//           \_/\_/   
//                                                                                                                                                                       
//  -   Kpipi style (two pion insertions) 
//                                K     
// *********
// Type 1                    
// *********
//                x
//            ___  _____ pi(b)  
//        K  /   \/____/          
//           \    \
//            \____\pi(a)   
//
//                        
//        W^W_sd V_s(x)  V^_d(xa) Wa(xa) Va^(x)    WW_bb' Vb Vb'(x)
//
//        (Kww.PIvw) VV ;   pi_ww.VV
// 
//  But Norman and Chris say g5 hermiticity not used, except on K
//
//        Kww    PIvw     PIvw       
//
//        (W^W_sd PIa_dd')   PIb_uu'  (v_s v_d' w_u v_u')(x)
// 
//  - Can form (Nmode^3)
//
//   (Kww . PIvw)_sd and then V_sV_d  tensor product contracting this
// 
//  - Can form 
//
//   (PIvw)_uu' W_uV_u'  tensor product.
//
// Both are lattice propagators.
//
// Contract with the same four quark type contraction as BK.
//
// *********
// Type 2
// *********
//            _  _____ pi 
//        K  / \/     |   
//           \_/\_____|pi 
//
// Norman/Chris would say
//
//  (Kww VV)(x) . (PIwv Piwv) VW(x)
//
// Full four quark via g5 hermiticity
//
//  (Kww VV)(x) . (PIww Pivw) VV(x)
//
// *********
// Type 3
// *********
//            ___  _____ pi 
//        K  /   \/     |   
//           \   /\     |
//            \  \/     |
//             \________|pi 
// 
// 
//   (V(x) . Kww . pi_vw . pivw . V(x)). WV(x)
//
// No difference possible to Norman, Chris, Diaqian
//
// *********
// Type 4
// *********
//            ___        pi 
//        K  /   \/\     |\  
//           \___/\/     |/
//                        pi
//            
//  (Kww VV ) WV        x   Tr(PIwv PIwv)
//        
// Could use alternate PI_ww PIvv for disco loop interesting comparison
//
// *********
// Primitives / Utilities for assembly
// *********
// 
// i)   Basic type for meson field - mode indexed object, whether WW, VW, VV etc..
// ii)  Multiply two meson fields (modes^3) ; use BLAS MKL via Eigen
// iii) Multiply and trace two meson fields ; use BLAS MKL via Eigen. The element wise product trick
// iv)  Contract a meson field (whether WW,VW, VV WV) with W and V fields to form LatticePropagator
// v)   L^3 sum a four quark contaction with two LatticePropagators.
//
// Use lambda functions to enable flexibility in v)
// Use lambda functions to unify the pion field / meson field contraction codes.
////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////
// DeltaF=2 contraction ; use exernal WW field for Kaon, anti Kaon sink
////////////////////////////////////////////////////////////////////////
//
// WW -- i vectors have adjoint, and j vectors not. 
//    -- Think of "i" as the strange quark, forward prop from 0
//    -- Think of "j" as the anti-down quark.
//
// WW_sd are w^dag_s w_d
//
// Hence VV vectors correspondingly are  v^dag_d,  v_s    from t=0
//                                  and  v^dag_d,  v_s    from t=dT
//
// There is an implicit g5 associated with each v^dag_d from use of g5 Hermiticity.
// The other gamma_5 lies in the WW external meson operator.
//
// From UKhadron wallbag.cc:
//
//   LatticePropagator anti_d0 =  adj( Gamma(G5) * Qd0 * Gamma(G5));
//   LatticePropagator anti_d1 =  adj( Gamma(G5) * Qd1 * Gamma(G5));
//
//   PR1 = Qs0 * Gamma(G5) * anti_d0;
//   PR2 = Qs1 * Gamma(G5) * anti_d1;
//
//   TR1 = trace( PR1 * G1 );
//   TR2 = trace( PR2 * G2 );
//   Wick1 = TR1 * TR2;
//
//   Wick2 = trace( PR1* G1 * PR2 * G2 );
//   // was      Wick2 = trace( Qs0 * Gamma(G5) * anti_d0 * G1 * Qs1 * Gamma(G5) * anti_d1 * G2 );
//
// TR TR(tx) = Wick1 = sum_x WW[t0]_sd < v^_d |g5 G| v_s>   WW[t1]_s'd' < v^_d' |g5 G| v_s'> |_{x,tx)
//           = sum_x [ Trace(WW[t0] VgV(t,x) )  x Trace( WW_[t1] VgV(t,x) ) ]
//           
//
// Calc all Nt Trace(WW VV) products at once, take Nt^2 products of these.
//
// Fig8(tx)  = Wick2 = sum_x WW[t0]_sd WW[t1]_s'd'  < v^_d |g5 G| v_s'> < v^_d' |g5 G| v_s> |_{x,tx}
//
//                   = sum_x Trace( WW[t0] VV[t,x] WW[t1] VV[t,x] )
// 
///////////////////////////////////////////////////////////////////////////////
//
// WW is w_s^dag (x) w_d       (G5 implicitly absorbed)
//
// WWVV will have spin-col (x) spin-col tensor.
//
// Want this to look like a strange fwd prop, anti-d prop.
// 
// Take WW_sd v^dag_d (x) v_s
// 

template<class FImpl>
void A2Autils<FImpl>::ContractWWVV(std::vector<PropagatorField> &WWVV,
				   const Eigen::Tensor<ComplexD,3> &WW_sd,
				   const FermionField *vs,
				   const FermionField *vd)
{
  GridBase *grid = vs[0]._grid;

  int nd    = grid->_ndimension;
  int Nsimd = grid->Nsimd();
  int N_t   = WW_sd.dimension(0);
  int N_s = WW_sd.dimension(1); 
  int N_d = WW_sd.dimension(2);

  int d_unroll = 32;// Empirical optimisation

  for(int t=0;t<N_t;t++){
    WWVV[t] = zero;
  }

  parallel_for(int ss=0;ss<grid->oSites();ss++){
    for(int d_o=0;d_o<N_d;d_o+=d_unroll){
      for(int t=0;t<N_t;t++){
      for(int s=0;s<N_s;s++){
	auto tmp1 = vs[s]._odata[ss];
	vobj tmp2 = zero;

	for(int d=d_o;d<MIN(d_o+d_unroll,N_d);d++){
	  Scalar_v coeff = WW_sd(t,s,d);
	  mac(&tmp2 ,& coeff, & vd[d]._odata[ss]);
	}

	//////////////////////////
	// Fast outer product of tmp1 with a sum of terms suppressed by d_unroll
	//////////////////////////
	tmp2 = conjugate(tmp2);
	for(int s1=0;s1<Ns;s1++){
	for(int s2=0;s2<Ns;s2++){
	  WWVV[t]._odata[ss]()(s1,s2)(0,0) += tmp1()(s1)(0)*tmp2()(s2)(0);
	  WWVV[t]._odata[ss]()(s1,s2)(0,1) += tmp1()(s1)(0)*tmp2()(s2)(1);
	  WWVV[t]._odata[ss]()(s1,s2)(0,2) += tmp1()(s1)(0)*tmp2()(s2)(2);
	  WWVV[t]._odata[ss]()(s1,s2)(1,0) += tmp1()(s1)(1)*tmp2()(s2)(0);
	  WWVV[t]._odata[ss]()(s1,s2)(1,1) += tmp1()(s1)(1)*tmp2()(s2)(1);
	  WWVV[t]._odata[ss]()(s1,s2)(1,2) += tmp1()(s1)(1)*tmp2()(s2)(2);
	  WWVV[t]._odata[ss]()(s1,s2)(2,0) += tmp1()(s1)(2)*tmp2()(s2)(0);
	  WWVV[t]._odata[ss]()(s1,s2)(2,1) += tmp1()(s1)(2)*tmp2()(s2)(1);
	  WWVV[t]._odata[ss]()(s1,s2)(2,2) += tmp1()(s1)(2)*tmp2()(s2)(2);
	}}

      }}
    }
  }
}


template<class FImpl>
void A2Autils<FImpl>::DeltaFeq2(int dt_min,int dt_max,
		      Eigen::Tensor<ComplexD,2> &dF2_fig8,
		      Eigen::Tensor<ComplexD,2> &dF2_trtr,
		      Eigen::Tensor<ComplexD,1> &den0,
		      Eigen::Tensor<ComplexD,1> &den1,
		      Eigen::Tensor<ComplexD,3> &WW_sd, 
		      const FermionField *vs,
		      const FermionField *vd,
		      int orthogdim)
{
  GridBase *grid = vs[0]._grid;
  LOG(Message) << "Computing A2A DeltaF=2 graph" << std::endl;

  auto G5 = Gamma(Gamma::Algebra::Gamma5);
  
  double nodes = grid->NodeCount();
  int nd    = grid->_ndimension;
  int Nsimd = grid->Nsimd();
  int N_t = WW_sd.dimension(0); 
  int N_s = WW_sd.dimension(1); 
  int N_d = WW_sd.dimension(2);

  assert(grid->GlobalDimensions()[orthogdim] == N_t);
  double vol         = 1.0;
  for(int dim=0;dim<nd;dim++){
    vol = vol * grid->GlobalDimensions()[dim];
  }

  double t_tot  = -usecond();
  std::vector<PropagatorField> WWVV (N_t,grid);

  double t_outer= -usecond();
  ContractWWVV(WWVV,WW_sd,&vs[0],&vd[0]);
  t_outer+=usecond();

  //////////////////////////////
  // Implicit gamma-5 
  //////////////////////////////
  for(int t=0;t<N_t;t++){
    WWVV[t] = WWVV[t]* G5 ; 
  }

  ////////////////////////////////////////////////////////
  // Contraction
  ////////////////////////////////////////////////////////
  dF2_trtr.resize(N_t,5);
  dF2_fig8.resize(N_t,5);
  den0.resize(N_t);
  den1.resize(N_t);
  for(int t=0;t<N_t;t++){
    for(int g=0;g<dF2_trtr.dimension(1);g++) dF2_trtr(t,g)= ComplexD(0.0);
    for(int g=0;g<dF2_fig8.dimension(1);g++) dF2_fig8(t,g)= ComplexD(0.0);
    den0(t) =ComplexD(0.0);
    den1(t) =ComplexD(0.0);
  }

  ComplexField D0(grid);   D0 = zero; // <P|A0> correlator from each wall
  ComplexField D1(grid);   D1 = zero;

  ComplexField O1_trtr(grid);  O1_trtr = zero;
  ComplexField O2_trtr(grid);  O2_trtr = zero;
  ComplexField O3_trtr(grid);  O3_trtr = zero;
  ComplexField O4_trtr(grid);  O4_trtr = zero;
  ComplexField O5_trtr(grid);  O5_trtr = zero;

  ComplexField O1_fig8(grid);  O1_fig8 = zero;
  ComplexField O2_fig8(grid);  O2_fig8 = zero;
  ComplexField O3_fig8(grid);  O3_fig8 = zero;
  ComplexField O4_fig8(grid);  O4_fig8 = zero;
  ComplexField O5_fig8(grid);  O5_fig8 = zero;

  //////////////////////////////////////////////////
  // Used to store appropriate correlation funcs
  //////////////////////////////////////////////////
  std::vector<TComplex>  C1;
  std::vector<TComplex>  C2;
  std::vector<TComplex>  C3;
  std::vector<TComplex>  C4;
  std::vector<TComplex>  C5;
  
  //////////////////////////////////////////////////////////
  // Could do AA, VV, SS, PP, TT and form linear combinations later.
  // Almost 2x. but for many modes, the above loop dominates.
  //////////////////////////////////////////////////////////
  double t_contr= -usecond();

  // Tr Tr Wick contraction
  auto VX = Gamma(Gamma::Algebra::GammaX);
  auto VY = Gamma(Gamma::Algebra::GammaY);
  auto VZ = Gamma(Gamma::Algebra::GammaZ);
  auto VT = Gamma(Gamma::Algebra::GammaT);
  
  auto AX = Gamma(Gamma::Algebra::GammaXGamma5);
  auto AY = Gamma(Gamma::Algebra::GammaYGamma5);
  auto AZ = Gamma(Gamma::Algebra::GammaZGamma5);
  auto AT = Gamma(Gamma::Algebra::GammaTGamma5);
  
  auto S  = Gamma(Gamma::Algebra::Identity);
  auto P  = Gamma(Gamma::Algebra::Gamma5);
  
  auto T0 = Gamma(Gamma::Algebra::SigmaXY);
  auto T1 = Gamma(Gamma::Algebra::SigmaXZ);
  auto T2 = Gamma(Gamma::Algebra::SigmaXT);
  auto T3 = Gamma(Gamma::Algebra::SigmaYZ);
  auto T4 = Gamma(Gamma::Algebra::SigmaYT);
  auto T5 = Gamma(Gamma::Algebra::SigmaZT);

  std::cout <<GridLogMessage << " dt " <<dt_min<<"..." <<dt_max<<std::endl;
  for(int t0=0;t0<N_t;t0++){
    std::cout <<GridLogMessage << " t0 " <<t0<<std::endl;
    //    for(int dt=dt_min;dt<dt_max;dt++){
    {     
      int dt  = dt_min;

      int t1 = (t0+dt)%N_t;
      std::cout <<GridLogMessage << " t1 " <<t1<<std::endl;

      parallel_for(int ss=0;ss<grid->oSites();ss++){

	auto VV0= WWVV[t0]._odata[ss];
	auto VV1= WWVV[t1]._odata[ss];

	O1_trtr._odata[ss] = trace(VX*VV0) * trace(VX*VV1)
	  +                  trace(VY*VV0) * trace(VY*VV1)
	  +                  trace(VZ*VV0) * trace(VZ*VV1)
	  +                  trace(VT*VV0) * trace(VT*VV1)
          +                  trace(AX*VV0) * trace(AX*VV1)
	  +                  trace(AY*VV0) * trace(AY*VV1)
	  +                  trace(AZ*VV0) * trace(AZ*VV1)
	  +                  trace(AT*VV0) * trace(AT*VV1);

	O2_trtr._odata[ss] = trace(VX*VV0) * trace(VX*VV1)
	  +                  trace(VY*VV0) * trace(VY*VV1)
	  +                  trace(VZ*VV0) * trace(VZ*VV1)
	  +                  trace(VT*VV0) * trace(VT*VV1)
          -                  trace(AX*VV0) * trace(AX*VV1)
	  -                  trace(AY*VV0) * trace(AY*VV1)
	  -                  trace(AZ*VV0) * trace(AZ*VV1)
	  -                  trace(AT*VV0) * trace(AT*VV1);

	O3_trtr._odata[ss] = trace(S*VV0) * trace(S*VV1)
	  +                  trace(P*VV0) * trace(P*VV1);

	O4_trtr._odata[ss] = trace(S*VV0) * trace(S*VV1)
	  -                  trace(P*VV0) * trace(P*VV1);

	O5_trtr._odata[ss] = trace(T0*VV0) * trace(T0*VV1)
	  +                  trace(T1*VV0) * trace(T1*VV1)
	  +                  trace(T2*VV0) * trace(T2*VV1)
	  +                  trace(T3*VV0) * trace(T3*VV1)
	  +                  trace(T4*VV0) * trace(T4*VV1)
	  +                  trace(T5*VV0) * trace(T5*VV1);

	////////////////////////////////////
	// Fig8 Wick contraction
	////////////////////////////////////

	O1_fig8._odata[ss]  = trace (VV0 * VX * VV1 * VX)
	  +                   trace (VV0 * VY * VV1 * VY)
	  +                   trace (VV0 * VZ * VV1 * VZ)
	  +                   trace (VV0 * VT * VV1 * VT)
          +                   trace (VV0 * AX * VV1 * AX)
	  +                   trace (VV0 * AY * VV1 * AY)
	  +                   trace (VV0 * AZ * VV1 * AZ)
	  +                   trace (VV0 * AT * VV1 * AT);

	O2_fig8._odata[ss]  = trace (VV0 * VX * VV1 * VX)
	  +                   trace (VV0 * VY * VV1 * VY)
	  +                   trace (VV0 * VZ * VV1 * VZ)
	  +                   trace (VV0 * VT * VV1 * VT)
          -                   trace (VV0 * AX * VV1 * AX)
	  -                   trace (VV0 * AY * VV1 * AY)
	  -                   trace (VV0 * AZ * VV1 * AZ)
	  -                   trace (VV0 * AT * VV1 * AT);

	O3_fig8._odata[ss]  = trace (VV0 * S * VV1 * S)
	  +                   trace (VV0 * P * VV1 * P);

	O4_fig8._odata[ss]  = trace (VV0 * S * VV1 * S)
	  -                   trace (VV0 * P * VV1 * P);

	O5_fig8._odata[ss] =  trace (VV0 * T0 * VV1 * T0)
	  +                   trace (VV0 * T1 * VV1 * T1)
	  +                   trace (VV0 * T2 * VV1 * T2)
	  +                   trace (VV0 * T3 * VV1 * T3)
	  +                   trace (VV0 * T4 * VV1 * T4)
	  +                   trace (VV0 * T5 * VV1 * T5);

	D0._odata[ss] = trace(AT*VV0);
	D1._odata[ss] = trace(AT*VV1);

      }

      sliceSum(O1_trtr,C1, orthogdim);
      sliceSum(O2_trtr,C2, orthogdim);
      sliceSum(O3_trtr,C3, orthogdim);
      sliceSum(O4_trtr,C4, orthogdim);
      sliceSum(O5_trtr,C5, orthogdim);


      for(int t=0;t<N_t;t++){// 2x from Wick contraction reordering
	dF2_trtr(t,0)+= 2.0*C1[(t+t0)%N_t]()()()/vol;
	dF2_trtr(t,1)+= 2.0*C2[(t+t0)%N_t]()()()/vol;
	dF2_trtr(t,2)+= 2.0*C3[(t+t0)%N_t]()()()/vol;
	dF2_trtr(t,3)+= 2.0*C4[(t+t0)%N_t]()()()/vol;
	dF2_trtr(t,4)+= 2.0*C5[(t+t0)%N_t]()()()/vol;
      }

      sliceSum(O1_fig8,C1, orthogdim);
      sliceSum(O2_fig8,C2, orthogdim);
      sliceSum(O3_fig8,C3, orthogdim);
      sliceSum(O4_fig8,C4, orthogdim);
      sliceSum(O5_fig8,C5, orthogdim);

      for(int t=0;t<N_t;t++){
	dF2_fig8(t,0)= 2.0*C1[(t+t0)%N_t]()()()/vol;
	dF2_fig8(t,1)= 2.0*C2[(t+t0)%N_t]()()()/vol;
	dF2_fig8(t,2)= 2.0*C3[(t+t0)%N_t]()()()/vol;
	dF2_fig8(t,3)= 2.0*C4[(t+t0)%N_t]()()()/vol;
	dF2_fig8(t,4)= 2.0*C5[(t+t0)%N_t]()()()/vol;
      }

      sliceSum(D0,C1, orthogdim);
      sliceSum(D1,C2, orthogdim);

      for(int t=0;t<N_t;t++){
	den0(t)+=C1[(t+t0)%N_t]()()()/vol;
	den1(t)+=C2[(t+t0)%N_t]()()()/vol;
      }

    }
  }
  t_contr +=usecond();
  
  t_tot+=usecond();
  double million=1.0e6;
  LOG(Message) << "Computing A2A DeltaF=2 graph t_tot      " << t_tot      /million << " s "<< std::endl;
  LOG(Message) << "Computing A2A DeltaF=2 graph t_outer    " << t_outer    /million << " s "<< std::endl;
  LOG(Message) << "Computing A2A DeltaF=2 graph t_contr    " << t_contr    /million << " s "<< std::endl;
}

END_HADRONS_NAMESPACE
