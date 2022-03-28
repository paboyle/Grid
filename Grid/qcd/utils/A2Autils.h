#pragma once
//#include <Grid/Hadrons/Global.hpp>
#include <Grid/Grid_Eigen_Tensor.h>

NAMESPACE_BEGIN(Grid);

#undef DELTA_F_EQ_2

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

  template <typename TensorType> // output: rank 5 tensor, e.g. Eigen::Tensor<ComplexD, 5>
  static void MesonField(TensorType &mat, 
			 const FermionField *lhs_wi,
			 const FermionField *rhs_vj,
			 std::vector<Gamma::Algebra> gammas,
			 const std::vector<ComplexField > &mom,
			 int orthogdim, double *t_kernel = nullptr, double *t_gsum = nullptr);

  template <typename TensorType> // output: rank 5 tensor, e.g. Eigen::Tensor<ComplexD, 5>
  static void StagMesonField(TensorType &mat,
                             const FermionField *lhs_wi,
                             const FermionField *rhs_vj,
                             std::vector<Gamma::Algebra> gammas,
                             const std::vector<ComplexField > &mom,
                             int orthogdim, double *t_kernel = nullptr, double *t_gsum = nullptr);
  template <typename TensorType> // output: rank 5 tensor, e.g. Eigen::Tensor<ComplexD, 5>
  static void StagMesonFieldMILC(TensorType &mat,
                             const FermionField *lhs_wi,
                             const FermionField *rhs_vj,
                             std::vector<Gamma::Algebra> gammas,
                             const std::vector<ComplexField > &mom,
                             int orthogdim, double *t_kernel = nullptr, double *t_gsum = nullptr);
  template <typename TensorType> // output: rank 5 tensor, e.g. Eigen::Tensor<ComplexD, 5>
  static void StagMesonFieldCC(TensorType &mat,
                               //const LatticeGaugeField &U,
                               const FermionField *lhs_wi,
                               const FermionField *rhs_vj,
                               std::vector<Gamma::Algebra> gammas,
                               const std::vector<ComplexField > &mom,
                               int orthogdim, double *t_kernel = nullptr, double *t_gsum = nullptr);

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
    
  static void StagGSPionField(Eigen::Tensor<ComplexD,3> &mat,
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

  template <typename TensorType> // output: rank 5 tensor, e.g. Eigen::Tensor<ComplexD, 5>
  static void AslashField(TensorType &mat, 
        const FermionField *lhs_wi,
        const FermionField *rhs_vj,
        const std::vector<ComplexField> &emB0,
        const std::vector<ComplexField> &emB1,
        int orthogdim, double *t_kernel = nullptr, double *t_gsum = nullptr);

  template <typename TensorType>
  typename std::enable_if<(std::is_same<Eigen::Tensor<ComplexD,3>, TensorType>::value ||
                           std::is_same<Eigen::TensorMap<Eigen::Tensor<Complex, 3, Eigen::RowMajor>>, TensorType>::value),
                           void>::type
  static ContractWWVV(std::vector<PropagatorField> &WWVV,
			   const TensorType &WW_sd,
			   const FermionField *vs,
			   const FermionField *vd);

  template <typename TensorType>
  typename std::enable_if<!(std::is_same<Eigen::Tensor<ComplexD,3>, TensorType>::value ||
                            std::is_same<Eigen::TensorMap<Eigen::Tensor<Complex, 3, Eigen::RowMajor>>, TensorType>::value),
                            void>::type
  static ContractWWVV(std::vector<PropagatorField> &WWVV,
			   const TensorType &WW_sd,
			   const FermionField *vs,
			   const FermionField *vd);

  static void ContractFourQuarkColourDiagonal(const PropagatorField &WWVV0,
					      const PropagatorField &WWVV1,
					      const std::vector<Gamma> &gamma0,
					      const std::vector<Gamma> &gamma1,
					      ComplexField &O_trtr,
					      ComplexField &O_fig8);

  static void ContractFourQuarkColourMix(const PropagatorField &WWVV0,
					 const PropagatorField &WWVV1,
					 const std::vector<Gamma> &gamma0,
					 const std::vector<Gamma> &gamma1,
					 ComplexField &O_trtr,
					 ComplexField &O_fig8);
#ifdef DELTA_F_EQ_2
  static void DeltaFeq2(int dt_min,int dt_max,
			Eigen::Tensor<ComplexD,2> &dF2_fig8,
			Eigen::Tensor<ComplexD,2> &dF2_trtr,
			Eigen::Tensor<ComplexD,2> &dF2_fig8_mix,
			Eigen::Tensor<ComplexD,2> &dF2_trtr_mix,
			Eigen::Tensor<ComplexD,1> &denom_A0,
			Eigen::Tensor<ComplexD,1> &denom_P,
			Eigen::Tensor<ComplexD,3> &WW_sd, 
			const FermionField *vs,
			const FermionField *vd,
			int orthogdim);
#endif
private:
  inline static void OuterProductWWVV(PropagatorField &WWVV,
                               const vobj &lhs,
                               const vobj &rhs,
                               const int Ns, const int ss);
};

template <class FImpl>
template <typename TensorType>
void A2Autils<FImpl>::MesonField(TensorType &mat, 
				 const FermionField *lhs_wi,
				 const FermionField *rhs_vj,
				 std::vector<Gamma::Algebra> gammas,
				 const std::vector<ComplexField > &mom,
				 int orthogdim, double *t_kernel, double *t_gsum) 
{
  typedef typename FImpl::SiteSpinor vobj;

  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iSpinMatrix<scalar_type> SpinMatrix_s;
  
  int Lblock = mat.dimension(3); 
  int Rblock = mat.dimension(4);

  GridBase *grid = lhs_wi[0].Grid();
  
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
  thread_for( r, MFrvol,{
    lvSum[r] = Zero();
  });

  Vector<SpinMatrix_s > lsSum(MFlvol);             
  thread_for(r,MFlvol,{
    lsSum[r]=scalar_type(0.0);
  });

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];

  // potentially wasting cores here if local time extent too small
  if (t_kernel) *t_kernel = -usecond();
  thread_for(r,rd,{

    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){

	int ss= so+n*stride+b;

	for(int i=0;i<Lblock;i++){

	  // Recreate view potentially expensive outside fo UVM mode
	  autoView(lhs_v,lhs_wi[i],CpuRead);
	  auto left = conjugate(lhs_v[ss]);
	  for(int j=0;j<Rblock;j++){

	    SpinMatrix_v vv;
	    // Recreate view potentially expensive outside fo UVM mode
	    autoView(rhs_v,rhs_vj[j],CpuRead);
	    auto right = rhs_v[ss];
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
	      autoView(mom_v,mom[m],CpuRead);
	      auto phase = mom_v[ss];
	      mac(&lvSum[idx],&vv,&phase);
	    }
	  }
	}
      }
    }
  });

  // Sum across simd lanes in the plane, breaking out orthog dir.
  thread_for(rt,rd,{

    Coordinate icoor(Nd);
    ExtractBuffer<SpinMatrix_s> extracted(Nsimd);               

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
  });
  if (t_kernel) *t_kernel += usecond();
  assert(mat.dimension(0) == Nmom);
  assert(mat.dimension(1) == Ngamma);
  assert(mat.dimension(2) == Nt);

  // ld loop and local only??
  int pd = grid->_processors[orthogdim];
  int pc = grid->_processor_coor[orthogdim];
  thread_for_collapse(2,lt,ld,{
    for(int pt=0;pt<pd;pt++){
      int t = lt + pt*ld;
      if (pt == pc){
	for(int i=0;i<Lblock;i++){
	  for(int j=0;j<Rblock;j++){
	    for(int m=0;m<Nmom;m++){
	      int ij_dx = m+Nmom*i + Nmom*Lblock * j + Nmom*Lblock * Rblock * lt;
	      for(int mu=0;mu<Ngamma;mu++){
		// this is a bit slow
		mat(m,mu,t,i,j) = trace(lsSum[ij_dx]*Gamma(gammas[mu]))()()();
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
  });

  ////////////////////////////////////////////////////////////////////
  // This global sum is taking as much as 50% of time on 16 nodes
  // Vector size is 7 x 16 x 32 x 16 x 16 x sizeof(complex) = 2MB - 60MB depending on volume
  // Healthy size that should suffice
  ////////////////////////////////////////////////////////////////////
  if (t_gsum) *t_gsum = -usecond();
  grid->GlobalSumVector(&mat(0,0,0,0,0),Nmom*Ngamma*Nt*Lblock*Rblock);
  if (t_gsum) *t_gsum += usecond();
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

  GridBase *grid = wi[0].Grid();
  
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
  thread_for(r,MFrvol,{
    lvSum[r] = Zero();
  });

  Vector<scalar_type > lsSum(MFlvol);             
  thread_for(r,MFlvol,{
    lsSum[r]=scalar_type(0.0);
  });

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];

  thread_for(r,rd,{

    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){

	int ss= so+n*stride+b;

	for(int i=0;i<Lblock;i++){

	  autoView(wi_v,wi[i],CpuRead);
	  auto w = conjugate(wi_v[ss]);
	  if (g5) {
	    w()(2)(0) = - w()(2)(0);
	    w()(2)(1) = - w()(2)(1);
	    w()(2)(2) = - w()(2)(2);
	    w()(3)(0) = - w()(3)(0);
	    w()(3)(1) = - w()(3)(1);
	    w()(3)(2) = - w()(3)(2);
	  }
	  for(int j=0;j<Rblock;j++){
	    
	    autoView(vj_v,vj[j],CpuRead);
	    auto v  = vj_v[ss];
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
  });

  // Sum across simd lanes in the plane, breaking out orthog dir.
  thread_for(rt,rd,{

      Coordinate icoor(nd);
    iScalar<vector_type> temp; 
    ExtractBuffer<iScalar<scalar_type> > extracted(Nsimd);               

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
  });

  assert(mat.dimension(0) == Nt);
  // ld loop and local only??
  int pd = grid->_processors[orthogdim];
  int pc = grid->_processor_coor[orthogdim];
  thread_for_collapse(2,lt,ld,{
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
  });

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

  GridBase *grid = wi[0].Grid();
  
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
  thread_for(r,MFrvol,{
    lvSum[r] = Zero();
  });

  Vector<scalar_type > lsSum(MFlvol);             
  thread_for(r,MFlvol,{
    lsSum[r]=scalar_type(0.0);
  });

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];

  thread_for(r,rd,{

    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){

	int ss= so+n*stride+b;

	for(int i=0;i<Lblock;i++){

	  autoView(wi_v,wi[i],CpuRead);
	  auto w = conjugate(wi_v[ss]);

	  for(int j=0;j<Rblock;j++){

	    autoView(vj_v,vj[j],CpuRead);
	    auto v = vj_v[ss];

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
	      autoView(mom_v,mom[m],CpuRead);
	      auto phase = mom_v[ss];
	      mac(&lvSum[idx],&vv,&phase()()());
	    }
	  }
	}
      }
    }
  });


  // Sum across simd lanes in the plane, breaking out orthog dir.
  thread_for(rt,rd,{

    Coordinate icoor(nd);
    iScalar<vector_type> temp; 
    ExtractBuffer<iScalar<scalar_type> > extracted(Nsimd);               

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
  });

  assert(mat.dimension(0) == Nmom);
  assert(mat.dimension(1) == Nt);
 
  int pd = grid->_processors[orthogdim];
  int pc = grid->_processor_coor[orthogdim];
  thread_for_collapse(2,lt,ld,{
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
  });

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

// "A-slash" field w_i(x)^dag * i * A_mu * gamma_mu * v_j(x)
//
// With:
//
// B_0 = A_0 + i A_1
// B_1 = A_2 + i A_3
// 
// then in spin space
// 
//                 ( 0          0          -conj(B_1) -B_0 )
// i * A_mu g_mu = ( 0          0          -conj(B_0)  B_1 )
//                 ( B_1        B_0        0          0    )
//                 ( conj(B_0)  -conj(B_1) 0          0    )
template <class FImpl>
template <typename TensorType>
void A2Autils<FImpl>::AslashField(TensorType &mat, 
          const FermionField *lhs_wi,
          const FermionField *rhs_vj,
          const std::vector<ComplexField> &emB0,
          const std::vector<ComplexField> &emB1,
          int orthogdim, double *t_kernel, double *t_gsum) 
{
    typedef typename FermionField::vector_object vobj;
    typedef typename vobj::scalar_object         sobj;
    typedef typename vobj::scalar_type           scalar_type;
    typedef typename vobj::vector_type           vector_type;

    typedef iSpinMatrix<vector_type> SpinMatrix_v;
    typedef iSpinMatrix<scalar_type> SpinMatrix_s;
    typedef iSinglet<vector_type>    Singlet_v;
    typedef iSinglet<scalar_type>    Singlet_s;
    
    int Lblock = mat.dimension(3); 
    int Rblock = mat.dimension(4);

    GridBase *grid = lhs_wi[0].Grid();
    
    const int    Nd = grid->_ndimension;
    const int Nsimd = grid->Nsimd();

    int Nt  = grid->GlobalDimensions()[orthogdim];
    int Nem = emB0.size();
    assert(emB1.size() == Nem);

    int fd=grid->_fdimensions[orthogdim];
    int ld=grid->_ldimensions[orthogdim];
    int rd=grid->_rdimensions[orthogdim];

    // will locally sum vectors first
    // sum across these down to scalars
    // splitting the SIMD
    int MFrvol = rd*Lblock*Rblock*Nem;
    int MFlvol = ld*Lblock*Rblock*Nem;

    Vector<vector_type> lvSum(MFrvol);
    thread_for(r,MFrvol,
    {
      lvSum[r] = Zero();
    });

    Vector<scalar_type> lsSum(MFlvol);             
    thread_for(r,MFlvol,
    {
        lsSum[r] = scalar_type(0.0);
    });

    int e1=    grid->_slice_nblock[orthogdim];
    int e2=    grid->_slice_block [orthogdim];
    int stride=grid->_slice_stride[orthogdim];

    // Nested parallelism would be ok
    // Wasting cores here. Test case r
    if (t_kernel) *t_kernel = -usecond();
    thread_for(r,rd,
    {
        int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

        for(int n=0;n<e1;n++)
        for(int b=0;b<e2;b++)
        {
            int ss= so+n*stride+b;

            for(int i=0;i<Lblock;i++)
            {
  	        autoView(wi_v,lhs_wi[i],CpuRead);
                auto left = conjugate(wi_v[ss]);

                for(int j=0;j<Rblock;j++)
                {
                    SpinMatrix_v vv;
		    autoView(vj_v,rhs_vj[j],CpuRead);
                    auto right = vj_v[ss];

                    for(int s1=0;s1<Ns;s1++)
                    for(int s2=0;s2<Ns;s2++)
                    {
		          vv()(s1,s2)() = left()(s2)(0) * right()(s1)(0)
                                        + left()(s2)(1) * right()(s1)(1)
                                        + left()(s2)(2) * right()(s1)(2);
                    }
                    
                    // After getting the sitewise product do the mom phase loop
                    int base = Nem*i+Nem*Lblock*j+Nem*Lblock*Rblock*r;

                    for ( int m=0;m<Nem;m++)
                    {
  		        autoView(emB0_v,emB0[m],CpuRead);
		        autoView(emB1_v,emB1[m],CpuRead);
                        int idx  = m+base;
                        auto b0  = emB0_v[ss];
                        auto b1  = emB1_v[ss];
                        auto cb0 = conjugate(b0);
                        auto cb1 = conjugate(b1);

                        lvSum[idx] += - vv()(3,0)()*b0()()()  - vv()(2,0)()*cb1()()()
                                      + vv()(3,1)()*b1()()()  - vv()(2,1)()*cb0()()()
                                      + vv()(0,2)()*b1()()()  + vv()(1,2)()*b0()()()
                                      + vv()(0,3)()*cb0()()() - vv()(1,3)()*cb1()()();
                    }
                }
            }
        }
    });

    // Sum across simd lanes in the plane, breaking out orthog dir.
    thread_for(rt,rd,
    {
        Coordinate icoor(Nd);
        ExtractBuffer<scalar_type> extracted(Nsimd);               

        for(int i=0;i<Lblock;i++)
        for(int j=0;j<Rblock;j++)
        for(int m=0;m<Nem;m++)
        {

            int ij_rdx = m+Nem*i+Nem*Lblock*j+Nem*Lblock*Rblock*rt;

            extract<vector_type,scalar_type>(lvSum[ij_rdx],extracted);
            for(int idx=0;idx<Nsimd;idx++)
            {
                grid->iCoorFromIindex(icoor,idx);

                int ldx    = rt+icoor[orthogdim]*rd;
                int ij_ldx = m+Nem*i+Nem*Lblock*j+Nem*Lblock*Rblock*ldx;

                lsSum[ij_ldx]=lsSum[ij_ldx]+extracted[idx];
            }
        }
    });
    if (t_kernel) *t_kernel += usecond();

    // ld loop and local only??
    int pd = grid->_processors[orthogdim];
    int pc = grid->_processor_coor[orthogdim];
    thread_for_collapse(2,lt,ld,
    {
        for(int pt=0;pt<pd;pt++)
        {
            int t = lt + pt*ld;
            if (pt == pc)
            {
                for(int i=0;i<Lblock;i++)
                for(int j=0;j<Rblock;j++)
                for(int m=0;m<Nem;m++)
                {
                    int ij_dx = m+Nem*i + Nem*Lblock * j + Nem*Lblock * Rblock * lt;

                    mat(m,0,t,i,j) = lsSum[ij_dx];
                }
            } 
            else 
            { 
                const scalar_type zz(0.0);

                for(int i=0;i<Lblock;i++)
                for(int j=0;j<Rblock;j++)
                for(int m=0;m<Nem;m++)
                {
                    mat(m,0,t,i,j) = zz;
                }
            }
        }
    });
    if (t_gsum) *t_gsum = -usecond();
    grid->GlobalSumVector(&mat(0,0,0,0,0),Nem*Nt*Lblock*Rblock);
    if (t_gsum) *t_gsum += usecond();
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

template <class FImpl>
template <typename TensorType>
typename std::enable_if<(std::is_same<Eigen::Tensor<ComplexD,3>, TensorType>::value ||
                         std::is_same<Eigen::TensorMap<Eigen::Tensor<Complex, 3, Eigen::RowMajor>>, TensorType>::value),
                         void>::type
A2Autils<FImpl>::ContractWWVV(std::vector<PropagatorField> &WWVV,
				   const TensorType &WW_sd,
				   const FermionField *vs,
				   const FermionField *vd)
{
  GridBase *grid = vs[0].Grid();

  //  int nd    = grid->_ndimension;
  int Nsimd = grid->Nsimd();
  int N_t   = WW_sd.dimension(0);
  int N_s = WW_sd.dimension(1); 
  int N_d = WW_sd.dimension(2);

  int d_unroll = 32;// Empirical optimisation

  for(int t=0;t<N_t;t++){
    WWVV[t] = Zero();
  }

  thread_for(ss,grid->oSites(),{
    for(int d_o=0;d_o<N_d;d_o+=d_unroll){
      for(int t=0;t<N_t;t++){
      for(int s=0;s<N_s;s++){
	autoView(vs_v,vs[s],CpuRead);
	auto tmp1 = vs_v[ss];
	vobj tmp2 = Zero();
	vobj tmp3 = Zero();
	for(int d=d_o;d<MIN(d_o+d_unroll,N_d);d++){
	  autoView(vd_v,vd[d],CpuRead);
	  Scalar_v coeff = WW_sd(t,s,d);
	  tmp3 = conjugate(vd_v[ss]);
	  mac(&tmp2, &coeff, &tmp3);
	}

	//////////////////////////
	// Fast outer product of tmp1 with a sum of terms suppressed by d_unroll
	//////////////////////////
	OuterProductWWVV(WWVV[t], tmp1, tmp2, Ns, ss);

      }}
    }
  });
}

template <class FImpl>
template <typename TensorType>
typename std::enable_if<!(std::is_same<Eigen::Tensor<ComplexD, 3>, TensorType>::value ||
                          std::is_same<Eigen::TensorMap<Eigen::Tensor<Complex, 3, Eigen::RowMajor>>, TensorType>::value),
                          void>::type
A2Autils<FImpl>::ContractWWVV(std::vector<PropagatorField> &WWVV,
                              const TensorType &WW_sd,
                              const FermionField *vs,
                              const FermionField *vd)
{
  GridBase *grid = vs[0].Grid();

  int nd    = grid->_ndimension;
  int Nsimd = grid->Nsimd();
  int N_t = WW_sd.dimensions()[0];
  int N_s = WW_sd.dimensions()[1];
  int N_d = WW_sd.dimensions()[2];

  int d_unroll = 32;// Empirical optimisation

  Eigen::Matrix<Complex, -1, -1, Eigen::RowMajor> buf;

  for(int t=0;t<N_t;t++){
    WWVV[t] = Zero();
  }

  for (int t = 0; t < N_t; t++){
    std::cout << GridLogMessage << "Contraction t = " << t << std::endl;
    buf = WW_sd[t];
    thread_for(ss,grid->oSites(),{
      for(int d_o=0;d_o<N_d;d_o+=d_unroll){
        for(int s=0;s<N_s;s++){
	  autoView(vs_v,vs[s],CpuRead);
	  auto tmp1 = vs_v[ss];
	  vobj tmp2 = Zero();
	  vobj tmp3 = Zero();
	  for(int d=d_o;d<MIN(d_o+d_unroll,N_d);d++){
	    autoView(vd_v,vd[d],CpuRead);
	    Scalar_v coeff = buf(s,d);
	    tmp3 = conjugate(vd_v[ss]);
	    mac(&tmp2, &coeff, &tmp3);
	  }
	  //////////////////////////
	  // Fast outer product of tmp1 with a sum of terms suppressed by d_unroll
	  //////////////////////////
	  OuterProductWWVV(WWVV[t], tmp1, tmp2, Ns, ss);
      }}
    });
  }
}

template <class FImpl>
inline void A2Autils<FImpl>::OuterProductWWVV(PropagatorField &WWVV,
                                             const vobj &lhs,
                                             const vobj &rhs,
                                             const int Ns, const int ss)
{
  autoView(WWVV_v,WWVV,CpuWrite);
  for (int s1 = 0; s1 < Ns; s1++){
    for (int s2 = 0; s2 < Ns; s2++){
      WWVV_v[ss]()(s1,s2)(0, 0) += lhs()(s1)(0) * rhs()(s2)(0);
      WWVV_v[ss]()(s1,s2)(0, 1) += lhs()(s1)(0) * rhs()(s2)(1);
      WWVV_v[ss]()(s1,s2)(0, 2) += lhs()(s1)(0) * rhs()(s2)(2);
      WWVV_v[ss]()(s1,s2)(1, 0) += lhs()(s1)(1) * rhs()(s2)(0);
      WWVV_v[ss]()(s1,s2)(1, 1) += lhs()(s1)(1) * rhs()(s2)(1);
      WWVV_v[ss]()(s1,s2)(1, 2) += lhs()(s1)(1) * rhs()(s2)(2);
      WWVV_v[ss]()(s1,s2)(2, 0) += lhs()(s1)(2) * rhs()(s2)(0);
      WWVV_v[ss]()(s1,s2)(2, 1) += lhs()(s1)(2) * rhs()(s2)(1);
      WWVV_v[ss]()(s1,s2)(2, 2) += lhs()(s1)(2) * rhs()(s2)(2);
    }
  }
}

template<class FImpl>
void A2Autils<FImpl>::ContractFourQuarkColourDiagonal(const PropagatorField &WWVV0,
						      const PropagatorField &WWVV1,
						      const std::vector<Gamma> &gamma0,
						      const std::vector<Gamma> &gamma1,
						      ComplexField &O_trtr,
						      ComplexField &O_fig8)
{
  assert(gamma0.size()==gamma1.size());
  int Ng = gamma0.size();

  GridBase *grid = WWVV0.Grid();

  autoView(WWVV0_v , WWVV0,CpuRead);
  autoView(WWVV1_v , WWVV1,CpuRead);
  autoView(O_trtr_v, O_trtr,CpuWrite);
  autoView(O_fig8_v, O_fig8,CpuWrite);
  thread_for(ss,grid->oSites(),{

    typedef typename ComplexField::vector_object vobj;

    vobj v_trtr;
    vobj v_fig8;

    auto VV0 = WWVV0_v[ss];
    auto VV1 = WWVV1_v[ss];
    
    for(int g=0;g<Ng;g++){

      v_trtr = trace(VV0 * gamma0[g])* trace(VV1*gamma1[g]);
      v_fig8 = trace(VV0 * gamma0[g] * VV1 * gamma1[g]);

      if ( g==0 ) {
	O_trtr_v[ss] = v_trtr; 
	O_fig8_v[ss] = v_fig8;
      } else { 
	O_trtr_v[ss]+= v_trtr; 
	O_fig8_v[ss]+= v_fig8;
      }
      
    }
  });
}

template<class FImpl>
void A2Autils<FImpl>::ContractFourQuarkColourMix(const PropagatorField &WWVV0,
						 const PropagatorField &WWVV1,
						 const std::vector<Gamma> &gamma0,
						 const std::vector<Gamma> &gamma1,
						 ComplexField &O_trtr,
						 ComplexField &O_fig8)
{
  assert(gamma0.size()==gamma1.size());
  int Ng = gamma0.size();

  GridBase *grid = WWVV0.Grid();

  autoView( WWVV0_v , WWVV0,CpuRead);
  autoView( WWVV1_v , WWVV1,CpuRead);
  autoView( O_trtr_v, O_trtr,CpuWrite);
  autoView( O_fig8_v, O_fig8,CpuWrite);

  thread_for(ss,grid->oSites(),{

    typedef typename ComplexField::vector_object vobj;

    auto VV0 = WWVV0_v[ss];
    auto VV1 = WWVV1_v[ss];
    
    for(int g=0;g<Ng;g++){

      auto VV0G = VV0 * gamma0[g];  // Spin multiply
      auto VV1G = VV1 * gamma1[g];

      vobj v_trtr=Zero();
      vobj v_fig8=Zero();

      /////////////////////////////////////////
      // Colour mixed
      /////////////////////////////////////////
      // _                 _
      // s_sa G_st d_tb    s_s'b  G_s't' d_t'a
      // 
      //                                         
      // Contracted with prop factor (VV0)_sd,ab (VV1)_s'd',ba
      //
      // Wick1 [ spin TR TR ]
      //
      //    (VV0*G0)_ss,ba .  (VV1*G1)_tt,ab
       //
      // Wick2 [ spin fig8 ]
      //
      //    (VV0*G0)_st,aa (VV1*G1)_ts,bb
      // 
      /////////////////////////////////////////

      for(int a=0;a<Nc;a++){
      for(int b=0;b<Nc;b++){
      for(int s=0;s<Ns;s++){
      for(int t=0;t<Ns;t++){
	// Mixed traces
	v_trtr()()() += VV0G()(s,s)(a,b)*VV1G()(t,t)(b,a); // Was the fig8 before Fierzing
	v_fig8()()() += VV0G()(s,t)(a,a)*VV1G()(t,s)(b,b); // Was the trtr before Fierzing

	/*
	 * CHECKS -- use Fierz identities as a strong test, 4 Oct 2018.
	 * 
BagMix [8,0]  fig8 (21.5596,-3.83908e-17) trtr (0.064326,2.51001e-17) // Fierz        -1 0   0 0 0    
BagMix [8,1]  fig8 (-1346.99,1.2481e-16) trtr (34.2501,-3.36935e-17)  //               0 0   2 0 0 
BagMix [8,2]  fig8 (13.7536,-6.04625e-19) trtr (-215.542,3.24326e-17) //               0 1/2 0 0 0
BagMix [8,3]  fig8 (555.878,-7.39942e-17) trtr (463.82,-4.73909e-17)  //               0 0   0 1/2 -1/2  
BagMix [8,4]  fig8 (-1602.48,9.08511e-17) trtr (-936.302,1.14156e-16) //               0 0   0 -3/2 -1/2

Bag [8,0]  fig8 (-0.064326,1.06281e-17) trtr (-21.5596,1.06051e-17)   
Bag [8,2]  fig8 (17.125,-3.40959e-17) trtr (-673.493,7.68134e-17)   
Bag [8,1]  fig8 (-431.084,2.76423e-17) trtr (27.5073,-5.76967e-18)    /////////// TR TR                           FIG8
Bag [8,3]  fig8 (700.061,-1.14925e-16) trtr (1079.18,-1.35476e-16)    //    555.878 = 0.5(1079.18+32.5776) ;   463.82     =0.5(700.061+227.58)
Bag [8,4]  fig8 (-227.58,3.58808e-17) trtr (-32.5776,1.83286e-17)     //  - 1602.48 = - 1.5*1079.18 + .5* 32.5776; 936.302=-1.5* 700+0.5*227
	*/
	
	//Unmixed debug check consistency
	//	v_trtr()()() += VV0G()(s,s)(a,a)*VV1G()(t,t)(b,b);
	//	v_fig8()()() += VV0G()(s,t)(a,b)*VV1G()(t,s)(b,a);
      }}}}

      if ( g==0 ) {
	O_trtr_v[ss] = v_trtr; 
	O_fig8_v[ss] = v_fig8;
      } else { 
	O_trtr_v[ss]+= v_trtr; 
	O_fig8_v[ss]+= v_fig8;
      }
      
    }
  });
}

#ifdef DELTA_F_EQ_2
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Perhaps this should move out of the utils and into Hadrons module
// Now makes use of the primitives above and doesn't touch inside 
// the lattice structures.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class FImpl>
void A2Autils<FImpl>::DeltaFeq2(int dt_min,int dt_max,
				Eigen::Tensor<ComplexD,2> &dF2_fig8,
				Eigen::Tensor<ComplexD,2> &dF2_trtr,
				Eigen::Tensor<ComplexD,2> &dF2_fig8_mix,
				Eigen::Tensor<ComplexD,2> &dF2_trtr_mix,
				Eigen::Tensor<ComplexD,1> &denom_A0,
				Eigen::Tensor<ComplexD,1> &denom_P,
				Eigen::Tensor<ComplexD,3> &WW_sd, 
				const FermionField *vs,
				const FermionField *vd,
				int orthogdim)
{
  GridBase *grid = vs[0].Grid();

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
  int Ng=5;
  dF2_trtr.resize(N_t,Ng);
  dF2_fig8.resize(N_t,Ng);
  dF2_trtr_mix.resize(N_t,Ng);
  dF2_fig8_mix.resize(N_t,Ng);
  denom_A0.resize(N_t);
  denom_P.resize(N_t);
  for(int t=0;t<N_t;t++){
    for(int g=0;g<Ng;g++) dF2_trtr(t,g)= ComplexD(0.0);
    for(int g=0;g<Ng;g++) dF2_fig8(t,g)= ComplexD(0.0);
    for(int g=0;g<Ng;g++) dF2_trtr_mix(t,g)= ComplexD(0.0);
    for(int g=0;g<Ng;g++) dF2_fig8_mix(t,g)= ComplexD(0.0);
    denom_A0(t) =ComplexD(0.0);
    denom_P(t) =ComplexD(0.0);
  }

  ComplexField D0(grid);   D0 = Zero(); // <P|A0> correlator from each wall
  ComplexField D1(grid);   D1 = Zero();

  ComplexField O1_trtr(grid);  O1_trtr = Zero();
  ComplexField O2_trtr(grid);  O2_trtr = Zero();
  ComplexField O3_trtr(grid);  O3_trtr = Zero();
  ComplexField O4_trtr(grid);  O4_trtr = Zero();
  ComplexField O5_trtr(grid);  O5_trtr = Zero();

  ComplexField O1_fig8(grid);  O1_fig8 = Zero();
  ComplexField O2_fig8(grid);  O2_fig8 = Zero();
  ComplexField O3_fig8(grid);  O3_fig8 = Zero();
  ComplexField O4_fig8(grid);  O4_fig8 = Zero();
  ComplexField O5_fig8(grid);  O5_fig8 = Zero();

  ComplexField VV_trtr(grid);  VV_trtr = Zero();
  ComplexField AA_trtr(grid);  AA_trtr = Zero();
  ComplexField SS_trtr(grid);  SS_trtr = Zero();
  ComplexField PP_trtr(grid);  PP_trtr = Zero();
  ComplexField TT_trtr(grid);  TT_trtr = Zero();

  ComplexField VV_fig8(grid);  VV_fig8 = Zero();
  ComplexField AA_fig8(grid);  AA_fig8 = Zero();
  ComplexField SS_fig8(grid);  SS_fig8 = Zero();
  ComplexField PP_fig8(grid);  PP_fig8 = Zero();
  ComplexField TT_fig8(grid);  TT_fig8 = Zero();

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
      std::vector<Gamma> VV({VX,VY,VZ,VT});
      std::vector<Gamma> AA({AX,AY,AZ,AT});
      std::vector<Gamma> SS({S});
      std::vector<Gamma> PP({P});
      std::vector<Gamma> TT({T0,T1,T2,T3,T4,T5});
      std::vector<Gamma> A0({AT});

      ContractFourQuarkColourDiagonal(WWVV[t0],	WWVV[t1],VV,VV,VV_trtr,VV_fig8); // VV
      ContractFourQuarkColourDiagonal(WWVV[t0],	WWVV[t1],AA,AA,AA_trtr,AA_fig8); // AA
      ContractFourQuarkColourDiagonal(WWVV[t0],	WWVV[t1],SS,SS,SS_trtr,SS_fig8); // SS
      ContractFourQuarkColourDiagonal(WWVV[t0],	WWVV[t1],PP,PP,PP_trtr,PP_fig8); // PP
      ContractFourQuarkColourDiagonal(WWVV[t0],	WWVV[t1],TT,TT,TT_trtr,TT_fig8); // TT

      O1_trtr = VV_trtr+AA_trtr;      O2_trtr = VV_trtr-AA_trtr;  // VV+AA,VV-AA
      O1_fig8 = VV_fig8+AA_fig8;      O2_fig8 = VV_fig8-AA_fig8;  

      O3_trtr = SS_trtr-PP_trtr;      O4_trtr = SS_trtr+PP_trtr;  // SS+PP,SS-PP
      O3_fig8 = SS_fig8-PP_fig8;      O4_fig8 = SS_fig8+PP_fig8;  

      O5_trtr = TT_trtr;
      O5_fig8 = TT_fig8;

      sliceSum(O1_trtr,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_trtr(t,0)+= 2.0*C1[(t+t0)%N_t]()()()/vol;
      sliceSum(O2_trtr,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_trtr(t,1)+= 2.0*C1[(t+t0)%N_t]()()()/vol;
      sliceSum(O3_trtr,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_trtr(t,2)+= 2.0*C1[(t+t0)%N_t]()()()/vol;
      sliceSum(O4_trtr,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_trtr(t,3)+= 2.0*C1[(t+t0)%N_t]()()()/vol;
      sliceSum(O5_trtr,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_trtr(t,4)+= 2.0*C1[(t+t0)%N_t]()()()/vol;

      sliceSum(O1_fig8,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_fig8(t,0)+= 2.0*C1[(t+t0)%N_t]()()()/vol;
      sliceSum(O2_fig8,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_fig8(t,1)+= 2.0*C1[(t+t0)%N_t]()()()/vol;
      sliceSum(O3_fig8,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_fig8(t,2)+= 2.0*C1[(t+t0)%N_t]()()()/vol;
      sliceSum(O4_fig8,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_fig8(t,3)+= 2.0*C1[(t+t0)%N_t]()()()/vol;
      sliceSum(O5_fig8,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_fig8(t,4)+= 2.0*C1[(t+t0)%N_t]()()()/vol;

      ContractFourQuarkColourDiagonal(WWVV[t0],	WWVV[t1],A0,A0,AA_trtr,AA_fig8); // A0 insertion

      sliceSum(AA_trtr,C1, orthogdim);
      sliceSum(PP_trtr,C2, orthogdim);

      for(int t=0;t<N_t;t++){
	denom_A0(t)+=C1[(t+t0)%N_t]()()()/vol;
	denom_P(t) +=C2[(t+t0)%N_t]()()()/vol;
      }

      ///////////////////////////////////////////////////////////////////////////
      // Colour mixed contractions
      ///////////////////////////////////////////////////////////////////////////

      ContractFourQuarkColourMix(WWVV[t0],	WWVV[t1],VV,VV,VV_trtr,VV_fig8); // VV
      ContractFourQuarkColourMix(WWVV[t0],	WWVV[t1],AA,AA,AA_trtr,AA_fig8); // AA
      ContractFourQuarkColourMix(WWVV[t0],	WWVV[t1],SS,SS,SS_trtr,SS_fig8); // SS
      ContractFourQuarkColourMix(WWVV[t0],	WWVV[t1],PP,PP,PP_trtr,PP_fig8); // PP
      ContractFourQuarkColourMix(WWVV[t0],	WWVV[t1],TT,TT,TT_trtr,TT_fig8); // TT

      O1_trtr = VV_trtr+AA_trtr;      O2_trtr = VV_trtr-AA_trtr;  // VV+AA,VV-AA
      O1_fig8 = VV_fig8+AA_fig8;      O2_fig8 = VV_fig8-AA_fig8;  

      O3_trtr = SS_trtr-PP_trtr;      O4_trtr = SS_trtr+PP_trtr;  // SS+PP,SS-PP
      O3_fig8 = SS_fig8-PP_fig8;      O4_fig8 = SS_fig8+PP_fig8;  

      O5_trtr = TT_trtr;
      O5_fig8 = TT_fig8;

      sliceSum(O1_trtr,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_trtr_mix(t,0)+= 2.0*C1[(t+t0)%N_t]()()()/vol;
      sliceSum(O2_trtr,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_trtr_mix(t,1)+= 2.0*C1[(t+t0)%N_t]()()()/vol;
      sliceSum(O3_trtr,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_trtr_mix(t,2)+= 2.0*C1[(t+t0)%N_t]()()()/vol;
      sliceSum(O4_trtr,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_trtr_mix(t,3)+= 2.0*C1[(t+t0)%N_t]()()()/vol;
      sliceSum(O5_trtr,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_trtr_mix(t,4)+= 2.0*C1[(t+t0)%N_t]()()()/vol;

      sliceSum(O1_fig8,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_fig8_mix(t,0)+= 2.0*C1[(t+t0)%N_t]()()()/vol;
      sliceSum(O2_fig8,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_fig8_mix(t,1)+= 2.0*C1[(t+t0)%N_t]()()()/vol;
      sliceSum(O3_fig8,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_fig8_mix(t,2)+= 2.0*C1[(t+t0)%N_t]()()()/vol;
      sliceSum(O4_fig8,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_fig8_mix(t,3)+= 2.0*C1[(t+t0)%N_t]()()()/vol;
      sliceSum(O5_fig8,C1, orthogdim); for(int t=0;t<N_t;t++) dF2_fig8_mix(t,4)+= 2.0*C1[(t+t0)%N_t]()()()/vol;

    }
  }
  t_contr +=usecond();
  
  t_tot+=usecond();
  double million=1.0e6;
  LOG(Message) << "Computing A2A DeltaF=2 graph t_tot      " << t_tot      /million << " s "<< std::endl;
  LOG(Message) << "Computing A2A DeltaF=2 graph t_outer    " << t_outer    /million << " s "<< std::endl;
  LOG(Message) << "Computing A2A DeltaF=2 graph t_contr    " << t_contr    /million << " s "<< std::endl;
}
#endif 

template <class FImpl>
template <typename TensorType>
void A2Autils<FImpl>::StagMesonFieldMILC(TensorType &mat,
                                     const FermionField *lhs_wi,
                                     const FermionField *rhs_vj,
                                     std::vector<Gamma::Algebra> gammas,
                                     const std::vector<ComplexField > &mom,
                                     int orthogdim, double *t_kernel, double *t_gsum)
{
  typedef decltype(coalescedRead(Scalar_v())) calcScalar;
  typedef decltype(coalescedRead(vobj())) calcSpinor;
  

  int Lblock = mat.dimension(3);
  int Rblock = mat.dimension(4);

  int sizeL = Lblock;
  int sizeR = Rblock;

  GridBase *grid  = mom[0].Grid();    
  GridBase *hgrid = nullptr;

  int cb;
  Vector<Integer> hcoor;

  int checkerL = lhs_wi[0].Grid()->_isCheckerBoarded, 
      checkerR = rhs_vj[0].Grid()->_isCheckerBoarded;

  if (checkerL) {
    sizeL = Lblock/2;
    hgrid = lhs_wi[0].Grid();
    cb = lhs_wi[0].Checkerboard();

    if (checkerR) {
      assert(cb == rhs_vj[0].Checkerboard());
    }
  }
  if(checkerR) {
    sizeR = Rblock/2;
    hgrid = rhs_vj[0].Grid();
    cb = rhs_vj[0].Checkerboard();
  }

  if (checkerR || checkerL)
    hcoor.resize(hgrid->oSites(),0);

  const int    Nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();
  
  int Nt     = grid->GlobalDimensions()[orthogdim];
  int Ngamma = gammas.size();
  int Nmom   = mom.size();
  
  // The full length of the lattice in the "time" dimension
  int fd=grid->_fdimensions[orthogdim];
  // The local number of time slices held by this processor
  int ld=grid->_ldimensions[orthogdim];
  // The reduced number of time slices after ignoring SIMD vector lanes
  int rd=grid->_rdimensions[orthogdim];
  
  int MFrvol = rd*Lblock*Rblock*Nmom*Ngamma;
  int MFlvol = ld*Lblock*Rblock*Nmom*Ngamma;
  
  Vector<Scalar_v> lvSum(MFrvol);
  Vector<Scalar_s> lsSum(MFlvol);

  thread_for(r, MFrvol,{
    lvSum[r] = Zero();
  });  
  thread_for(r, MFlvol,{
    lsSum[r]=scalar_type(0.0);
  });

  int nBlocks = grid->_slice_nblock[orthogdim];
  int vecsPerSlicePerBlock = grid->_slice_block[orthogdim];
  int blockStride  = grid->_slice_stride[orthogdim];
  int rtStride   = grid->_ostride[orthogdim];

  if (checkerR || checkerL) {
    // Global sum fails at end if checkerboard dim is orthogdim
    assert(hgrid->CheckerBoarded(orthogdim) != 1);

    nBlocks = hgrid->_slice_nblock[orthogdim];
    vecsPerSlicePerBlock = hgrid->_slice_block[orthogdim];
    blockStride  = hgrid->_slice_stride[orthogdim];
    rtStride   = hgrid->_ostride[orthogdim];

    thread_for(ss,grid->oSites(),{
      int cbos;
      Coordinate coor;

      grid->oCoorFromOindex(coor,ss);
      cbos=hgrid->CheckerBoard(coor);
        
      if (cbos==cb) {
        int ssh=hgrid->oIndex(coor);
        hcoor[ssh]=ss;
      }
    });
  }


  // potentially wasting cores here if local time extent too small
  if (t_kernel) *t_kernel = -usecond();
  
  // based on phases in Grid/qcd/action/fermion/FermionOperatorImpl.h
  // Pretty cute implementation, if I may say so myself (!) (-PAB)
  // Staggered Phases for local, taste non-singlet meson operators.
  // See Degrand and Detar, Ch 11.2
  
  Lattice<iScalar<vInteger> > x(grid); LatticeCoordinate(x,0);
  Lattice<iScalar<vInteger> > y(grid); LatticeCoordinate(y,1);
  Lattice<iScalar<vInteger> > z(grid); LatticeCoordinate(z,2);
  Lattice<iScalar<vInteger> > t(grid); LatticeCoordinate(t,3);
  
  Lattice<iScalar<vInteger> > lin_x(grid); lin_x=y+z+t;
  Lattice<iScalar<vInteger> > lin_y(grid); lin_y=x+z+t;
  Lattice<iScalar<vInteger> > lin_z(grid); lin_z=x+y+t;
  Lattice<iScalar<vInteger> > lin_5(grid); lin_5=x+y+z+t;

  std::vector<ComplexField> stagphase(Ngamma,grid);   
  for (int mu = 0; mu < Ngamma; mu++) {

      stagphase[mu]=1.0;        
      if ( gammas[mu] == Gamma::Algebra::Gamma5 ) stagphase[mu] = where( mod(lin_5,2)==(Integer)0, stagphase[mu],-stagphase[mu]);
      else if ( gammas[mu] == Gamma::Algebra::GammaX ) stagphase[mu] = where( mod(lin_x,2)==(Integer)0, stagphase[mu],-stagphase[mu]);
      else if ( gammas[mu] == Gamma::Algebra::GammaY ) stagphase[mu] = where( mod(lin_y,2)==(Integer)0, stagphase[mu],-stagphase[mu]);
      else if ( gammas[mu] == Gamma::Algebra::GammaZ ) stagphase[mu] = where( mod(lin_z,2)==(Integer)0, stagphase[mu],-stagphase[mu]);
      else {
          std::cout << gammas[mu] << " not implemented for staggered fermion meson field" << std::endl;
          assert(0);
      }
  }

  // Setup lists of pointers to share with accelerators
  typedef LatticeView<vobj> FermView;
  typedef LatticeView<typename ComplexField::vector_object> CView;

  Vector<FermView> AcceleratorViewContainerV, AcceleratorViewContainerW;
  Vector<CView> AcceleratorViewContainerStag, AcceleratorViewContainerMom;
  Vector<Coordinate> icoorContainer(Nsimd,Coordinate(Nd));
  for(int p=0;p<sizeL;p++) {
    AcceleratorViewContainerW.push_back(lhs_wi[p].View(AcceleratorRead));
  }
  for(int p=0;p<sizeR;p++) {
    AcceleratorViewContainerV.push_back(rhs_vj[p].View(AcceleratorRead));
  } 
  for(int p=0;p<Ngamma;p++) {
    AcceleratorViewContainerStag.push_back(stagphase[p].View(AcceleratorRead));
  }
  for(int p=0;p<Nmom;p++) {
    AcceleratorViewContainerMom.push_back(mom[p].View(AcceleratorRead));
  }
  for(int p=0;p<Nsimd;p++) {
    grid->iCoorFromIindex(icoorContainer[p],p);
  }

  // These are hacks to pass data to the device without using std::vector
  FermView *view_pW = & AcceleratorViewContainerW[0];
  FermView *view_pV = & AcceleratorViewContainerV[0];
  CView    *view_pS = & AcceleratorViewContainerStag[0];
  CView    *view_pM = & AcceleratorViewContainerMom[0];
  Integer    *hcoor_p = & hcoor[0];

  Scalar_v *lvSum_p = &lvSum[0];
  Scalar_s *lsSum_p = &lsSum[0];

  Coordinate *icoor_p = &icoorContainer[0];

  // Hack to prefetch fields onto device memory. There's probably a more explicit way to do it?
  accelerator_forNB(index,sizeL+sizeR,1,{
    if (index < sizeL)
      auto a = view_pW[index][0];
    else
      auto a = view_pV[index-sizeL][0];                                                                                                                                                                                
  });  

  accelerator_for2d(l_index,sizeL,r_index,sizeR,Nsimd,{

    // Spawn threads for each time slice on the local processor
    for (int rt=0;rt<rd;rt++) {
      int so=rt*rtStride; // base offset for start of the local plane
      int base = Ngamma*Nmom*(checkerL?2:1)*l_index+Ngamma*Nmom*Lblock*(checkerR?2:1)*r_index+Ngamma*Nmom*Lblock*Rblock*rt;

      calcSpinor left, right;
      calcScalar temp_site, temp_phase, temp_phase_neg, gamma_phase, momentum_phase;
      for(int n=0;n<nBlocks;n++){
        for(int b=0;b<vecsPerSlicePerBlock;b++){
          int ss = so+n*blockStride+b;
          int fullss = ss;

          if (checkerL || checkerR)
            fullss = hcoor_p[ss];

          if (checkerL)
            left  = coalescedRead(view_pW[l_index][ss]);
          else
            left  = coalescedRead(view_pW[l_index][fullss]);

          if (checkerR)
            right = coalescedRead(view_pV[r_index][ss]);
          else
            right = coalescedRead(view_pV[r_index][fullss]);
 
          temp_site  = innerProduct(left,right);

          for (int mu = 0; mu < Ngamma; mu++) {
            gamma_phase = coalescedRead(view_pS[mu][fullss]);
            for ( int m=0;m<Nmom;m++){
              int idx = base+mu*Nmom+m;

              momentum_phase = coalescedRead(view_pM[m][fullss]);

              temp_phase = gamma_phase*momentum_phase;
              temp_phase_neg = -temp_phase;

              // Add contribution to 4 combinations of evec pairs
              for (int lshifts = 0; lshifts < (checkerL?2:1);lshifts ++){
                for (int rshifts = 0; rshifts < (checkerR?2:1);rshifts ++){

                  int shift_idx = idx + Ngamma*Nmom*(lshifts+rshifts*Lblock);
                  auto sum = coalescedRead(lvSum_p[shift_idx]);

                  // The Mdag eigenvectors have a negative odd checkerboard
                  if ( (lshifts != rshifts) && cb == Odd)
                    mac(&sum,&temp_site,&temp_phase_neg);
                  else
                    mac(&sum,&temp_site,&temp_phase);

                  coalescedWrite(lvSum_p[shift_idx],sum);
                }
              }
            }
          }
        }
      }
    }
  });


  accelerator_for(l_index,Lblock,acceleratorThreads(),{
    // Sum across simd lanes in the plane, breaking out orthog dir.
    ExtractBuffer<Scalar_s> extracted(Nsimd);

    for(int r_index=0;r_index<Rblock;r_index++){

      for (int rt=0;rt<rd;rt++) {

        int base = Ngamma*Nmom*l_index+Ngamma*Nmom*Lblock*r_index+Ngamma*Nmom*Lblock*Rblock*rt;
        for (int mu = 0; mu < Ngamma; mu++) {
          for ( int m=0;m<Nmom;m++){
            
            // Final matrix layout is mat[Nmom,Ngamma,Nt,Nw,Nv]
            // Thus, ignoring gamma and time we get the index of lvSum
            int ij_rdx = base+mu*Nmom+m;
            
            extract(lvSum_p[ij_rdx],extracted);
            
            for(int idx=0;idx<Nsimd;idx++){
                
                // Get the *internal* coord that the SIMD lane idx represents
                
                // We are looking to add each element of a SIMD lane to the appropriate time slice in lsSum.
                // ldx is the actual time slice on the lattice
                int ldx    = rt+icoor_p[idx][orthogdim]*rd;
                
                // The offset for lsSum is thus lsSum[timeSlice,Lblock,Rblock,gamma,momentum] in row major
                int ij_ldx = m+mu*Nmom+Ngamma*Nmom*l_index+Ngamma*Nmom*Lblock*r_index+Ngamma*Nmom*Lblock*Rblock*ldx;
                
                lsSum_p[ij_ldx]=lsSum_p[ij_ldx]+extracted[idx];
            }
          }
        }
      }
    }
  });

  assert(mat.dimension(0) == Nmom);
  assert(mat.dimension(1) == Ngamma);
  assert(mat.dimension(2) == Nt);

  int pd = grid->_processors[orthogdim];
  int pc = grid->_processor_coor[orthogdim];

  thread_for_collapse(2,lt,ld,{
    for(int pt=0;pt<pd;pt++){
      int t = lt + pt*ld;
      // Fill mat only with the data you have locally
      if (pt == pc){
        for(int i=0;i<Lblock;i++){
          for(int j=0;j<Rblock;j++){
            for (int mu = 0; mu < Ngamma; mu++) {
              for(int m=0;m<Nmom;m++){
                int ij_dx = m+Nmom*mu+Ngamma*Nmom*i + Ngamma*Nmom*Lblock * j + Ngamma*Nmom*Lblock * Rblock * lt;
                mat(m,mu,t,i,j) = lsSum[ij_dx];
              }
            }
          }
        }
      // Otherwise fill with zeros
      } else {
        const scalar_type zz(0.0);
        for(int i=0;i<Lblock;i++){
          for(int j=0;j<Rblock;j++){
            for (int mu = 0; mu < Ngamma; mu++) {
              for(int m=0;m<Nmom;m++){
                mat(m,mu,t,i,j) = zz;
              }
            }
          }
        }
      }
    }
  });

  if (t_kernel) *t_kernel += usecond();

  // Combine local mat objects from all processors so that it's completely filled in.
  ////////////////////////////////////////////////////////////////////
  // This global sum is taking as much as 50% of time on 16 nodes
  // Vector size is 7 x 16 x 32 x 16 x 16 x sizeof(complex) = 2MB - 60MB depending on volume
  // Healthy size that should suffice
  ////////////////////////////////////////////////////////////////////
  if (t_gsum) *t_gsum = -usecond();
  grid->GlobalSumVector(&mat(0,0,0,0,0),Nmom*Ngamma*Nt*Lblock*Rblock);
  if (t_gsum) *t_gsum += usecond();

  // Free all Lattice Views
  for(int p=0;p<sizeL;p++) AcceleratorViewContainerW[p].ViewClose();
  for(int p=0;p<sizeR;p++) AcceleratorViewContainerV[p].ViewClose();
  for(int p=0;p<Ngamma;p++) AcceleratorViewContainerStag[p].ViewClose();
  for(int p=0;p<Nmom;p++) AcceleratorViewContainerMom[p].ViewClose();
}

// local, taste-nonsignlet staggered mesons
// only the goldstone pion and gamma_{x,y,z} implemented for now.
template <class FImpl>
template <typename TensorType>
void A2Autils<FImpl>::StagMesonField(TensorType &mat,
                                     const FermionField *lhs_wi,
                                     const FermionField *rhs_vj,
                                     std::vector<Gamma::Algebra> gammas,
                                     const std::vector<ComplexField > &mom,
                                     int orthogdim, double *t_kernel, double *t_gsum)
{
    typedef typename FImpl::SiteSpinor vobj;
    
    typedef typename vobj::scalar_object sobj;
    typedef typename vobj::scalar_type scalar_type;
    typedef typename vobj::vector_type vector_type;
    
    typedef iSinglet<vector_type> Singlet_v;
    typedef iSinglet<scalar_type> Singlet_s;
    
    int Lblock = mat.dimension(3);
    int Rblock = mat.dimension(4);
    
    //GridBase *grid = lhs_wi[0]._grid;
    GridBase *grid = lhs_wi[0].Grid();    

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
    
    Vector<Singlet_v > lvSum(MFrvol);
    // thread_for(r, MFrvol,{
    //   lvSum[r] = Zero();
    // });
    //parallel_for (int r = 0; r < MFrvol; r++){
    //    lvSum[r] = zero;
    //}
    
    Vector<Singlet_s > lsSum(MFlvol);
    // thread_for(r, MFlvol,{
    //   lsSum[r]=scalar_type(0.0);
    // });
    //parallel_for (int r = 0; r < MFlvol; r++){
    //    lsSum[r]=scalar_type(0.0);
    //}
    
    int e1=    grid->_slice_nblock[orthogdim];
    int e2=    grid->_slice_block [orthogdim];
    int stride=grid->_slice_stride[orthogdim];
    
    // potentially wasting cores here if local time extent too small
    if (t_kernel) *t_kernel = -usecond();
    
    // based on phases in Grid/qcd/action/fermion/FermionOperatorImpl.h
    // Pretty cute implementation, if I may say so myself (!) (-PAB)
    // Staggered Phases for local, taste non-singlet meson operators.
    // See Degrand and Detar, Ch 11.2
    
    //Lattice<iScalar<vInteger> > coor(grid);
    Lattice<iScalar<vInteger> > x(grid); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(grid); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > z(grid); LatticeCoordinate(z,2);
    Lattice<iScalar<vInteger> > t(grid); LatticeCoordinate(t,3);
    
    Lattice<iScalar<vInteger> > lin_x(grid); lin_x=y+z+t;
    Lattice<iScalar<vInteger> > lin_y(grid); lin_y=x+z+t;
    Lattice<iScalar<vInteger> > lin_z(grid); lin_z=x+y+t;
    Lattice<iScalar<vInteger> > lin_5(grid); lin_5=x+y+z+t;

    for (int mu = 0; mu < Ngamma; mu++) {
        // Re-initialize working variables before starting on a new gamma
        thread_for(r, MFrvol,{
          lvSum[r] = Zero();
        });  
        thread_for(r, MFlvol,{
          lsSum[r]=scalar_type(0.0);
        });
        ComplexField stagphase(grid);   stagphase=1.0;
        
        if ( gammas[mu] == Gamma::Algebra::Gamma5 ) stagphase = where( mod(lin_5,2)==(Integer)0, stagphase,-stagphase);
        else if ( gammas[mu] == Gamma::Algebra::GammaX ) stagphase = where( mod(lin_x,2)==(Integer)0, stagphase,-stagphase);
        else if ( gammas[mu] == Gamma::Algebra::GammaY ) stagphase = where( mod(lin_y,2)==(Integer)0, stagphase,-stagphase);
        else if ( gammas[mu] == Gamma::Algebra::GammaZ ) stagphase = where( mod(lin_z,2)==(Integer)0, stagphase,-stagphase);
        else {
            std::cout << gammas[mu] << " not implemented for staggered fermion meson field" << std::endl;
            assert(0);
        }
        thread_for(r, rd,{
            
            int so=r*grid->_ostride[orthogdim]; // base offset for start of plane
            
            for(int n=0;n<e1;n++){
                for(int b=0;b<e2;b++){
                    
                    int ss= so+n*stride+b;
                    for(int i=0;i<Lblock;i++){
                        
                        //auto left = conjugate(lhs_wi[i]._odata[ss]);
                        //auto wi_v  = lhs_wi[i].View();
    		        autoView(wi_v, lhs_wi[i], CpuRead);
                        auto left = conjugate(wi_v[ss]);

                        for(int j=0;j<Rblock;j++){
                            
                            Singlet_v vv;
                            //auto right = rhs_vj[j]._odata[ss];
                  			    autoView(vj_v,rhs_vj[j],CpuRead);
                            //auto vj_v = rhs_vj[j].View();
                            auto right = vj_v[ss];
 
                            vv()()() = left()()(0) * right()()(0)
                                     + left()()(1) * right()()(1)
                                     + left()()(2) * right()()(2);
                            
                            // After getting the sitewise product do the mom phase loop
                            int base = Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*r;
                            for ( int m=0;m<Nmom;m++){
                                int idx = m+base;
                                //auto phase = mom[m]._odata[ss] * stagphase._odata[ss];
				autoView(mom_v,mom[m],CpuRead);
				// auto mom_v = mom[m].View();
				autoView(stagphase_v,stagphase,CpuRead);
                                //auto stagphase_v = stagphase.View();
                                auto phase = mom_v[ss] * stagphase_v[ss];
                                mac(&lvSum[idx],&vv,&phase);
                            }
                        }
                    }
                }
            }
        });
    
    
        // Sum across simd lanes in the plane, breaking out orthog dir.
        thread_for(rt, rd,{
        //parallel_for(int rt=0;rt<rd;rt++){
            Coordinate icoor(Nd);
            //std::vector<int> icoor(Nd);
            //std::vector<Singlet_s> extracted(Nsimd);
            ExtractBuffer<Singlet_s> extracted(Nsimd);
            
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
        });
        if (t_kernel) *t_kernel += usecond();
        assert(mat.dimension(0) == Nmom);
        assert(mat.dimension(1) == Ngamma);
        assert(mat.dimension(2) == Nt);
        
        // ld loop and local only??
        int pd = grid->_processors[orthogdim];
        int pc = grid->_processor_coor[orthogdim];
        thread_for_collapse(2,lt,ld,{
            
            for(int pt=0;pt<pd;pt++){
                int t = lt + pt*ld;
                if (pt == pc){
                    for(int i=0;i<Lblock;i++){
                        for(int j=0;j<Rblock;j++){
                            for(int m=0;m<Nmom;m++){
                                int ij_dx = m+Nmom*i + Nmom*Lblock * j + Nmom*Lblock * Rblock * lt;
                                mat(m,mu,t,i,j) = lsSum[ij_dx];
                            }
                        }
                    }
                } else {
                    const scalar_type zz(0.0);
                    for(int i=0;i<Lblock;i++){
                        for(int j=0;j<Rblock;j++){
                            for(int m=0;m<Nmom;m++){
                                mat(m,mu,t,i,j) =zz;
                            }
                        }
                    }
                }
            }
        });
    } // end loop on gamma
    ////////////////////////////////////////////////////////////////////
    // This global sum is taking as much as 50% of time on 16 nodes
    // Vector size is 7 x 16 x 32 x 16 x 16 x sizeof(complex) = 2MB - 60MB depending on volume
    // Healthy size that should suffice
    ////////////////////////////////////////////////////////////////////
    if (t_gsum) *t_gsum = -usecond();
    grid->GlobalSumVector(&mat(0,0,0,0,0),Nmom*Ngamma*Nt*Lblock*Rblock);
    if (t_gsum) *t_gsum += usecond();
}

// Conserved vector current for staggered mesons
// gamma_{x,y,z} implemented for now.
// links could be thin or fat.
// KS phases in links
// no Naik term
//#include <Grid/qcd/action/fermion/ImprovedStaggeredFermion.h>
template <class FImpl>
template <typename TensorType>
void A2Autils<FImpl>::StagMesonFieldCC(TensorType &mat,
                                       //const LatticeGaugeField &Umu,
                                       const FermionField *lhs_wi,
                                       const FermionField *rhs_vj,
                                       std::vector<Gamma::Algebra> gammas,
                                       const std::vector<ComplexField > &mom,
                                       int orthogdim, double *t_kernel, double *t_gsum)
{
    typedef typename FImpl::SiteSpinor vobj;
    
    typedef typename vobj::scalar_object sobj;
    typedef typename vobj::scalar_type scalar_type;
    typedef typename vobj::vector_type vector_type;
    
    typedef iSinglet<vector_type> Singlet_v;
    typedef iSinglet<scalar_type> Singlet_s;
    
    int Lblock = mat.dimension(3);
    int Rblock = mat.dimension(4);
    
    GridBase *grid = lhs_wi[0].Grid();
    
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
    
    Vector<Singlet_v > lvSum(MFrvol);
    Vector<Singlet_s > lsSum(MFlvol);
    
    //do shift and mult outside of A2Autils
    // U_mu(x) right(x+mu)
    //FermionField Urpl(grid);
    //std::vector<LatticeColourMatrix> U(4,grid);
    //for(int mu=0;mu<Nd;mu++){
    //    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
    //}
    
    int e1=    grid->_slice_nblock[orthogdim];
    int e2=    grid->_slice_block [orthogdim];
    int stride=grid->_slice_stride[orthogdim];
    
    // potentially wasting cores here if local time extent too small
    if (t_kernel) *t_kernel = -usecond();
    
    // do x,y,z dirs

    //std::vector<FermionField> temp(Rblock, grid);
    //int mu;
    int ng=0; //doing one gamma mu (point split dir) at a time
    //for (int n = 0; n < Ngamma; n++) {

        //if ( gammas[n] == Gamma::Algebra::GammaX ) mu=0;
        //else if ( gammas[n] == Gamma::Algebra::GammaY ) mu=1;
        //else if ( gammas[n] == Gamma::Algebra::GammaZ ) mu=2;
        //else {
          //std::cout << gammas[n] << " not implemented for staggered fermion conserverd current meson field" << std::endl;
          //assert(0);
        //}

        // Re-initialize working variables before starting on a new gamma
        thread_for(r, MFrvol,{
            lvSum[r] = Zero();
        });
        thread_for(r, MFlvol,{
            lsSum[r]=scalar_type(0.0);
        });
       
 
        // do shift outside of A2AUtils.h
        //std::cout << GridLogMessage << "Cshift * Umu " << std::endl;
        //cshift must be outside thread loop
        //for(int j=0;j<Rblock;j++)
        //   temp[j] = U[mu]*Cshift(rhs_vj[j], mu, 1);
        //std::cout << GridLogMessage << "Cshift * Umu finished " << std::endl;

        thread_for(r, rd,{
            
            int so=r*grid->_ostride[orthogdim]; // base offset for start of plane
            
            for(int n=0;n<e1;n++){
                for(int b=0;b<e2;b++){
                    
                    int ss= so+n*stride+b;
                    for(int i=0;i<Lblock;i++){
                        
 		        autoView(wi_v, lhs_wi[i], CpuRead);
		        //auto wi_v  = lhs_wi[i].View();
                        auto left = conjugate(wi_v[ss]);
                        
                        for(int j=0;j<Rblock;j++){

                            //Urpl = temp[j];
                            
			    autoView(vjplU_v,rhs_vj[j],CpuRead);
			    // auto vjplU_v = rhs_vj[j].View();
                            //auto vjplU = temp[j].View();
                            auto right = vjplU_v[ss];
                            
                            Singlet_v vv;
                            vv()()() = left()()(0) * right()()(0)
                                     + left()()(1) * right()()(1)
                                     + left()()(2) * right()()(2);
                            
                            // After getting the sitewise product do the mom phase loop
                            int base = Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*r;
                            for ( int m=0;m<Nmom;m++){
                                int idx = m+base;
				autoView(mom_phase_v,mom[m],CpuRead);
                                //auto mom_phase_v = mom[m].View();
                                auto phase = mom_phase_v[ss];
                                mac(&lvSum[idx],&vv,&phase);
                            }
                        }
                    }
                }
            }
        });
        
        
        // Sum across simd lanes in the plane, breaking out orthog dir.
        thread_for(rt, rd,{
            //parallel_for(int rt=0;rt<rd;rt++){
            Coordinate icoor(Nd);
            //std::vector<int> icoor(Nd);
            //std::vector<Singlet_s> extracted(Nsimd);
            ExtractBuffer<Singlet_s> extracted(Nsimd);
            
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
        });
        if (t_kernel) *t_kernel += usecond();
        assert(mat.dimension(0) == Nmom);
        assert(mat.dimension(1) == Ngamma);
        assert(mat.dimension(2) == Nt);
        
        // ld loop and local only??
        int pd = grid->_processors[orthogdim];
        int pc = grid->_processor_coor[orthogdim];
        thread_for_collapse(2,lt,ld,{
            
            for(int pt=0;pt<pd;pt++){
                int t = lt + pt*ld;
                if (pt == pc){
                    for(int i=0;i<Lblock;i++){
                        for(int j=0;j<Rblock;j++){
                            for(int m=0;m<Nmom;m++){
                                int ij_dx = m+Nmom*i + Nmom*Lblock * j + Nmom*Lblock * Rblock * lt;
                                mat(m,ng,t,i,j) = lsSum[ij_dx];
                            }
                        }
                    }
                } else {
                    const scalar_type zz(0.0);
                    for(int i=0;i<Lblock;i++){
                        for(int j=0;j<Rblock;j++){
                            for(int m=0;m<Nmom;m++){
                                mat(m,ng,t,i,j) =zz;
                            }
                        }
                    }
                }
            }
        }); 
    //} // end loop on gamma
    ////////////////////////////////////////////////////////////////////
    // This global sum is taking as much as 50% of time on 16 nodes
    // Vector size is 7 x 16 x 32 x 16 x 16 x sizeof(complex) = 2MB - 60MB depending on volume
    // Healthy size that should suffice
    ////////////////////////////////////////////////////////////////////
    if (t_gsum) *t_gsum = -usecond();
    grid->GlobalSumVector(&mat(0,0,0,0,0),Nmom*Ngamma*Nt*Lblock*Rblock);
    if (t_gsum) *t_gsum += usecond();
}

NAMESPACE_END(Grid);

