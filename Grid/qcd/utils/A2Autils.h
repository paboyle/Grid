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

  static void PionFieldWVmom(Eigen::Tensor<ComplexD,4> &mat, 
			     const FermionField *wi,
			     const FermionField *vj,
			     const std::vector<ComplexField > &mom,
			     int orthogdim);


  static void NucleonFieldMom(Eigen::Tensor<ComplexD,6> &mat, 
                              const std::vector<FermionField> &one,
                              const std::vector<FermionField> &two,
                              const std::vector<FermionField> &three,
				     const std::vector<ComplexField > &mom,
				     int parity,
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

  template <typename TensorType> // output: rank 5 tensor, e.g. Eigen::Tensor<ComplexD, 5>
  static void AslashField(TensorType &mat, 
        const FermionField *lhs_wi,
        const FermionField *rhs_vj,
        const std::vector<ComplexField> &emB0,
        const std::vector<ComplexField> &emB1,
        int orthogdim, double *t_kernel = nullptr, double *t_gsum = nullptr);

  static void ContractWWVV(std::vector<PropagatorField> &WWVV,
			   const Eigen::Tensor<ComplexD,3> &WW_sd,
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
};

template<class FImpl>
void A2Autils<FImpl>::NucleonFieldMom(Eigen::Tensor<ComplexD,6> &mat, 
				     const std::vector<FermionField> &one,
				     const std::vector<FermionField> &two,
				     const std::vector<FermionField> &three,
				     const std::vector<ComplexField > &mom,
				     int parity,
				     int orthogdim) 
{
  assert(parity == 1 || parity == -1);

  typedef typename FImpl::SiteSpinor vobj;

  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  typedef iSpinVector<vector_type> SpinVector_v;
  typedef iSpinVector<scalar_type> SpinVector_s;
  
  int oneBlock = mat.dimension(2); 
  int twoBlock = mat.dimension(3);
  int threeBlock = mat.dimension(4);

  GridBase *grid = one[0].Grid();

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
  int MFrvol = rd*oneBlock*twoBlock*threeBlock*Nmom;
  int MFlvol = ld*oneBlock*twoBlock*threeBlock*Nmom;

  Vector<SpinVector_v > lvSum(MFrvol);
  accelerator_for (r, MFrvol, Nsimd, {
    lvSum[r] = 0;
  } );

  Vector<SpinVector_s > lsSum(MFlvol);             
  accelerator_for (r, MFlvol, Nsimd, {
    lsSum[r] = scalar_type(0.0);
  } );

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];

  accelerator_for(r, rd, Nsimd, {
    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane

    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){

	int ss= so+n*stride+b;

	for(int i=0;i<oneBlock;i++){

	  auto v1 = one[i]._odata[ss];
	  auto gv1 = Gamma(Gamma::Algebra::GammaT)*v1;
	  auto pv1 = 0.5*(v1 + (double)parity*gv1);

       	  for(int j=0;j<twoBlock;j++){

	    //auto v2 = conjugate(two[j]._odata[ss]);
	    auto v2 = two[j]._odata[ss];
            // C = i gamma_2 gamma_4 => C gamma_5 = - i gamma_1 gamma_3
	    //auto v2g = v2*Gamma(Gamma::Algebra::SigmaXZ);
            //auto v2g=v2;  

	    for(int k=0;k<threeBlock;k++){

	      auto v3 = three[k]._odata[ss];
	      auto gv3 = Gamma(Gamma::Algebra::SigmaXZ)*v3;
	      SpinVector_v vv;

              for(int s1=0;s1<Ns;s1++){
              vv()(s1)() = 0;
	      for(int s2=0;s2<Ns;s2++){            
	     /*   vv()(s1)() =  pv1()(s1)(0) * v2g()(s2)(1) * v3()(s2)(2)   //Cross product
                  -           pv1()(s1)(0) * v2g()(s2)(2) * v3()(s2)(1)    
                  +           pv1()(s1)(1) * v2g()(s2)(2) * v3()(s2)(0)    
                  -           pv1()(s1)(1) * v2g()(s2)(0) * v3()(s2)(2)    
                  +           pv1()(s1)(2) * v2g()(s2)(0) * v3()(s2)(1)    
                  -           pv1()(s1)(2) * v2g()(s2)(1) * v3()(s2)(0);    */
	        vv()(s1)() +=  pv1()(s1)(0) * v2()(s2)(1) * gv3()(s2)(2)   //Cross product
                  -           pv1()(s1)(0) * v2()(s2)(2) * gv3()(s2)(1)    
                  +           pv1()(s1)(1) * v2()(s2)(2) * gv3()(s2)(0)    
                  -           pv1()(s1)(1) * v2()(s2)(0) * gv3()(s2)(2)    
                  +           pv1()(s1)(2) * v2()(s2)(0) * gv3()(s2)(1)    
                  -           pv1()(s1)(2) * v2()(s2)(1) * gv3()(s2)(0);    
              }}

              /*if (i+j+k == 0) {
                Serializable::WriteMember(std::cout, pv1);
                Serializable::WriteMember(std::cout, v2);
                Serializable::WriteMember(std::cout, gv3);
                Serializable::WriteMember(std::cout, vv);
              }*/
	    
	      // After getting the sitewise product do the mom phase loop
              int base = Nmom*i+Nmom*oneBlock*j+Nmom*oneBlock*twoBlock*k+Nmom*oneBlock*twoBlock*threeBlock*r;
	      for ( int m=0;m<Nmom;m++){
	        int idx = m+base;
	        auto phase = mom[m]._odata[ss];
                for(int is=0;is<Ns;is++){
	          mac(&lvSum[idx]()(is)(),&vv()(is)(),&phase()()());
		}
              }
	    }
	  }
	}
      }
    }
  } );


  // Sum across simd lanes in the plane, breaking out orthog dir.
  accelerator_for(rt, rd, Nsimd, {

    Coordinate icoor(nd);
    std::vector<SpinVector_s>  extracted(Nsimd);               

    for(int i=0;i<oneBlock;i++){
    for(int j=0;j<twoBlock;j++){
    for(int k=0;k<threeBlock;k++){
    for(int m=0;m<Nmom;m++){

      int ij_rdx = m+Nmom*i + Nmom*oneBlock * j + Nmom*oneBlock * twoBlock * k + Nmom*oneBlock * twoBlock *threeBlock * rt;

      extract(lvSum[ij_rdx],extracted);

      for(int idx=0;idx<Nsimd;idx++){

	grid->iCoorFromIindex(icoor,idx);

	int ldx    = rt+icoor[orthogdim]*rd;

	int ij_ldx = m+Nmom*i + Nmom*oneBlock * j + Nmom*oneBlock * twoBlock * k + Nmom*oneBlock * twoBlock *threeBlock * ldx;

	lsSum[ij_ldx]=lsSum[ij_ldx]+extracted[idx];

      }
    }}}}
  } );

  assert(mat.dimension(0) == Nmom);
  assert(mat.dimension(1) == Nt);

  int pd = grid->_processors[orthogdim];
  int pc = grid->_processor_coor[orthogdim];
  accelerator_for(lt, ld, Nsimd,
  {
    for(int pt=0;pt<pd;pt++){
      int t = lt + pt*ld;
      if (pt == pc){
        for(int i=0;i<oneBlock;i++){
          for(int j=0;j<twoBlock;j++){
            for(int k=0;k<threeBlock;k++){
	      for(int m=0;m<Nmom;m++){
	        int ij_dx = m+Nmom*i + Nmom*oneBlock * j + Nmom*oneBlock * twoBlock * k + Nmom*oneBlock * twoBlock *threeBlock * lt;
	        for(int is=0;is<4;is++){
	          mat(m,t,i,j,k,is) = lsSum[ij_dx]()(is)();
                }
              }
	    }
	  }
	}
      } else { 
	const scalar_type zz(0.0);
        for(int i=0;i<oneBlock;i++){
          for(int j=0;j<twoBlock;j++){
            for(int k=0;k<threeBlock;k++){
	      for(int m=0;m<Nmom;m++){
	        for(int is=0;is<4;is++){
	          mat(m,t,i,j,k,is) =zz;
                }
              }
	    }
	  }
	}
      }
    }
  } );

  grid->GlobalSumVector(&mat(0,0,0,0,0,0),Nmom*Nt*oneBlock*twoBlock*threeBlock*4);
}





/*
template <class FImpl>
template <typename TensorType>
void A2Autils<FImpl>::BaryonField(TensorType &mat, 
				 const FermionField *one,
				 const FermionField *two,
				 const FermionField *three,
				 std::vector<Gamma::Algebra> gammaA,
				 std::vector<Gamma::Algebra> gammaB,
				 const std::vector<ComplexField > &mom,
				 int orthogdim, double *t_kernel, double *t_gsum) 
{
  typedef typename FImpl::SiteSpinor vobj;

  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iSpinMatrix<scalar_type> SpinMatrix_s;
  
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;

  int oneBlock = mat.dimension(3); 
  int twoBlock = mat.dimension(4);
  int threeBlock = mat.dimension(5);

  GridBase *grid = one[0].Grid();
  
  const int    Nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();

  int Nt     = grid->GlobalDimensions()[orthogdim];
  int Ngamma = gammaA.size();
  assert (Ngamma = gammaB.size()); // Only combinatin of two gammas gives correct operator!
  int Nmom   = mom.size();

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  // will locally sum vectors first
  // sum across these down to scalars
  // splitting the SIMD
  int MFrvol = rd*twoBlock*threeBlock*Nmom;
  int MFlvol = ld*twoBlock*threeBlock*Nmom;

  Vector<Vector<SpinMatrix_v >> lvSum(3);
  for (int ic=0;ic<3;ic++){
    lvSum[ic].resize(MFrvol);
    parallel_for (int r = 0; r < MFrvol; r++){
      lvSum[ic][r] = zero;
    }
  }

  Vector<Vector<SpinMatrix_s >> lsSum(3);
  for (int ic=0;ic<3;ic++){
    lsSum[ic].resize(MFlvol);
    parallel_for (int r = 0; r < MFlvol; r++){
      lsSum[ic][r] = scalar_type(0.0);
    }
  }

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];

  // potentially wasting cores here if local time extent too small
  if (t_kernel) *t_kernel = -usecond();
  parallel_for(int r=0;r<rd;r++){

    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

    // first, the diquark two*gammaB*three
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){

	int ss= so+n*stride+b;

	  for(int j=0;j<twoBlock;j++){

	    auto two_j = two[j]._odata[ss];

	  for(int k=0;k<threeBlock;k++){

	    auto three_k = three[j]._odata[ss];

	    Vector<SpinMatrix_v > vv(3);

	    for(int s1=0;s1<Ns;s1++){
	    for(int s2=0;s2<Ns;s2++){
	      vv[0]()(s1,s2)() = two_j()(s2)(1) * three_k()(s1)(2)           //ideal would be SpinMatrix but ColourVector...
                -                two_j()(s2)(2) * three_k()(s1)(1);          //this is the cross product (two x three)^i
	      vv[1]()(s1,s2)() = two_j()(s2)(2) * three_k()(s1)(0)           
                -                two_j()(s2)(0) * three_k()(s1)(2);         
	      vv[2]()(s1,s2)() = two_j()(s2)(0) * three_k()(s1)(1)           
                -                two_j()(s2)(1) * three_k()(s1)(0);         
	    }}
	    
	    // After getting the sitewise product do the mom phase loop
	    int base = Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*r;
	    for ( int m=0;m<Nmom;m++){
	      for ( int ic=0;ic<3;ic++){
	        int idx = m+base;
	        auto phase = mom[m]._odata[ss];
	        mac(&lvSum[ic][idx],&vv,&phase);
	      }
	    }
	 
	  }
	}
      }
    }
  }


  for ( int ic=0;ic<3;ic++){
  // Sum across simd lanes in the plane, breaking out orthog dir.
  parallel_for(int rt=0;rt<rd;rt++){

    std::vector<int> icoor(Nd);
    std::vector<SpinMatrix_s> extracted(Nsimd);               

    for(int i=0;i<twoBlock;i++){
    for(int j=0;j<threeBlock;j++){
    for(int m=0;m<Nmom;m++){

      int ij_rdx = m+Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*rt;

      extract(lvSum[ic][ij_rdx],extracted);

      for(int idx=0;idx<Nsimd;idx++){

	grid->iCoorFromIindex(icoor,idx);

	int ldx    = rt+icoor[orthogdim]*rd;

	int ij_ldx = m+Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*ldx;

	lsSum[ic][ij_ldx]=lsSum[ic][ij_ldx]+extracted[idx];

      }
    }}}
  }
 

  if (t_kernel) *t_kernel += usecond();
  assert(mat.dimension(0) == Nmom);
  assert(mat.dimension(1) == Ngamma);
  assert(mat.dimension(2) == Nt);

  TensorType diquark; // Need this instead of mat!!!

  // ld loop and local only??
  int pd = grid->_processors[orthogdim];
  int pc = grid->_processor_coor[orthogdim];
  parallel_for_nest2(int lt=0;lt<ld;lt++)
  {
    for(int pt=0;pt<pd;pt++){
      int t = lt + pt*ld;
      if (pt == pc){
	for(int i=0;i<twoBlock;i++){
	  for(int j=0;j<threeBlock;j++){
	    for(int m=0;m<Nmom;m++){
	      int ij_dx = m+Nmom*i + Nmom*Lblock * j + Nmom*Lblock * Rblock * lt;
	      for(int mu=0;mu<Ngamma;mu++){
		// this is a bit slow
		mat(m,mu,t,i,j) = trace(lsSum[ic][ij_dx]*Gamma(gammaB[mu]));
	      }
	    }
	  }
	}
      } else { 
	const scalar_type zz(0.0);
	for(int i=0;i<twoBlock;i++){
	  for(int j=0;j<threeBlock;j++){
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
  }
  ////////////////////////////////////////////////////////////////////
  // This global sum is taking as much as 50% of time on 16 nodes
  // Vector size is 7 x 16 x 32 x 16 x 16 x sizeof(complex) = 2MB - 60MB depending on volume
  // Healthy size that should suffice
  ////////////////////////////////////////////////////////////////////
  if (t_gsum) *t_gsum = -usecond();
  grid->GlobalSumVector(&mat(0,0,0,0,0),Nmom*Ngamma*Nt*Lblock*Rblock);
  if (t_gsum) *t_gsum += usecond();
}
*/

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

	  auto lhs_v = lhs_wi[i].View();
	  auto left = conjugate(lhs_v[ss]);

	  for(int j=0;j<Rblock;j++){

	    SpinMatrix_v vv;
	    auto rhs_v = rhs_vj[j].View();
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
	      auto mom_v = mom[m].View();
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

	  auto wi_v = wi[i].View();
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
	    
	    auto vj_v=vj[j].View();
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

	  auto wi_v = wi[i].View();
	  auto w = conjugate(wi_v[ss]);

	  for(int j=0;j<Rblock;j++){
	    
	    auto vj_v = vj[j].View();
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
	      auto mom_v = mom[m].View();
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
  	        auto wi_v = lhs_wi[i].View();
                auto left = conjugate(wi_v[ss]);

                for(int j=0;j<Rblock;j++)
                {
                    SpinMatrix_v vv;
		    auto vj_v  = rhs_vj[j].View();
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
  		        auto emB0_v = emB0[m].View();
  		        auto emB1_v = emB1[m].View();
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

template<class FImpl>
void A2Autils<FImpl>::ContractWWVV(std::vector<PropagatorField> &WWVV,
				   const Eigen::Tensor<ComplexD,3> &WW_sd,
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
	auto vs_v = vs[s].View();
	auto tmp1 = vs_v[ss];
	vobj tmp2 = Zero();
	vobj tmp3 = Zero();
	for(int d=d_o;d<MIN(d_o+d_unroll,N_d);d++){
	  auto vd_v = vd[d].View();
	  Scalar_v coeff = WW_sd(t,s,d);
	  tmp3 = conjugate(vd_v[ss]);
	  mac(&tmp2, &coeff, &tmp3);
	}

	//////////////////////////
	// Fast outer product of tmp1 with a sum of terms suppressed by d_unroll
	//////////////////////////
	auto WWVV_v = WWVV[t].View();
	for(int s1=0;s1<Ns;s1++){
	for(int s2=0;s2<Ns;s2++){
	  WWVV_v[ss]()(s1,s2)(0,0) += tmp1()(s1)(0)*tmp2()(s2)(0);
	  WWVV_v[ss]()(s1,s2)(0,1) += tmp1()(s1)(0)*tmp2()(s2)(1);
	  WWVV_v[ss]()(s1,s2)(0,2) += tmp1()(s1)(0)*tmp2()(s2)(2);
	  WWVV_v[ss]()(s1,s2)(1,0) += tmp1()(s1)(1)*tmp2()(s2)(0);
	  WWVV_v[ss]()(s1,s2)(1,1) += tmp1()(s1)(1)*tmp2()(s2)(1);
	  WWVV_v[ss]()(s1,s2)(1,2) += tmp1()(s1)(1)*tmp2()(s2)(2);
	  WWVV_v[ss]()(s1,s2)(2,0) += tmp1()(s1)(2)*tmp2()(s2)(0);
	  WWVV_v[ss]()(s1,s2)(2,1) += tmp1()(s1)(2)*tmp2()(s2)(1);
	  WWVV_v[ss]()(s1,s2)(2,2) += tmp1()(s1)(2)*tmp2()(s2)(2);
	}}

      }}
    }
  });
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

  auto WWVV0_v = WWVV0.View();
  auto WWVV1_v = WWVV1.View();
  auto O_trtr_v= O_trtr.View();
  auto O_fig8_v= O_fig8.View();
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

  auto WWVV0_v = WWVV0.View();
  auto WWVV1_v = WWVV1.View();
  auto O_trtr_v= O_trtr.View();
  auto O_fig8_v= O_fig8.View();

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

NAMESPACE_END(Grid);

