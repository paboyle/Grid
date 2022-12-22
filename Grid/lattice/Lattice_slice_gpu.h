NAMESPACE_BEGIN(Grid);

// If NOT CUDA or HIP -- we should provide
// -- atomicAdd(float *,float)
// -- atomicAdd(double *,double)
// 
// Augment CUDA with complex atomics
#if !defined(GRID_HIP) || !defined(GRID_CUDA)
inline void atomicAdd(float *acc,float elem)
{
  *acc += elem;
}
inline void atomicAdd(double *acc,double elem)
{
  *acc += elem;
}
#endif
inline void atomicAdd(ComplexD *accum,ComplexD & elem)
{
  double *a_p = (double *)accum;
  double *e_p = (double *)&elem;
  for(int w=0;w<2;w++){
    atomicAdd(&a_p[w],e_p[w]);
  }
}
inline void atomicAdd(ComplexF *accum,ComplexF & elem)
{
  float *a_p = (float *)accum;
  float *e_p = (float *)&elem;
  for(int w=0;w<2;w++){
    atomicAdd(&a_p[w],e_p[w]);
  }
}
// Augment CUDA with vobj atomics
template<class vobj> accelerator_inline void atomicAdd(vobj *accum, vobj & elem)
{
  typedef typename vobj::scalar_type scalar_type;
  scalar_type *a_p= (scalar_type *)accum;
  scalar_type *e_p= (scalar_type *)& elem;
  for(int w=0;w<vobj::Nsimd();w++){
    atomicAdd(&a_p[w],e_p[w]);
  }
}
// Atomics based slice sum
template<class vobj> inline void sliceSumGpu(const Lattice<vobj> &Data,std::vector<typename vobj::scalar_object> &result,int orthogdim)
{
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_object::scalar_type scalar_type;
  GridBase  *grid = Data.Grid();
  assert(grid!=NULL);

  const int    Nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();

  assert(orthogdim >= 0);
  assert(orthogdim < Nd);

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  // Move to device memory and copy in / out
  Vector<vobj> lvSum(rd); // will locally sum vectors first
  Vector<sobj> lsSum(ld,Zero());                    // sum across these down to scalars
  ExtractBuffer<sobj> extracted(Nsimd);                  // splitting the SIMD

  result.resize(fd); // And then global sum to return the same vector to every node 
  for(int r=0;r<rd;r++){
    lvSum[r]=Zero();
  }

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];

  // sum over reduced dimension planes, breaking out orthog dir
  // Parallel over orthog direction
  autoView( Data_v, Data, AcceleratorRead);
  auto lvSum_p=&lvSum[0];
  int ostride = grid->_ostride[orthogdim]; 
  accelerator_for( ree,rd*e1*e2,1, {
    int b = ree%e2;
    int re= ree/e2;
    int n=re%e1;
    int r=re/e1;
    int so=r*ostride;
    int ss=so+n*stride+b;
    atomicAdd(&lvSum_p[r],Data_v[ss]);
  });

  // Sum across simd lanes in the plane, breaking out orthog dir.
  Coordinate icoor(Nd);

  for(int rt=0;rt<rd;rt++){

    extract(lvSum[rt],extracted);

    for(int idx=0;idx<Nsimd;idx++){

      grid->iCoorFromIindex(icoor,idx);

      int ldx =rt+icoor[orthogdim]*rd;

      lsSum[ldx]=lsSum[ldx]+extracted[idx];

    }
  }
  
  // sum over nodes.
  for(int t=0;t<fd;t++){
    int pt = t/ld; // processor plane
    int lt = t%ld;
    if ( pt == grid->_processor_coor[orthogdim] ) {
      result[t]=lsSum[lt];
    } else {
      result[t]=Zero();
    }

  }
  scalar_type * ptr = (scalar_type *) &result[0];
  int words = fd*sizeof(sobj)/sizeof(scalar_type);
  grid->GlobalSumVector(ptr, words);
}


NAMESPACE_END(Grid);
