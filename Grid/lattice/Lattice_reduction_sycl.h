NAMESPACE_BEGIN(Grid);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Possibly promote to double and sum
/////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class vobj>
inline typename vobj::scalar_objectD sumD_gpu_tensor(const vobj *lat, Integer osites) 
{
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_objectD sobjD;
  sobj *mysum =(sobj *) malloc_shared(sizeof(sobj),*theGridAccelerator);
  sobj identity; zeroit(identity);
  sobj ret ; 

  Integer nsimd= vobj::Nsimd();
  
  theGridAccelerator->submit([&](cl::sycl::handler &cgh) {
     auto Reduction = cl::sycl::reduction(mysum,identity,std::plus<>());
     cgh.parallel_for(cl::sycl::range<1>{osites},
		      Reduction,
		      [=] (cl::sycl::id<1> item, auto &sum) {
      auto osite   = item[0];
      sum +=Reduce(lat[osite]);
     });
   });
  theGridAccelerator->wait();
  ret = mysum[0];
  free(mysum,*theGridAccelerator);
  sobjD dret; convertType(dret,ret);
  return dret;
}

template <class vobj>
inline typename vobj::scalar_objectD sumD_gpu_large(const vobj *lat, Integer osites)
{
  return sumD_gpu_tensor(lat,osites);
}
template <class vobj>
inline typename vobj::scalar_objectD sumD_gpu_small(const vobj *lat, Integer osites)
{
  return sumD_gpu_large(lat,osites);
}

template <class vobj>
inline typename vobj::scalar_objectD sumD_gpu(const vobj *lat, Integer osites)
{
  return sumD_gpu_large(lat,osites);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Return as same precision as input performing reduction in double precision though
/////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class vobj>
inline typename vobj::scalar_object sum_gpu(const vobj *lat, Integer osites) 
{
  typedef typename vobj::scalar_object sobj;
  sobj result;
  result = sumD_gpu(lat,osites);
  return result;
}

template <class vobj>
inline typename vobj::scalar_object sum_gpu_large(const vobj *lat, Integer osites)
{
  typedef typename vobj::scalar_object sobj;
  sobj result;
  result = sumD_gpu_large(lat,osites);
  return result;
}

NAMESPACE_END(Grid);

/*
template<class Double> Double svm_reduce(Double *vec,uint64_t L)
{
  Double sumResult; zeroit(sumResult);
  Double *d_sum =(Double *)cl::sycl::malloc_shared(sizeof(Double),*theGridAccelerator);
  Double identity;  zeroit(identity);
  theGridAccelerator->submit([&](cl::sycl::handler &cgh) {
     auto Reduction = cl::sycl::reduction(d_sum,identity,std::plus<>());
     cgh.parallel_for(cl::sycl::range<1>{L},
		      Reduction,
		      [=] (cl::sycl::id<1> index, auto &sum) {
	 sum +=vec[index];
     });
   });
  theGridAccelerator->wait();
  Double ret = d_sum[0];
  free(d_sum,*theGridAccelerator);
  std::cout << " svm_reduce finished "<<L<<" sites sum = " << ret <<std::endl;
  return ret;
}

template <class vobj>
inline typename vobj::scalar_objectD sumD_gpu_repack(const vobj *lat, Integer osites)
{
  typedef typename vobj::vector_type  vector;
  typedef typename vobj::scalar_type  scalar;

  typedef typename vobj::scalar_typeD scalarD;
  typedef typename vobj::scalar_objectD sobjD;

  sobjD ret;
  scalarD *ret_p = (scalarD *)&ret;
  
  const int nsimd = vobj::Nsimd();
  const int words = sizeof(vobj)/sizeof(vector);

  Vector<scalar> buffer(osites*nsimd);
  scalar *buf = &buffer[0];
  vector *dat = (vector *)lat;

  for(int w=0;w<words;w++) {

    accelerator_for(ss,osites,nsimd,{
	int lane = acceleratorSIMTlane(nsimd);
	buf[ss*nsimd+lane] = dat[ss*words+w].getlane(lane);
    });
    //Precision change at this point is to late to gain precision
    ret_p[w] = svm_reduce(buf,nsimd*osites);
  }
  return ret;
}
*/
