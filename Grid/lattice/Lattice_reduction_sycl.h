NAMESPACE_BEGIN(Grid);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Possibly promote to double and sum
/////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class vobj>
inline typename vobj::scalar_objectD sumD_gpu_tensor(const vobj *lat, Integer osites) 
{
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_objectD sobjD;
  static Vector<sobj> mysum;
  mysum.resize(1);
  sobj *mysum_p = & mysum[0];
  sobj identity; zeroit(identity);
  mysum[0] = identity;
  sobj ret ; 

  Integer nsimd= vobj::Nsimd();

  const cl::sycl::property_list PropList ({ cl::sycl::property::reduction::initialize_to_identity() });
  theGridAccelerator->submit([&](cl::sycl::handler &cgh) {
    auto Reduction = cl::sycl::reduction(mysum_p,identity,std::plus<>(),PropList);
     cgh.parallel_for(cl::sycl::range<1>{osites},
		      Reduction,
		      [=] (cl::sycl::id<1> item, auto &sum) {
      auto osite   = item[0];
      sum +=Reduce(lat[osite]);
     });
   });
  theGridAccelerator->wait();
  ret = mysum[0];
  //  free(mysum,*theGridAccelerator);
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


template<class Word> Word svm_xor(Word *vec,uint64_t L)
{
  Word xorResult; xorResult = 0;
  static Vector<Word> d_sum;
  d_sum.resize(1);
  Word *d_sum_p=&d_sum[0];
  Word identity;  identity=0;
  d_sum[0] = identity;
  const cl::sycl::property_list PropList ({ cl::sycl::property::reduction::initialize_to_identity() });
  theGridAccelerator->submit([&](cl::sycl::handler &cgh) {
    auto Reduction = cl::sycl::reduction(d_sum_p,identity,std::bit_xor<>(),PropList);
     cgh.parallel_for(cl::sycl::range<1>{L},
		      Reduction,
		      [=] (cl::sycl::id<1> index, auto &sum) {
	 sum^=vec[index];
     });
   });
  theGridAccelerator->wait();
  Word ret = d_sum[0];
  //  free(d_sum,*theGridAccelerator);
  return ret;
}

NAMESPACE_END(Grid);

/*

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
