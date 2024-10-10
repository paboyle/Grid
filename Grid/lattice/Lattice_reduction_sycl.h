NAMESPACE_BEGIN(Grid);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Possibly promote to double and sum
/////////////////////////////////////////////////////////////////////////////////////////////////////////


template <class vobj>
inline typename vobj::scalar_objectD sumD_gpu_tensor(const vobj *lat, Integer osites) 
{
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_objectD sobjD;

  sobj identity; zeroit(identity);
  sobj ret; zeroit(ret);
  Integer nsimd= vobj::Nsimd();
  { 
    sycl::buffer<sobj, 1> abuff(&ret, {1});
    theGridAccelerator->submit([&](sycl::handler &cgh) {
      auto Reduction = sycl::reduction(abuff,cgh,identity,std::plus<>());
      cgh.parallel_for(sycl::range<1>{osites},
                      Reduction,
                      [=] (sycl::id<1> item, auto &sum) {
                        auto osite   = item[0];
                        sum +=Reduce(lat[osite]);
                      });
    });
  }
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
  Word identity;  identity=0;
  Word ret = 0;
  { 
    sycl::buffer<Word, 1> abuff(&ret, {1});
    theGridAccelerator->submit([&](sycl::handler &cgh) {
      auto Reduction = sycl::reduction(abuff,cgh,identity,std::bit_xor<>());
      cgh.parallel_for(sycl::range<1>{L},
                      Reduction,
                      [=] (sycl::id<1> index, auto &sum) {
                        sum ^=vec[index];
                      });
    });
  }
  theGridAccelerator->wait();
  return ret;
}

NAMESPACE_END(Grid);

