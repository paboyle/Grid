namespace Grid { 

/*
  BlockProjector

  If _HP_BLOCK_PROJECTORS_ is defined, we assume that _evec is a basis that is not
  fully orthonormalized (to the precision of the coarse field) and we allow for higher-precision
  coarse field than basis field.

*/
//#define _HP_BLOCK_PROJECTORS_

template<typename Field>
class BlockProjector {
public:

  BasisFieldVector<Field>& _evec;
  BlockedGrid<Field>& _bgrid;

  BlockProjector(BasisFieldVector<Field>& evec, BlockedGrid<Field>& bgrid) : _evec(evec), _bgrid(bgrid) {
  }

  void createOrthonormalBasis(RealD thres = 0.0) {

    GridStopWatch sw;
    sw.Start();

    int cnt = 0;

#pragma omp parallel shared(cnt)
    {
      int lcnt = 0;

#pragma omp for
      for (int b=0;b<_bgrid._o_blocks;b++) {
	
	for (int i=0;i<_evec._Nm;i++) {
	  
	  auto nrm0 = _bgrid.block_sp(b,_evec._v[i],_evec._v[i]);
	  
	  // |i> -= <j|i> |j>
	  for (int j=0;j<i;j++) {
	    _bgrid.block_caxpy(b,_evec._v[i],-_bgrid.block_sp(b,_evec._v[j],_evec._v[i]),_evec._v[j],_evec._v[i]);
	  }
	  
	  auto nrm = _bgrid.block_sp(b,_evec._v[i],_evec._v[i]);
	  
	  auto eps = nrm/nrm0;
	  if (Reduce(eps).real() < thres) {
	    lcnt++;
	  }
	  
	  // TODO: if norm is too small, remove this eigenvector/mark as not needed; in practice: set it to zero norm here and return a mask
	  // that is then used later to decide not to write certain eigenvectors to disk (add a norm calculation before subtraction step and look at nrm/nrm0 < eps to decide)
	  _bgrid.block_cscale(b,1.0 / sqrt(nrm),_evec._v[i]);
	  
	}
	
      }

#pragma omp critical
      {
	cnt += lcnt;
      }
    }
    sw.Stop();
    std::cout << GridLogMessage << "Gram-Schmidt to create blocked basis took " << sw.Elapsed() << " (" << ((RealD)cnt / (RealD)_bgrid._o_blocks / (RealD)_evec._Nm) 
	      << " below threshold)" << std::endl;

  }

  template<typename CoarseField>
  void coarseToFine(const CoarseField& in, Field& out) {

    out = zero;
    out.checkerboard = _evec._v[0].checkerboard;

    int Nbasis = sizeof(in._odata[0]._internal._internal) / sizeof(in._odata[0]._internal._internal[0]);
    assert(Nbasis == _evec._Nm);
    
#pragma omp parallel for
    for (int b=0;b<_bgrid._o_blocks;b++) {
      for (int j=0;j<_evec._Nm;j++) {
	_bgrid.block_caxpy(b,out,in._odata[b]._internal._internal[j],_evec._v[j],out);
      }
    }

  }

  template<typename CoarseField>
  void fineToCoarse(const Field& in, CoarseField& out) {

    out = zero;

    int Nbasis = sizeof(out._odata[0]._internal._internal) / sizeof(out._odata[0]._internal._internal[0]);
    assert(Nbasis == _evec._Nm);


    Field tmp(_bgrid._grid);
    tmp = in;
    
#pragma omp parallel for
    for (int b=0;b<_bgrid._o_blocks;b++) {
      for (int j=0;j<_evec._Nm;j++) {
	// |rhs> -= <j|rhs> |j>
	auto c = _bgrid.block_sp(b,_evec._v[j],tmp);
	_bgrid.block_caxpy(b,tmp,-c,_evec._v[j],tmp); // may make this more numerically stable
	out._odata[b]._internal._internal[j] = c;
      }
    }

  }

  template<typename CoarseField>
    void deflateFine(BasisFieldVector<CoarseField>& _coef,const std::vector<RealD>& eval,int N,const Field& src_orig,Field& result) {
    result = zero;
    for (int i=0;i<N;i++) {
      Field tmp(result._grid);
      coarseToFine(_coef._v[i],tmp);
      axpy(result,TensorRemove(innerProduct(tmp,src_orig)) / eval[i],tmp,result);
    }
  }

  template<typename CoarseField>
    void deflateCoarse(BasisFieldVector<CoarseField>& _coef,const std::vector<RealD>& eval,int N,const Field& src_orig,Field& result) {
    CoarseField src_coarse(_coef._v[0]._grid);
    CoarseField result_coarse = src_coarse;
    result_coarse = zero;
    fineToCoarse(src_orig,src_coarse);
    for (int i=0;i<N;i++) {
      axpy(result_coarse,TensorRemove(innerProduct(_coef._v[i],src_coarse)) / eval[i],_coef._v[i],result_coarse);
    }
    coarseToFine(result_coarse,result);
  }

  template<typename CoarseField>
    void deflate(BasisFieldVector<CoarseField>& _coef,const std::vector<RealD>& eval,int N,const Field& src_orig,Field& result) {
    // Deflation on coarse Grid is much faster, so use it by default.  Deflation on fine Grid is kept for legacy reasons for now.
    deflateCoarse(_coef,eval,N,src_orig,result);
  }

};
}
