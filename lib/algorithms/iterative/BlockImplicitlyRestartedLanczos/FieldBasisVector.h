namespace Grid { 

template<class Field>
class BasisFieldVector {
 public:
  int _Nm;

  typedef typename Field::scalar_type Coeff_t;
  typedef typename Field::vector_type vCoeff_t;
  typedef typename Field::vector_object vobj;
  typedef typename vobj::scalar_object sobj;

  std::vector<Field> _v; // _Nfull vectors

  void report(int n,GridBase* value) {

    std::cout << GridLogMessage << "BasisFieldVector allocated:\n";
    std::cout << GridLogMessage << " Delta N = " << n << "\n";
    std::cout << GridLogMessage << " Size of full vectors (size) = " << 
      ((double)n*sizeof(vobj)*value->oSites() / 1024./1024./1024.) << " GB\n";
    std::cout << GridLogMessage << " Size = " << _v.size() << " Capacity = " << _v.capacity() << std::endl;

    value->Barrier();

    if (value->IsBoss()) {
      system("cat /proc/meminfo");
    }

    value->Barrier();

  }

  BasisFieldVector(int Nm,GridBase* value) : _Nm(Nm), _v(Nm,value) {
    report(Nm,value);
  }
  
  ~BasisFieldVector() {
  }

  Field& operator[](int i) {
    return _v[i];
  }

  void orthogonalize(Field& w, int k) {
    for(int j=0; j<k; ++j){
      Coeff_t ip = (Coeff_t)innerProduct(_v[j],w);
      w = w - ip*_v[j];
    }
  }

  void rotate(std::vector<RealD>& Qt,int j0, int j1, int k0,int k1,int Nm) {
    
    GridBase* grid = _v[0]._grid;
      
#pragma omp parallel
    {
      std::vector < vobj > B(Nm);
      
#pragma omp for
      for(int ss=0;ss < grid->oSites();ss++){
	for(int j=j0; j<j1; ++j) B[j]=0.;
	
	for(int j=j0; j<j1; ++j){
	  for(int k=k0; k<k1; ++k){
	    B[j] +=Qt[k+Nm*j] * _v[k]._odata[ss];
	  }
	}
	for(int j=j0; j<j1; ++j){
	  _v[j]._odata[ss] = B[j];
	}
      }
    }

  }

  size_t size() const {
    return _Nm;
  }

  void resize(int n) {
    if (n > _Nm)
      _v.reserve(n);
    
    _v.resize(n,_v[0]._grid);

    if (n < _Nm)
      _v.shrink_to_fit();

    report(n - _Nm,_v[0]._grid);

    _Nm = n;
  }

  std::vector<int> getIndex(std::vector<RealD>& sort_vals) {

    std::vector<int> idx(sort_vals.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
	 [&sort_vals](int i1, int i2) {return ::fabs(sort_vals[i1]) < ::fabs(sort_vals[i2]);});

    return idx;
  }

  void reorderInPlace(std::vector<RealD>& sort_vals, std::vector<int>& idx) {
    GridStopWatch gsw;
    gsw.Start();

    int nswaps = 0;
    for (size_t i=0;i<idx.size();i++) {
      if (idx[i] != i) {

	// find proper place (this could be done in logarithmic time, don't bother for now)
	size_t j;
	for (j=i;j<idx.size();j++)
	  if (idx[j]==i)
	    break;
	assert(j!=idx.size());
	
	Field _t(_v[0]._grid);
	_t = _v[idx[j]];
	_v[idx[j]] = _v[idx[i]];
	_v[idx[i]] = _t;

	RealD _td = sort_vals[idx[j]];
	sort_vals[idx[j]] = sort_vals[idx[i]];
	sort_vals[idx[i]] = _td;

	int _tt = idx[i];
	idx[i] = idx[j];
	idx[j] = _tt;
	
	nswaps++;
      }
    }

    // sort values
    gsw.Stop();
    std::cout << GridLogMessage << "Sorted eigenspace in place in " << gsw.Elapsed() << " using " << nswaps << " swaps" << std::endl;
  }

  void sortInPlace(std::vector<RealD>& sort_vals, bool reverse) {

    std::vector<int> idx = getIndex(sort_vals);
    if (reverse)
      std::reverse(idx.begin(), idx.end());

    reorderInPlace(sort_vals,idx);

  }

  void deflate(const std::vector<RealD>& eval,const Field& src_orig,Field& result) {
    result = zero;
    int N = (int)_v.size();
    for (int i=0;i<N;i++) {
      Field& tmp = _v[i];
      axpy(result,TensorRemove(innerProduct(tmp,src_orig)) / eval[i],tmp,result);
    }
  }

 }; 
}
