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

#ifdef __linux
    if (value->IsBoss()) {
      system("cat /proc/meminfo");
    }
#endif

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
    basisOrthogonalize(_v,w,k);
  }

  void rotate(Eigen::MatrixXd& Qt,int j0, int j1, int k0,int k1,int Nm) {
    basisRotate(_v,Qt,j0,j1,k0,k1,Nm);
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

  void sortInPlace(std::vector<RealD>& sort_vals, bool reverse) {
    basisSortInPlace(_v,sort_vals,reverse);
  }

  void deflate(const std::vector<RealD>& eval,const Field& src_orig,Field& result) {
    basisDeflate(_v,eval,src_orig,result);
  }

 }; 
}
