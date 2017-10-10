namespace Grid {

template<typename Field>
class BlockedGrid {
public:
  GridBase* _grid;
  typedef typename Field::scalar_type  Coeff_t;
  typedef typename Field::vector_type vCoeff_t;
  
  std::vector<int> _bs; // block size
  std::vector<int> _nb; // number of blocks
  std::vector<int> _l;  // local dimensions irrespective of cb
  std::vector<int> _l_cb;  // local dimensions of checkerboarded vector
  std::vector<int> _l_cb_o;  // local dimensions of inner checkerboarded vector
  std::vector<int> _bs_cb; // block size in checkerboarded vector
  std::vector<int> _nb_o; // number of blocks of simd o-sites

  int _nd, _blocks, _cf_size, _cf_block_size, _cf_o_block_size, _o_blocks, _block_sites;
  
  BlockedGrid(GridBase* grid, const std::vector<int>& block_size) :
    _grid(grid), _bs(block_size), _nd((int)_bs.size()), 
      _nb(block_size), _l(block_size), _l_cb(block_size), _nb_o(block_size),
      _l_cb_o(block_size), _bs_cb(block_size) {

    _blocks = 1;
    _o_blocks = 1;
    _l = grid->FullDimensions();
    _l_cb = grid->LocalDimensions();
    _l_cb_o = grid->_rdimensions;

    _cf_size = 1;
    _block_sites = 1;
    for (int i=0;i<_nd;i++) {
      _l[i] /= grid->_processors[i];

      assert(!(_l[i] % _bs[i])); // lattice must accommodate choice of blocksize

      int r = _l[i] / _l_cb[i];
      assert(!(_bs[i] % r)); // checkerboarding must accommodate choice of blocksize
      _bs_cb[i] = _bs[i] / r;
      _block_sites *= _bs_cb[i];
      _nb[i] = _l[i] / _bs[i];
      _nb_o[i] = _nb[i] / _grid->_simd_layout[i];
      if (_nb[i] % _grid->_simd_layout[i]) { // simd must accommodate choice of blocksize
	std::cout << GridLogMessage << "Problem: _nb[" << i << "] = " << _nb[i] << " _grid->_simd_layout[" << i << "] = " << _grid->_simd_layout[i] << std::endl;
	assert(0);
      }
      _blocks *= _nb[i];
      _o_blocks *= _nb_o[i];
      _cf_size *= _l[i];
    }

    _cf_size *= 12 / 2;
    _cf_block_size = _cf_size / _blocks;
    _cf_o_block_size = _cf_size / _o_blocks;

    std::cout << GridLogMessage << "BlockedGrid:" << std::endl;
    std::cout << GridLogMessage << " _l     = " << _l << std::endl;
    std::cout << GridLogMessage << " _l_cb     = " << _l_cb << std::endl;
    std::cout << GridLogMessage << " _l_cb_o     = " << _l_cb_o << std::endl;
    std::cout << GridLogMessage << " _bs    = " << _bs << std::endl;
    std::cout << GridLogMessage << " _bs_cb    = " << _bs_cb << std::endl;

    std::cout << GridLogMessage << " _nb    = " << _nb << std::endl;
    std::cout << GridLogMessage << " _nb_o    = " << _nb_o << std::endl;
    std::cout << GridLogMessage << " _blocks = " << _blocks << std::endl;
    std::cout << GridLogMessage << " _o_blocks = " << _o_blocks << std::endl;
    std::cout << GridLogMessage << " sizeof(vCoeff_t) = " << sizeof(vCoeff_t) << std::endl;
    std::cout << GridLogMessage << " _cf_size = " << _cf_size << std::endl;
    std::cout << GridLogMessage << " _cf_block_size = " << _cf_block_size << std::endl;
    std::cout << GridLogMessage << " _block_sites = " << _block_sites << std::endl;
    std::cout << GridLogMessage << " _grid->oSites() = " << _grid->oSites() << std::endl;

    //    _grid->Barrier();
    //abort();
  }

    void block_to_coor(int b, std::vector<int>& x0) {

      std::vector<int> bcoor;
      bcoor.resize(_nd);
      x0.resize(_nd);
      assert(b < _o_blocks);
      Lexicographic::CoorFromIndex(bcoor,b,_nb_o);
      int i;

      for (i=0;i<_nd;i++) {
	x0[i] = bcoor[i]*_bs_cb[i];
      }

      //std::cout << GridLogMessage << "Map block b -> " << x0 << std::endl;

    }

    void block_site_to_o_coor(const std::vector<int>& x0, std::vector<int>& coor, int i) {
      Lexicographic::CoorFromIndex(coor,i,_bs_cb);
      for (int j=0;j<_nd;j++)
	coor[j] += x0[j];
    }

    int block_site_to_o_site(const std::vector<int>& x0, int i) {
      std::vector<int> coor;  coor.resize(_nd);
      block_site_to_o_coor(x0,coor,i);
      Lexicographic::IndexFromCoor(coor,i,_l_cb_o);
      return i;
    }

    vCoeff_t block_sp(int b, const Field& x, const Field& y) {

      std::vector<int> x0;
      block_to_coor(b,x0);

      vCoeff_t ret = 0.0;
      for (int i=0;i<_block_sites;i++) { // only odd sites
	int ss = block_site_to_o_site(x0,i);
	ret += TensorRemove(innerProduct(x._odata[ss],y._odata[ss]));
      }

      return ret;

    }

    vCoeff_t block_sp(int b, const Field& x, const std::vector< ComplexD >& y) {

      std::vector<int> x0;
      block_to_coor(b,x0);

      constexpr int nsimd = sizeof(vCoeff_t) / sizeof(Coeff_t);
      int lsize = _cf_o_block_size / _block_sites;

      std::vector< ComplexD > ret(nsimd);
      for (int i=0;i<nsimd;i++)
	ret[i] = 0.0;

      for (int i=0;i<_block_sites;i++) { // only odd sites
	int ss = block_site_to_o_site(x0,i);

	int n = lsize / nsimd;
	for (int l=0;l<n;l++) {
	  for (int j=0;j<nsimd;j++) {
	    int t = lsize * i + l*nsimd + j;

	    ret[j] += conjugate(((Coeff_t*)&x._odata[ss]._internal)[l*nsimd + j]) * y[t];
	  }
	}
      }

      vCoeff_t vret;
      for (int i=0;i<nsimd;i++)
	((Coeff_t*)&vret)[i] = (Coeff_t)ret[i];

      return vret;

    }

    template<class T>
      void vcaxpy(iScalar<T>& r,const vCoeff_t& a,const iScalar<T>& x,const iScalar<T>& y) {
      vcaxpy(r._internal,a,x._internal,y._internal);
    }

    template<class T,int N>
      void vcaxpy(iVector<T,N>& r,const vCoeff_t& a,const iVector<T,N>& x,const iVector<T,N>& y) {
      for (int i=0;i<N;i++)
	vcaxpy(r._internal[i],a,x._internal[i],y._internal[i]);
    }

    void vcaxpy(vCoeff_t& r,const vCoeff_t& a,const vCoeff_t& x,const vCoeff_t& y) {
      r = a*x + y;
    }

    void block_caxpy(int b, Field& ret, const vCoeff_t& a, const Field& x, const Field& y) {

      std::vector<int> x0;
      block_to_coor(b,x0);

      for (int i=0;i<_block_sites;i++) { // only odd sites
	int ss = block_site_to_o_site(x0,i);
	vcaxpy(ret._odata[ss],a,x._odata[ss],y._odata[ss]);
      }

    }

    void block_caxpy(int b, std::vector< ComplexD >& ret, const vCoeff_t& a, const Field& x, const std::vector< ComplexD >& y) {
      std::vector<int> x0;
      block_to_coor(b,x0);

      constexpr int nsimd = sizeof(vCoeff_t) / sizeof(Coeff_t);
      int lsize = _cf_o_block_size / _block_sites;

      for (int i=0;i<_block_sites;i++) { // only odd sites
	int ss = block_site_to_o_site(x0,i);

	int n = lsize / nsimd;
	for (int l=0;l<n;l++) {
	  vCoeff_t r = a* ((vCoeff_t*)&x._odata[ss]._internal)[l];

	  for (int j=0;j<nsimd;j++) {
	    int t = lsize * i + l*nsimd + j;
	    ret[t] = y[t] + ((Coeff_t*)&r)[j];
	  }
	}
      }

    }

    void block_set(int b, Field& ret, const std::vector< ComplexD >& x) {
      std::vector<int> x0;
      block_to_coor(b,x0);

      int lsize = _cf_o_block_size / _block_sites;

      for (int i=0;i<_block_sites;i++) { // only odd sites
	int ss = block_site_to_o_site(x0,i);

	for (int l=0;l<lsize;l++)
	  ((Coeff_t*)&ret._odata[ss]._internal)[l] = (Coeff_t)x[lsize * i + l]; // convert precision
      }

    }

    void block_get(int b, const Field& ret, std::vector< ComplexD >& x) {
      std::vector<int> x0;
      block_to_coor(b,x0);

      int lsize = _cf_o_block_size / _block_sites;

      for (int i=0;i<_block_sites;i++) { // only odd sites
	int ss = block_site_to_o_site(x0,i);

	for (int l=0;l<lsize;l++)
	  x[lsize * i + l] = (ComplexD)((Coeff_t*)&ret._odata[ss]._internal)[l];
      }

    }

    template<class T>
    void vcscale(iScalar<T>& r,const vCoeff_t& a,const iScalar<T>& x) {
      vcscale(r._internal,a,x._internal);
    }

    template<class T,int N>
    void vcscale(iVector<T,N>& r,const vCoeff_t& a,const iVector<T,N>& x) {
      for (int i=0;i<N;i++)
	vcscale(r._internal[i],a,x._internal[i]);
    }

    void vcscale(vCoeff_t& r,const vCoeff_t& a,const vCoeff_t& x) {
      r = a*x;
    }

    void block_cscale(int b, const vCoeff_t& a, Field& ret) {

      std::vector<int> x0;
      block_to_coor(b,x0);
      
      for (int i=0;i<_block_sites;i++) { // only odd sites
	int ss = block_site_to_o_site(x0,i);
	vcscale(ret._odata[ss],a,ret._odata[ss]);
      }
    }

    void getCanonicalBlockOffset(int cb, std::vector<int>& x0) {
      const int ndim = 5;
      assert(_nb.size() == ndim);
      std::vector<int> _nbc = { _nb[1], _nb[2], _nb[3], _nb[4], _nb[0] };
      std::vector<int> _bsc = { _bs[1], _bs[2], _bs[3], _bs[4], _bs[0] };
      x0.resize(ndim);

      assert(cb >= 0);
      assert(cb < _nbc[0]*_nbc[1]*_nbc[2]*_nbc[3]*_nbc[4]);

      Lexicographic::CoorFromIndex(x0,cb,_nbc);
      int i;

      for (i=0;i<ndim;i++) {
	x0[i] *= _bsc[i];
      }

      //if (cb < 2)
      //	std::cout << GridLogMessage << "Map: " << cb << " To: " << x0 << std::endl;
    }

    void pokeBlockOfVectorCanonical(int cb,Field& v,const std::vector<float>& buf) {
      std::vector<int> _bsc = { _bs[1], _bs[2], _bs[3], _bs[4], _bs[0] };
      std::vector<int> ldim = v._grid->LocalDimensions();
      std::vector<int> cldim = { ldim[1], ldim[2], ldim[3], ldim[4], ldim[0] };
      const int _nbsc = _bs_cb[0]*_bs_cb[1]*_bs_cb[2]*_bs_cb[3]*_bs_cb[4];
      // take canonical block cb of v and put it in canonical ordering in buf
      std::vector<int> cx0;
      getCanonicalBlockOffset(cb,cx0);

#pragma omp parallel
      {
	std::vector<int> co0,cl0;
	co0=cx0; cl0=cx0;

#pragma omp for
	for (int i=0;i<_nbsc;i++) {
	  Lexicographic::CoorFromIndex(co0,2*i,_bsc); // 2* for eo
	  for (int j=0;j<(int)_bsc.size();j++)
	    cl0[j] = cx0[j] + co0[j];
	  
	  std::vector<int> l0 = { cl0[4], cl0[0], cl0[1], cl0[2], cl0[3] };
	  int oi = v._grid->oIndex(l0);
	  int ii = v._grid->iIndex(l0);
	  int lti = i;

	  //if (cb < 2 && i<2)
	  //  std::cout << GridLogMessage << "Map: " << cb << ", " << i << " To: " << cl0 << ", " << cx0 << ", " << oi << ", " << ii << std::endl;
	  
	  for (int s=0;s<4;s++)
	    for (int c=0;c<3;c++) {
	      Coeff_t& ld = ((Coeff_t*)&v._odata[oi]._internal._internal[s]._internal[c])[ii];
	      int ti = 12*lti + 3*s + c;
	      ld = Coeff_t(buf[2*ti+0], buf[2*ti+1]);
	    }
	}
      }
    }

    void peekBlockOfVectorCanonical(int cb,const Field& v,std::vector<float>& buf) {
      std::vector<int> _bsc = { _bs[1], _bs[2], _bs[3], _bs[4], _bs[0] };
      std::vector<int> ldim = v._grid->LocalDimensions();
      std::vector<int> cldim = { ldim[1], ldim[2], ldim[3], ldim[4], ldim[0] };
      const int _nbsc = _bs_cb[0]*_bs_cb[1]*_bs_cb[2]*_bs_cb[3]*_bs_cb[4];
      // take canonical block cb of v and put it in canonical ordering in buf
      std::vector<int> cx0;
      getCanonicalBlockOffset(cb,cx0);

      buf.resize(_cf_block_size * 2);

#pragma omp parallel
      {
	std::vector<int> co0,cl0;
	co0=cx0; cl0=cx0;

#pragma omp for
	for (int i=0;i<_nbsc;i++) {
	  Lexicographic::CoorFromIndex(co0,2*i,_bsc); // 2* for eo
	  for (int j=0;j<(int)_bsc.size();j++)
	    cl0[j] = cx0[j] + co0[j];
	  
	  std::vector<int> l0 = { cl0[4], cl0[0], cl0[1], cl0[2], cl0[3] };
	  int oi = v._grid->oIndex(l0);
	  int ii = v._grid->iIndex(l0);
	  int lti = i;
	  
	  //if (cb < 2 && i<2)
	  //  std::cout << GridLogMessage << "Map: " << cb << ", " << i << " To: " << cl0 << ", " << cx0 << ", " << oi << ", " << ii << std::endl;

	  for (int s=0;s<4;s++)
	    for (int c=0;c<3;c++) {
	      Coeff_t& ld = ((Coeff_t*)&v._odata[oi]._internal._internal[s]._internal[c])[ii];
	      int ti = 12*lti + 3*s + c;
	      buf[2*ti+0] = ld.real();
	      buf[2*ti+1] = ld.imag();
	    }
	}
      }
    }

    int globalToLocalCanonicalBlock(int slot,const std::vector<int>& src_nodes,int nb) {
      // processor coordinate
      int _nd = (int)src_nodes.size();
      std::vector<int> _src_nodes = src_nodes;
      std::vector<int> pco(_nd);
      Lexicographic::CoorFromIndex(pco,slot,_src_nodes);
      std::vector<int> cpco = { pco[1], pco[2], pco[3], pco[4], pco[0] };

      // get local block
      std::vector<int> _nbc = { _nb[1], _nb[2], _nb[3], _nb[4], _nb[0] };
      assert(_nd == 5);
      std::vector<int> c_src_local_blocks(_nd);
      for (int i=0;i<_nd;i++) {
	assert(_grid->_fdimensions[i] % (src_nodes[i] * _bs[i]) == 0);
	c_src_local_blocks[(i+4) % 5] = _grid->_fdimensions[i] / src_nodes[i] / _bs[i];
      }
      std::vector<int> cbcoor(_nd); // coordinate of block in slot in canonical form
      Lexicographic::CoorFromIndex(cbcoor,nb,c_src_local_blocks);

      // cpco, cbcoor
      std::vector<int> clbcoor(_nd);
      for (int i=0;i<_nd;i++) {
	int cgcoor = cpco[i] * c_src_local_blocks[i] + cbcoor[i]; // global block coordinate
	int pcoor = cgcoor / _nbc[i]; // processor coordinate in my Grid
	int tpcoor = _grid->_processor_coor[(i+1)%5];
	if (pcoor != tpcoor)
	  return -1;
	clbcoor[i] = cgcoor - tpcoor * _nbc[i]; // canonical local block coordinate for canonical dimension i
      }

      int lnb;
      Lexicographic::IndexFromCoor(clbcoor,lnb,_nbc);
      //std::cout << "Mapped slot = " << slot << " nb = " << nb << " to " << lnb << std::endl;
      return lnb;
    }


 };

}
