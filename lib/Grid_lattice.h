#ifndef GRID_LATTICE_H
#define GRID_LATTICE_H

namespace Grid {

// TODO: Indexing ()
//       mac,real,imag
//
// Functionality required:
//     -=,+=,*=,()
//     add,+,sub,-,mult,mac,*
//     adj,conj
//     real,imag
//     transpose,transposeIndex  
//     trace,traceIndex
//     peekIndex
//     innerProduct,outerProduct,
//     localNorm2
//     localInnerProduct
//     

extern int GridCshiftPermuteMap[4][16];

template<class vobj>
class Lattice
{
public:
    GridBase *_grid;
    int checkerboard;
    std::vector<vobj,alignedAllocator<vobj> > _odata;
public:

    typedef typename vobj::scalar_type scalar_type;
    typedef typename vobj::vector_type vector_type;

    Lattice(GridBase *grid) : _grid(grid) {
        _odata.reserve(_grid->oSites());
        assert((((uint64_t)&_odata[0])&0xF) ==0);
        checkerboard=0;
    }

#include <Grid_cshift.h>
   
    template<class obj1,class obj2>
    friend void conformable(const Lattice<obj1> &lhs,const Lattice<obj2> &rhs);

    // FIXME Performance difference between operator * and mult is troubling.
    // Auto move constructor seems to lose surprisingly much.

    // Site wise binary operations
    // We eliminate a temporary object assignment if use the mult,add,sub routines.
    // For the operator versions we rely on move constructor to eliminate the
    // vector copy back.
    template<class obj1,class obj2,class obj3>
    friend void mult(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs);

    template<class obj1,class obj2,class obj3>
    friend void mac(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs);

    template<class obj1,class obj2,class obj3>
    friend void sub(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs);

    template<class obj1,class obj2,class obj3>
    friend void add(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs);

    
    friend void axpy(Lattice<vobj> &ret,double a,const Lattice<vobj> &lhs,const Lattice<vobj> &rhs){
        conformable(lhs,rhs);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            axpy(&ret._odata[ss],a,&lhs._odata[ss],&rhs._odata[ss]);
        }
    }
    friend void axpy(Lattice<vobj> &ret,std::complex<double> a,const Lattice<vobj> &lhs,const Lattice<vobj> &rhs){
        conformable(lhs,rhs);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            axpy(&ret._odata[ss],a,&lhs._odata[ss],&rhs._odata[ss]);
        }
    }
    inline friend Lattice<vobj> operator / (const Lattice<vobj> &lhs,const Lattice<vobj> &rhs){
        conformable(lhs,rhs);
        Lattice<vobj> ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = lhs._odata[ss]/rhs._odata[ss];
        }
        return ret;
    };

    template<class sobj>
    inline Lattice<vobj> & operator = (const sobj & r){
#pragma omp parallel for
        for(int ss=0;ss<_grid->oSites();ss++){
            this->_odata[ss]=r;
        }
        return *this;
    }
    
    // Poke a scalar object into the SIMD array
    template<class sobj>
    friend void pokeSite(const sobj &s,Lattice<vobj> &l,std::vector<int> &site){

      GridBase *grid=l._grid;

      typedef typename vobj::scalar_type scalar_type;
      typedef typename vobj::vector_type vector_type;

      int Nsimd = grid->Nsimd();

      assert( l.checkerboard== l._grid->CheckerBoard(site));
      assert( sizeof(sobj)*Nsimd == sizeof(vobj));

      int rank,odx,idx;
      grid->GlobalCoorToRankIndex(rank,odx,idx,site);

      // Optional to broadcast from node 0.
      grid->Broadcast(0,s);

      std::vector<sobj> buf(Nsimd);
      std::vector<scalar_type *> pointers(Nsimd);  
      for(int i=0;i<Nsimd;i++) pointers[i] = (scalar_type *)&buf[i];

      // extract-modify-merge cycle is easiest way and this is not perf critical
      extract(l._odata[odx],pointers);
      
      buf[idx] = s;

      for(int i=0;i<Nsimd;i++) pointers[i] = (scalar_type *)&buf[i];
      merge(l._odata[odx],pointers);

      return;
    };
    
    // Peek a scalar object from the SIMD array
    template<class sobj>
    friend void peekSite(sobj &s,Lattice<vobj> &l,std::vector<int> &site){
        
      GridBase *grid=l._grid;

      typedef typename vobj::scalar_type scalar_type;
      typedef typename vobj::vector_type vector_type;

      int Nsimd = grid->Nsimd();

      assert( l.checkerboard== l._grid->CheckerBoard(site));
      assert( sizeof(sobj)*Nsimd == sizeof(vobj));

      int rank,odx,idx;
      grid->GlobalCoorToRankIndex(rank,odx,idx,site);
      std::vector<sobj> buf(Nsimd);
      std::vector<scalar_type *> pointers(Nsimd);  
      for(int i=0;i<Nsimd;i++) pointers[i] = (scalar_type *)&buf[i];

      extract(l._odata[odx],pointers);
      
      s = buf[idx];
      grid->Broadcast(rank,s);

      return;
    };
    
    // FIXME Randomise; deprecate this
    friend void random(Lattice<vobj> &l){
        Real *v_ptr = (Real *)&l._odata[0];
        size_t v_len = l._grid->oSites()*sizeof(vobj);
        size_t d_len = v_len/sizeof(Real);
	
        for(int i=0;i<d_len;i++){

            v_ptr[i]=drand48();
        }
    };

    // FIXME for debug; deprecate this; made obscelete by 
    // LatticeCoordinate();
    friend void lex_sites(Lattice<vobj> &l){
      Real *v_ptr = (Real *)&l._odata[0];
      size_t o_len = l._grid->oSites();
      size_t v_len = sizeof(vobj)/sizeof(vRealF);
      size_t vec_len = vRealF::Nsimd();

      for(int i=0;i<o_len;i++){
	for(int j=0;j<v_len;j++){
          for(int vv=0;vv<vec_len;vv+=2){
	    v_ptr[i*v_len*vec_len+j*vec_len+vv  ]= i+vv*500;
	    v_ptr[i*v_len*vec_len+j*vec_len+vv+1]= i+vv*500;
	  }
	}}
    }
    
    // FIXME Implement a consistent seed management strategy
    friend void gaussian(Lattice<vobj> &l){
        // Zero mean, unit variance.
        std::normal_distribution<double> distribution(0.0,1.0);
        Real *v_ptr = (Real *)&l._odata[0];
        size_t v_len = l._grid->oSites()*sizeof(vobj);
        size_t d_len = v_len/sizeof(Real);

        for(int i=0;i<d_len;i++){
	  v_ptr[i]= drand48();
        }
    };

    // Unary functions and Unops
    friend inline Lattice<vobj> operator -(const Lattice<vobj> &r) {
        Lattice<vobj> ret(r._grid);
#pragma omp parallel for
        for(int ss=0;ss<r._grid->oSites();ss++){
            ret._odata[ss]= -r._odata[ss];
        }
        return ret;
    }
    // *=,+=,-= operators inherit behvour from correspond */+/- operation
    template<class T>
    inline Lattice<vobj> &operator *=(const T &r) {
        *this = (*this)*r;
        return *this;
    }
    template<class T>
    inline Lattice<vobj> &operator -=(const T &r) {
        *this = (*this)-r;
        return *this;
    }
    template<class T>
    inline Lattice<vobj> &operator +=(const T &r) {
        *this = (*this)+r;
        return *this;
    }
    
    inline friend Lattice<vobj> adj(const Lattice<vobj> &lhs){
        Lattice<vobj> ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = adj(lhs._odata[ss]);
        }
        return ret;
    };

    inline friend Lattice<vobj> transpose(const Lattice<vobj> &lhs){
        Lattice<vobj> ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = transpose(lhs._odata[ss]);
        }
        return ret;
    };


    inline friend Lattice<vobj> conj(const Lattice<vobj> &lhs){
        Lattice<vobj> ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = conj(lhs._odata[ss]);
        }
        return ret;
    };

    // remove and insert a half checkerboard
    friend void pickCheckerboard(int cb,Lattice<vobj> &half,const Lattice<vobj> &full){
      half.checkerboard = cb;
      int ssh=0;
#pragma omp parallel for
      for(int ss=0;ss<full._grid->oSites();ss++){
	std::vector<int> coor;
	int cbos;
	
	full._grid->oCoorFromOindex(coor,ss);
	cbos=half._grid->CheckerBoard(coor);
	
	if (cbos==cb) {

	  half._odata[ssh] = full._odata[ss];
	  ssh++;
	}
      }
    }
    friend void setCheckerboard(Lattice<vobj> &full,const Lattice<vobj> &half){
      int cb = half.checkerboard;
      int ssh=0;
#pragma omp parallel for
      for(int ss=0;ss<full._grid->oSites();ss++){
	std::vector<int> coor;
	int cbos;
	
	full._grid->oCoorFromOindex(coor,ss);
	cbos=half._grid->CheckerBoard(coor);

	if (cbos==cb) {
	  full._odata[ss]=half._odata[ssh];
	  ssh++;
	}
      }
    }
}; // class Lattice

    template<class obj1,class obj2>
    void conformable(const Lattice<obj1> &lhs,const Lattice<obj2> &rhs)
    {
        assert(lhs._grid == rhs._grid);
        assert(lhs.checkerboard == rhs.checkerboard);
    }

    template<class obj1,class obj2,class obj3>
    void mult(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
        conformable(lhs,rhs);
	uint32_t vec_len = lhs._grid->oSites();
#pragma omp parallel for
        for(int ss=0;ss<vec_len;ss++){
	  mult(&ret._odata[ss],&lhs._odata[ss],&rhs._odata[ss]);
        }
    }

    template<class obj1,class obj2,class obj3>
    void mac(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
        conformable(lhs,rhs);
	uint32_t vec_len = lhs._grid->oSites();
#pragma omp parallel for
        for(int ss=0;ss<vec_len;ss++){
	  mac(&ret._odata[ss],&lhs._odata[ss],&rhs._odata[ss]);
        }
    }

    template<class obj1,class obj2,class obj3>
    void sub(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
        conformable(lhs,rhs);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            sub(&ret._odata[ss],&lhs._odata[ss],&rhs._odata[ss]);
        }
    }
    template<class obj1,class obj2,class obj3>
    void add(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
        conformable(lhs,rhs);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            add(&ret._odata[ss],&lhs._odata[ss],&rhs._odata[ss]);
        }
    }

    // Lattice BinOp Lattice,
    template<class left,class right>
    inline auto operator * (const Lattice<left> &lhs,const Lattice<right> &rhs)-> Lattice<decltype(lhs._odata[0]*rhs._odata[0])>
    {
        //NB mult performs conformable check. Do not reapply here for performance.
      Lattice<decltype(lhs._odata[0]*rhs._odata[0])> ret(rhs._grid);
      mult(ret,lhs,rhs);
      return ret;
    }
    template<class left,class right>
    inline auto operator + (const Lattice<left> &lhs,const Lattice<right> &rhs)-> Lattice<decltype(lhs._odata[0]*rhs._odata[0])>
    {
        //NB mult performs conformable check. Do not reapply here for performance.
        Lattice<decltype(lhs._odata[0]*rhs._odata[0])> ret(rhs._grid);
        add(ret,lhs,rhs);
        return ret;
    }
    template<class left,class right>
    inline auto operator - (const Lattice<left> &lhs,const Lattice<right> &rhs)-> Lattice<decltype(lhs._odata[0]*rhs._odata[0])>
    {
        //NB mult performs conformable check. Do not reapply here for performance.
        Lattice<decltype(lhs._odata[0]*rhs._odata[0])> ret(rhs._grid);
        sub(ret,lhs,rhs);
        return ret;
    }

    // Scalar BinOp Lattice ;generate return type
    template<class left,class right>
    inline auto operator * (const left &lhs,const Lattice<right> &rhs) -> Lattice<decltype(lhs*rhs._odata[0])>
    {
        Lattice<decltype(lhs*rhs._odata[0])> ret(rhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<rhs._grid->oSites(); ss++){
            ret._odata[ss]=lhs*rhs._odata[ss];
        }
        return ret;
    }
    template<class left,class right>
    inline auto operator + (const left &lhs,const Lattice<right> &rhs) -> Lattice<decltype(lhs*rhs._odata[0])>
    {
        Lattice<decltype(lhs*rhs._odata[0])> ret(rhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<rhs._grid->oSites(); ss++){
            ret._odata[ss]=lhs+rhs._odata[ss];
        }
        return ret;
    }
    template<class left,class right>
    inline auto operator - (const left &lhs,const Lattice<right> &rhs) -> Lattice<decltype(lhs*rhs._odata[0])>
    {
        Lattice<decltype(lhs*rhs._odata[0])> ret(rhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<rhs._grid->oSites(); ss++){
            ret._odata[ss]=lhs-rhs._odata[ss];
        }
        return ret;
    }
    template<class left,class right>
    inline auto operator * (const Lattice<left> &lhs,const right &rhs) -> Lattice<decltype(lhs._odata[0]*rhs)>
    {
        Lattice<decltype(lhs._odata[0]*rhs)> ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites(); ss++){
            ret._odata[ss]=lhs._odata[ss]*rhs;
        }
        return ret;
    }
    template<class left,class right>
    inline auto operator + (const Lattice<left> &lhs,const right &rhs) -> Lattice<decltype(lhs._odata[0]*rhs)>
    {
        Lattice<decltype(lhs._odata[0]*rhs)> ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<rhs._grid->oSites(); ss++){
            ret._odata[ss]=lhs._odata[ss]+rhs;
        }
        return ret;
    }
    template<class left,class right>
    inline auto operator - (const Lattice<left> &lhs,const right &rhs) -> Lattice<decltype(lhs._odata[0]*rhs)>
    {
        Lattice<decltype(lhs._odata[0]*rhs)> ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<rhs._grid->oSites(); ss++){
            ret._odata[ss]=lhs._odata[ss]-rhs;
        }
        return ret;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Trace
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    template<class vobj>
    inline auto trace(const Lattice<vobj> &lhs)
      -> Lattice<decltype(trace(lhs._odata[0]))>
    {
      Lattice<decltype(trace(lhs._odata[0]))> ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = trace(lhs._odata[ss]);
        }
        return ret;
    };
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Index level dependent operations
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    template<int Index,class vobj>
    inline auto traceIndex(const Lattice<vobj> &lhs)
      -> Lattice<decltype(traceIndex<Index>(lhs._odata[0]))>
    {
      Lattice<decltype(traceIndex<Index>(lhs._odata[0]))> ret(lhs._grid);
#pragma omp parallel for
      for(int ss=0;ss<lhs._grid->oSites();ss++){
	ret._odata[ss] = traceIndex<Index>(lhs._odata[ss]);
      }
      return ret;
    };
    template<int Index,class vobj>
    inline auto transposeIndex(const Lattice<vobj> &lhs)
      -> Lattice<decltype(transposeIndex<Index>(lhs._odata[0]))>
    {
      Lattice<decltype(transposeIndex<Index>(lhs._odata[0]))> ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = transposeIndex<Index>(lhs._odata[ss]);
        }
        return ret;
    };

    // Fixme; this is problematic since the number of args is variable and 
    // may mismatch...
    template<int Index,class vobj>
    inline auto peekIndex(const Lattice<vobj> &lhs)
      -> Lattice<decltype(peekIndex<Index>(lhs._odata[0]))>
    {
      Lattice<decltype(peekIndex<Index>(lhs._odata[0]))> ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = peekIndex<Index>(lhs._odata[ss]);
        }
        return ret;
    };
    template<int Index,class vobj>
      inline auto peekIndex(const Lattice<vobj> &lhs,int i)
      -> Lattice<decltype(peekIndex<Index>(lhs._odata[0],i))>
    {
      Lattice<decltype(peekIndex<Index>(lhs._odata[0],i))> ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
	  ret._odata[ss] = peekIndex<Index>(lhs._odata[ss],i);
        }
        return ret;
    };
    template<int Index,class vobj>
      inline auto peekIndex(const Lattice<vobj> &lhs,int i,int j)
      -> Lattice<decltype(peekIndex<Index>(lhs._odata[0],i,j))>
    {
      Lattice<decltype(peekIndex<Index>(lhs._odata[0],i,j))> ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
	  ret._odata[ss] = peekIndex<Index>(lhs._odata[ss],i,j);
        }
        return ret;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Poke internal indices of a Lattice object
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    template<int Index,class vobj> inline
    void pokeIndex(Lattice<vobj> &lhs,const Lattice<decltype(peekIndex<Index>(lhs._odata[0]))> & rhs)
    {
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
	  pokeIndex<Index>(lhs._odata[ss],rhs._odata[ss]);
	}      
    }
    template<int Index,class vobj> inline
    void pokeIndex(Lattice<vobj> &lhs,const Lattice<decltype(peekIndex<Index>(lhs._odata[0],0))> & rhs,int i)
    {
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
	  pokeIndex<Index>(lhs._odata[ss],rhs._odata[ss],i);
	}      
    }
    template<int Index,class vobj> inline
    void pokeIndex(Lattice<vobj> &lhs,const Lattice<decltype(peekIndex<Index>(lhs._odata[0],0,0))> & rhs,int i,int j)
    {
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
	  pokeIndex<Index>(lhs._odata[ss],rhs._odata[ss],i,j);
	}      
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Reduction operations
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    template<class vobj>
    inline RealD norm2(const Lattice<vobj> &arg){

      typedef typename vobj::scalar_type scalar;
      typedef typename vobj::vector_type vector;
      decltype(innerProduct(arg._odata[0],arg._odata[0])) vnrm=zero;
      scalar nrm;
      //FIXME make this loop parallelisable
      vnrm=zero;
      for(int ss=0;ss<arg._grid->oSites(); ss++){
	vnrm = vnrm + innerProduct(arg._odata[ss],arg._odata[ss]);
      }
      vector vvnrm =TensorRemove(vnrm) ;
      nrm = Reduce(vvnrm);
      arg._grid->GlobalSum(nrm);
      return real(nrm);
    }

    template<class vobj>
    inline auto innerProduct(const Lattice<vobj> &left,const Lattice<vobj> &right) ->decltype(innerProduct(left._odata[0],right._odata[0]))
    {
      typedef typename vobj::scalar_type scalar;
      decltype(innerProduct(left._odata[0],right._odata[0])) vnrm=zero;

      scalar nrm;
      //FIXME make this loop parallelisable
      for(int ss=0;ss<left._grid->oSites(); ss++){
	vnrm = vnrm + innerProduct(left._odata[ss],right._odata[ss]);
      }
      nrm = Reduce(vnrm);
      right._grid->GlobalSum(nrm);
      return nrm;
    }

    /////////////////////////////////////////////////////
    // Non site reduced routines
    /////////////////////////////////////////////////////

    // localNorm2,
    template<class vobj>
    inline auto localNorm2 (const Lattice<vobj> &rhs)-> Lattice<typename vobj::tensor_reduced>
    {
      Lattice<typename vobj::tensor_reduced> ret(rhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<rhs._grid->oSites(); ss++){
	  ret._odata[ss]=innerProduct(rhs._odata[ss],rhs._odata[ss]);
        }
        return ret;
    }
    
    template<class vobj>
    inline auto real(const Lattice<vobj> &z) -> Lattice<decltype(real(z._odata[0]))>
    {
      Lattice<decltype(real(z._odata[0]))> ret(z._grid);
#pragma omp parallel for
        for(int ss=0;ss<z._grid->oSites();ss++){
            ret._odata[ss] = real(z._odata[ss]);
        }
      return ret;
    }

    template<class vobj>
    inline auto imag(const Lattice<vobj> &z) -> Lattice<decltype(imag(z._odata[0]))>
    {
      Lattice<decltype(imag(z._odata[0]))> ret(z._grid);
#pragma omp parallel for
        for(int ss=0;ss<z._grid->oSites();ss++){
            ret._odata[ss] = imag(z._odata[ss]);
        }
      return ret;
    }

    // localInnerProduct
    template<class vobj>
    inline auto localInnerProduct (const Lattice<vobj> &lhs,const Lattice<vobj> &rhs)
      -> Lattice<typename vobj::tensor_reduced>
    {
      Lattice<typename vobj::tensor_reduced> ret(rhs._grid);
#pragma omp parallel for
      for(int ss=0;ss<rhs._grid->oSites(); ss++){
	ret._odata[ss]=innerProduct(lhs._odata[ss],rhs._odata[ss]);
      }
      return ret;
    }
    
    // outerProduct Scalar x Scalar -> Scalar
    //              Vector x Vector -> Matrix
    template<class ll,class rr>
    inline auto outerProduct (const Lattice<ll> &lhs,const Lattice<rr> &rhs) -> Lattice<decltype(outerProduct(lhs._odata[0],rhs._odata[0]))>
    {
        Lattice<decltype(outerProduct(lhs._odata[0],rhs._odata[0]))> ret(rhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<rhs._grid->oSites(); ss++){
            ret._odata[ss]=outerProduct(lhs._odata[ss],rhs._odata[ss]);
        }
        return ret;
     }


}
#endif
