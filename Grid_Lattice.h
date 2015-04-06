#ifndef GRID_LATTICE_H
#define GRID_LATTICE_H

#include "Grid.h"



namespace Grid {

// Permute the pointers 32bitx16 = 512
static int permute_map[4][16] = { 
  { 1,0,3,2,5,4,7,6,9,8,11,10,13,12,15,14},
  { 2,3,0,1,6,7,4,5,10,11,8,9,14,15,12,13},
  { 4,5,6,7,0,1,2,3,12,13,14,15,8,9,10,11},
  { 9,10,11,12,13,14,15,0,1,2,3,4,5,6,7,8}
};


template<class vobj>
class Lattice
{
public:
    GridBase *_grid;
    int checkerboard;
    std::vector<vobj,alignedAllocator<vobj> > _odata;
    typedef typename vobj::scalar_type scalar_type;
    typedef typename vobj::vector_type vector_type;
public:


    Lattice(GridBase *grid) : _grid(grid) {
        _odata.reserve(_grid->oSites());
        assert((((uint64_t)&_odata[0])&0xF) ==0);
        checkerboard=0;
    }
    

#include <Grid_cshift.h>
    
    // overloading Grid::conformable but no conformable in Grid ...?:w
    template<class obj1,class obj2>
    friend void conformable(const Lattice<obj1> &lhs,const Lattice<obj2> &rhs);

    // Performance difference between operator * and mult is troubling.
    // Auto move constructor seems to lose surprisingly much.

    // Site wise binary operations
    // We eliminate a temporary object assignment if use the mult,add,sub routines.
    // For the operator versions we rely on move constructor to eliminate the
    // vector copy back.
    template<class obj1,class obj2,class obj3>
    friend void mult(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs);

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
    
    // FIXME for debug; deprecate this
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
    };

    friend void gaussian(Lattice<vobj> &l){
        // Zero mean, unit variance.
        std::normal_distribution<double> distribution(0.0,1.0);
        Real *v_ptr = (Real *)&l._odata[0];
        size_t v_len = l._grid->oSites()*sizeof(vobj);
        size_t d_len = v_len/sizeof(Real);
        
        // Not a parallel RNG. Could make up some seed per 4d site, seed
        // per hypercube type scheme.
        for(int i=0;i<d_len;i++){
	  v_ptr[i]= drand48();
        }
    };


    // Unary functions and Unops
    // Unary negation
    friend inline Lattice<vobj> operator -(const Lattice<vobj> &r) {
        Lattice<vobj> ret(r._grid);
#pragma omp parallel for
        for(int ss=0;ss<r._grid->oSites();ss++){
            ret._odata[ss]= -r._odata[ss];
        }
        return ret;
    }
    // *=,+=,-= operators
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
    
    inline friend Lattice<iScalar<vComplex> > _trace(const Lattice<vobj> &lhs){
        Lattice<iScalar<vComplex> > ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = trace(lhs._odata[ss]);
        }
        return ret;
    };
    
    inline friend Lattice<iScalar<iScalar< vComplex > > > trace(const Lattice<vobj> &lhs){
        Lattice<iScalar< iScalar<vComplex> >  > ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = trace(lhs._odata[ss]);
        }
        return ret;
    };

    
    inline friend Lattice<vobj> adj(const Lattice<vobj> &lhs){
        Lattice<vobj> ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = adj(lhs._odata[ss]);
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
        for(int ss=0;ss<rhs._grid->oSites(); ss++){
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

}
#endif
