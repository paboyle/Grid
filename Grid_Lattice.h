#ifndef GRID_LATTICE_H
#define GRID_LATTICE_H

#include "Grid.h"

namespace dpo {

template<class vobj>
class Lattice
{
public:
    Grid *_grid;
    int checkerboard;
    std::vector<vobj,alignedAllocator<vobj> > _odata;

public:
    
    Lattice(Grid *grid) : _grid(grid) {
        _odata.reserve(_grid->oSites());
        if ( ((uint64_t)&_odata[0])&0xF) {
            exit(-1);
        }
        checkerboard=0;
    }
    
    
    // overloading dpo::conformable but no conformable in dpo ...?:w
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

#if 0 
    // Collapse doesn't appear to work the way I think it should in icpc
    friend Lattice<vobj> Cshift(Lattice<vobj> &rhs,int dimension,int shift)
    {
        Lattice<vobj> ret(rhs._grid);
        
        ret.checkerboard = rhs._grid->CheckerBoardDestination(rhs.checkerboard,shift);
        shift = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift);
        int sx,so,o;
        int rd = rhs._grid->_rdimensions[dimension];
        int ld = rhs._grid->_dimensions[dimension];

        // Map to always positive shift.
        shift = (shift+ld)%ld;

        // Work out whether to permute and the permute type
        // ABCDEFGH ->   AE BF CG DH       permute
        // Shift 0       AE BF CG DH       0 0 0 0    ABCDEFGH
        // Shift 1       DH AE BF CG       1 0 0 0    HABCDEFG
        // Shift 2       CG DH AE BF       1 1 0 0    GHABCDEF
        // Shift 3       BF CG DH AE       1 1 1 0    FGHACBDE
        // Shift 4       AE BF CG DH       1 1 1 1    EFGHABCD
        // Shift 5       DH AE BF CG       0 1 1 1    DEFGHABC
        // Shift 6       CG DH AE BF       0 0 1 1    CDEFGHAB
        // Shift 7       BF CG DH AE       0 0 0 1    BCDEFGHA
        int permute_dim =rhs._grid->_layout[dimension]>1 ;
        int permute_type=0;
        for(int d=0;d<dimension;d++)
            if (rhs._grid->_layout[d]>1 ) permute_type++;

        
        // loop over perp slices.
        // Threading considerations:
        //   Need to map thread_num to
        //
        //               x_min,x_max for Loop-A
        //               n_min,n_max for Loop-B
        //               b_min,b_max for Loop-C
        //  In a way that maximally load balances.
        //
        //  Optimal:
        //      There are rd*n_block*block items of work.
        //      These serialise as item "w"
        //      b=w%block; w=w/block
        //      n=w%nblock; x=w/nblock. Perhaps 20 cycles?
        //
        //  Logic:
        //      x_chunk = (rd+thread)/nthreads simply divide work across nodes.
        //
        //      rd=5 , threads = 8;
        //      0 1 2 3 4 5 6 7
        //      0 0 0 1 1 1 1 1
        for(int x=0;x<rd;x++){         // Loop A
            sx = (x-shift+ld)%rd;
            o  = x*rhs._grid->_ostride[dimension];
            so =sx*rhs._grid->_ostride[dimension];

            int permute_slice=0;
            if ( permute_dim ) {
                permute_slice = shift/rd;
                if ( x<shift%rd ) permute_slice = 1-permute_slice;
            }


#if 0
            if ( permute_slice ) {
                
                int internal=sizeof(vobj)/sizeof(vComplex);
                int num =rhs._grid->_slice_block[dimension]*internal;
                
                for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
                    vComplex *optr = (vComplex *)&ret._odata[o];
                    vComplex *iptr = (vComplex *)&rhs._odata[so];
                    for(int b=0;b<num;b++){
                        permute(optr[b],iptr[b],permute_type);
                    }
                    o+=rhs._grid->_slice_stride[dimension];
                    so+=rhs._grid->_slice_stride[dimension];
                }
            } else {
                for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
                    for(int i=0;i<rhs._grid->_slice_block[dimension];i++){
                        ret._odata[o+i]=rhs._odata[so+i];
                    }
                    o+=rhs._grid->_slice_stride[dimension];
                    so+=rhs._grid->_slice_stride[dimension];
                }
            
            }
#else

            if ( permute_slice ) {

	        int internal=sizeof(vobj)/sizeof(vComplex);
                int num =rhs._grid->_slice_block[dimension]*internal;


#pragma omp parallel for collapse(2)
                for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
                    for(int b=0;b<num;b++){
		      vComplex *optr = (vComplex *)&ret._odata[o +n*rhs._grid->_slice_stride[dimension]];
		      vComplex *iptr = (vComplex *)&rhs._odata[so+n*rhs._grid->_slice_stride[dimension]];
                      permute(optr[b],iptr[b],permute_type);
                    }
                }

            } else {


#pragma omp parallel for collapse(2)
                for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
                    for(int i=0;i<rhs._grid->_slice_block[dimension];i++){
		      int oo = o+ n*rhs._grid->_slice_stride[dimension];
		      int soo=so+ n*rhs._grid->_slice_stride[dimension];
		      ret._odata[oo+i]=rhs._odata[soo+i];
                    }
                }
            
            }
#endif
        }
        return ret;
    }
#else
    friend Lattice<vobj> Cshift(Lattice<vobj> &rhs,int dimension,int shift)
    {
        Lattice<vobj> ret(rhs._grid);
        
        ret.checkerboard = rhs._grid->CheckerBoardDestination(rhs.checkerboard,shift);
        shift = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift);
        int rd = rhs._grid->_rdimensions[dimension];
        int ld = rhs._grid->_dimensions[dimension];
        
        // Map to always positive shift.
        shift = (shift+ld)%ld;
        
        // Work out whether to permute and the permute type
        // ABCDEFGH ->   AE BF CG DH       permute
        // Shift 0       AE BF CG DH       0 0 0 0    ABCDEFGH
        // Shift 1       DH AE BF CG       1 0 0 0    HABCDEFG
        // Shift 2       CG DH AE BF       1 1 0 0    GHABCDEF
        // Shift 3       BF CG DH AE       1 1 1 0    FGHACBDE
        // Shift 4       AE BF CG DH       1 1 1 1    EFGHABCD
        // Shift 5       DH AE BF CG       0 1 1 1    DEFGHABC
        // Shift 6       CG DH AE BF       0 0 1 1    CDEFGHAB
        // Shift 7       BF CG DH AE       0 0 0 1    BCDEFGHA
        int permute_dim =rhs._grid->_layout[dimension]>1 ;
        int permute_type=0;
        for(int d=0;d<dimension;d++)
            if (rhs._grid->_layout[d]>1 ) permute_type++;
        
        
        // loop over all work
        int work =rd*rhs._grid->_slice_nblock[dimension]*rhs._grid->_slice_block[dimension];

#pragma omp parallel for
        for(int ww=0;ww<work;ww++){

            
            // can optimise this if know w moves congtiguously for a given thread.
            // b=(b+1);
            // if (b==_slice_block) {b=0; n=n+1;}
            // if (n==_slice_nblock) { n=0; x=x+1}
            //
            // Perhaps a five cycle iterator, or so.
            int w=ww;
            int b = w%rhs._grid->_slice_block[dimension] ; w=w/rhs._grid->_slice_block[dimension];
            int n = w%rhs._grid->_slice_nblock[dimension]; w=w/rhs._grid->_slice_nblock[dimension];
            int x = w;

	    int sx,so,o;
            sx = (x-shift+ld)%rd;
            o  = x*rhs._grid->_ostride[dimension]+n*rhs._grid->_slice_stride[dimension]; // common sub expression alert.
            so =sx*rhs._grid->_ostride[dimension]+n*rhs._grid->_slice_stride[dimension];
            
            int permute_slice=0;
            if ( permute_dim ) {
                permute_slice = shift/rd;
                if ( x<shift%rd ) permute_slice = 1-permute_slice;
            }
            
            if ( permute_slice ) {
                
                int internal=sizeof(vobj)/sizeof(vComplex);
                vComplex *optr = (vComplex *)&ret._odata[o+b];
                vComplex *iptr = (vComplex *)&rhs._odata[so+b];
		const char *pf = (const char *)iptr;
		for(int i=0;i<sizeof(vobj);i+=64){
		  _mm_prefetch(pf+i,_MM_HINT_T0);
		}

                for(int i=0;i<internal;i++){
                    permute(optr[i],iptr[i],permute_type);
                }
            } else {
	        const char *pf = (const char *) &rhs._odata[so+b];
		for(int i=0;i<sizeof(vobj);i+=64){
		  _mm_prefetch(pf+i,_MM_HINT_T0);
		}
                ret._odata[o+b]=rhs._odata[so+b];
            }
        }
        return ret;
    }
#endif

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

        if ( l._grid.checkerboard != l._grid->Checkerboard(site)){
            printf("Poking wrong checkerboard\n");
            exit(EXIT_FAILURE);
        }

        int o_index = l._grid->oIndex(site);
        int i_index = l._grid->iIndex(site);
        
        Real *v_ptr = (Real *)&l._odata[o_index];
        Real *s_ptr = (Real *)&s;
        v_ptr = v_ptr + 2*i_index;
        
        for(int i=0;i<sizeof(sobj);i+=2*sizeof(Real)){
            v_ptr[0] = s_ptr[0];
            v_ptr[1] = s_ptr[1];
            v_ptr+=2*vComplex::Nsimd();
            s_ptr+=2;
        }
        return;
    };
    
    
    // Peek a scalar object from the SIMD array
    template<class sobj>
    friend void peekSite(sobj &s,const Lattice<vobj> &l,std::vector<int> &site){
        
        // FIXME : define exceptions set and throw up.
        if ( l.checkerboard != l._grid->CheckerBoard(site)){
            printf("Peeking wrong checkerboard\n");
            exit(EXIT_FAILURE);
        }
        int o_index = l._grid->oIndex(site);
        int i_index = l._grid->iIndex(site);
        
        Real *v_ptr = (Real *)&l._odata[o_index];
        Real *s_ptr = (Real *)&s;
        v_ptr = v_ptr + 2*i_index;
        
        for(int i=0;i<sizeof(sobj);i+=2*sizeof(Real)){
            s_ptr[0] = v_ptr[0];
            s_ptr[1] = v_ptr[1];
            v_ptr+=2*vComplex::Nsimd();
            s_ptr+=2;
        }
        return;
    };
    
    // Randomise
    friend void random(Lattice<vobj> &l){

        Real *v_ptr = (Real *)&l._odata[0];
        size_t v_len = l._grid->oSites()*sizeof(vobj);
        size_t d_len = v_len/sizeof(Real);
	
        for(int i=0;i<d_len;i++){

            v_ptr[i]=drand48();
        }
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
#pragma omp parallel for
        for(int ss=0;ss<half._grid->oSites();ss++){
            half._odata[ss] = full._odata[ss*2+cb];
        }
        half.checkerboard = cb;
    }
    friend void setCheckerboard(Lattice<vobj> &full,const Lattice<vobj> &half){
        int cb = half.checkerboard;
#pragma omp parallel for
        for(int ss=0;ss<half._grid->oSites();ss++){
            full._odata[ss*2+cb]=half._odata[ss];
        }
    }
}; // class Lattice

    /* Need to implement the multiplication return type matching S S -> S, S M -> M, M S -> M through
     all nested possibilities.
     template<template<class> class lhs,template<class> class rhs>
     class MultTypeSelector {
     template<typename vtype> using ltype = lhs
     typedef lhs type;
     };
     */
    
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

	  const char * ptr =(const char*)&lhs._odata[ss];
#ifdef PREFETCH
          v_prefetch0(sizeof(obj2), ptr);
#endif

	  for(int i=0;i<sizeof(obj2);i+=64){
	    _mm_prefetch(ptr+i+4096,_MM_HINT_T1);
	    _mm_prefetch(ptr+i+256,_MM_HINT_T0);
	  }

	  ptr =(const char*)&rhs._odata[ss];
#ifdef PREFETCH
          v_prefetch0(sizeof(obj3), ptr);
#endif

	  for(int i=0;i<sizeof(obj3);i+=64){
	    _mm_prefetch(ptr+i+4096,_MM_HINT_T1);
	    _mm_prefetch(ptr+i+256,_MM_HINT_T0);
	  }


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
	std::cerr <<"Oscalar * Lattice calling mult"<<std::endl;
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
	std::cerr <<"Lattice * Oscalar calling mult"<<std::endl;
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
