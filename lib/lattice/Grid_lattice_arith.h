#ifndef GRID_LATTICE_ARITH_H
#define GRID_LATTICE_ARITH_H

namespace Grid {


  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //  avoid copy back routines for mult, mac, sub, add
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class obj1,class obj2,class obj3> strong_inline
    void mult(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
    conformable(lhs,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<lhs._grid->oSites();ss++){
#ifdef STREAMING_STORES
      obj1 tmp;
      mult(&tmp,&lhs._odata[ss],&rhs._odata[ss]);
      vstream(ret._odata[ss],tmp);
#else
      mult(&ret._odata[ss],&lhs._odata[ss],&rhs._odata[ss]);
#endif
    }
  }
  
  template<class obj1,class obj2,class obj3> strong_inline
    void mac(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
    conformable(lhs,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<lhs._grid->oSites();ss++){
#ifdef STREAMING_STORES
      obj1 tmp;
      mac(&tmp,&lhs._odata[ss],&rhs._odata[ss]);
      vstream(ret._odata[ss],tmp);
#else
      mac(&ret._odata[ss],&lhs._odata[ss],&rhs._odata[ss]);
#endif
    }
  }
  
  template<class obj1,class obj2,class obj3> strong_inline
    void sub(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
    conformable(lhs,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<lhs._grid->oSites();ss++){
#ifdef STREAMING_STORES
      obj1 tmp;
      sub(&tmp,&lhs._odata[ss],&rhs._odata[ss]);
      vstream(ret._odata[ss],tmp);
#else
      sub(&ret._odata[ss],&lhs._odata[ss],&rhs._odata[ss]);
#endif
    }
  }
  template<class obj1,class obj2,class obj3> strong_inline
    void add(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
    conformable(lhs,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<lhs._grid->oSites();ss++){
#ifdef STREAMING_STORES
      obj1 tmp;
      add(&tmp,&lhs._odata[ss],&rhs._odata[ss]);
      vstream(ret._odata[ss],tmp);
#else
      add(&ret._odata[ss],&lhs._odata[ss],&rhs._odata[ss]);
#endif
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //  avoid copy back routines for mult, mac, sub, add
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class obj1,class obj2,class obj3> strong_inline
    void mult(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
    conformable(lhs,ret);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<lhs._grid->oSites();ss++){
      obj1 tmp;
      mult(&tmp,&lhs._odata[ss],&rhs);
      vstream(ret._odata[ss],tmp);
    }
  }
  
  template<class obj1,class obj2,class obj3> strong_inline
    void mac(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
    conformable(lhs,ret);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<lhs._grid->oSites();ss++){
      obj1 tmp;
      mac(&tmp,&lhs._odata[ss],&rhs);
      vstream(ret._odata[ss],tmp);
    }
  }
  
  template<class obj1,class obj2,class obj3> strong_inline
    void sub(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
    conformable(lhs,ret);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<lhs._grid->oSites();ss++){
#ifdef STREAMING_STORES
      obj1 tmp;
      sub(&tmp,&lhs._odata[ss],&rhs);
      vstream(ret._odata[ss],tmp);
#else 
      sub(&ret._odata[ss],&lhs._odata[ss],&rhs);
#endif
    }
  }
  template<class obj1,class obj2,class obj3> strong_inline
    void add(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
    conformable(lhs,ret);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<lhs._grid->oSites();ss++){
#ifdef STREAMING_STORES
      obj1 tmp;
      add(&tmp,&lhs._odata[ss],&rhs);
      vstream(ret._odata[ss],tmp);
#else 
      add(&ret._odata[ss],&lhs._odata[ss],&rhs);
#endif
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //  avoid copy back routines for mult, mac, sub, add
  //////////////////////////////////////////////////////////////////////////////////////////////////////
    template<class obj1,class obj2,class obj3> strong_inline
    void mult(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
    conformable(ret,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<rhs._grid->oSites();ss++){
#ifdef STREAMING_STORES
      obj1 tmp;
      mult(&tmp,&lhs,&rhs._odata[ss]);
      vstream(ret._odata[ss],tmp);
#else 
      mult(&ret._odata[ss],&lhs,&rhs._odata[ss]);
#endif
    }
  }
  
  template<class obj1,class obj2,class obj3> strong_inline
    void mac(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
    conformable(ret,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<rhs._grid->oSites();ss++){
#ifdef STREAMING_STORES
      obj1 tmp;
      mac(&tmp,&lhs,&rhs._odata[ss]);
      vstream(ret._odata[ss],tmp);
#else 
      mac(&ret._odata[ss],&lhs,&rhs._odata[ss]);
#endif
    }
  }
  
  template<class obj1,class obj2,class obj3> strong_inline
    void sub(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
    conformable(ret,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<rhs._grid->oSites();ss++){
#ifdef STREAMING_STORES
      obj1 tmp;
      sub(&tmp,&lhs,&rhs._odata[ss]);
      vstream(ret._odata[ss],tmp);
#else 
      sub(&ret._odata[ss],&lhs,&rhs._odata[ss]);
#endif
    }
  }
  template<class obj1,class obj2,class obj3> strong_inline
    void add(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
    conformable(ret,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<rhs._grid->oSites();ss++){
#ifdef STREAMING_STORES
      obj1 tmp;
      add(&tmp,&lhs,&rhs._odata[ss]);
      vstream(ret._odata[ss],tmp);
#else 
      add(&ret._odata[ss],&lhs,&rhs._odata[ss]);
#endif
    }
  }
  
  template<class sobj,class vobj> strong_inline
  void axpy(Lattice<vobj> &ret,sobj a,const Lattice<vobj> &x,const Lattice<vobj> &y){
    conformable(x,y);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<x._grid->oSites();ss++){
#ifdef STREAMING_STORES
      vobj tmp = a*x._odata[ss]+y._odata[ss];
      vstream(ret._odata[ss],tmp);
#else
      ret._odata[ss]=a*x._odata[ss]+y._odata[ss];
#endif
    }
  }
  template<class sobj,class vobj> strong_inline
  void axpby(Lattice<vobj> &ret,sobj a,sobj b,const Lattice<vobj> &x,const Lattice<vobj> &y){
    conformable(x,y);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<x._grid->oSites();ss++){
#ifdef STREAMING_STORES
      vobj tmp = a*x._odata[ss]+b*y._odata[ss];
      vstream(ret._odata[ss],tmp);
#else
      ret._odata[ss]=a*x._odata[ss]+b*y._odata[ss];
#endif
    }
  }

  template<class sobj,class vobj> strong_inline
  RealD axpy_norm(Lattice<vobj> &ret,sobj a,const Lattice<vobj> &x,const Lattice<vobj> &y){
    conformable(x,y);
    axpy(ret,a,x,y);
    return norm2(ret);
  }
  template<class sobj,class vobj> strong_inline
  RealD axpby_norm(Lattice<vobj> &ret,sobj a,sobj b,const Lattice<vobj> &x,const Lattice<vobj> &y){
    conformable(x,y);
    axpby(ret,a,b,x,y);
    return norm2(ret); // FIXME implement parallel norm in ss loop
  }

}
#endif
