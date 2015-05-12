#ifndef GRID_LATTICE_ARITH_H
#define GRID_LATTICE_ARITH_H

namespace Grid {


  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //  avoid copy back routines for mult, mac, sub, add
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class obj1,class obj2,class obj3>
    void mult(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
    conformable(lhs,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<lhs._grid->oSites();ss++){
      obj1 tmp;
      mult(&tmp,&lhs._odata[ss],&rhs._odata[ss]);
      vstream(ret._odata[ss],tmp);
    }
  }
  
  template<class obj1,class obj2,class obj3>
    void mac(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
    conformable(lhs,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<lhs._grid->oSites();ss++){
      obj1 tmp;
      mac(&tmp,&lhs._odata[ss],&rhs._odata[ss]);
      vstream(ret._odata[ss],tmp);
    }
  }
  
  template<class obj1,class obj2,class obj3>
    void sub(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
    conformable(lhs,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<lhs._grid->oSites();ss++){
      obj1 tmp;
      sub(&tmp,&lhs._odata[ss],&rhs._odata[ss]);
      vstream(ret._odata[ss],tmp);
    }
  }
  template<class obj1,class obj2,class obj3>
    void add(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
    conformable(lhs,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<lhs._grid->oSites();ss++){
      obj1 tmp;
      add(&tmp,&lhs._odata[ss],&rhs._odata[ss]);
      vstream(ret._odata[ss],tmp);
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //  avoid copy back routines for mult, mac, sub, add
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class obj1,class obj2,class obj3>
    void mult(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
    conformable(lhs,ret);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<lhs._grid->oSites();ss++){
      obj1 tmp;
      mult(&tmp,&lhs._odata[ss],&rhs);
      vstream(ret._odata[ss],tmp);
    }
  }
  
  template<class obj1,class obj2,class obj3>
    void mac(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
    conformable(lhs,ret);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<lhs._grid->oSites();ss++){
      obj1 tmp;
      mac(&tmp,&lhs._odata[ss],&rhs);
      vstream(ret._odata[ss],tmp);
    }
  }
  
  template<class obj1,class obj2,class obj3>
    void sub(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
    conformable(lhs,ret);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<lhs._grid->oSites();ss++){
      obj1 tmp;
      sub(&tmp,&lhs._odata[ss],&rhs);
      vstream(ret._odata[ss],tmp);
    }
  }
  template<class obj1,class obj2,class obj3>
    void add(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
    conformable(lhs,ret);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<lhs._grid->oSites();ss++){
      obj1 tmp;
      add(&tmp,&lhs._odata[ss],&rhs);
      vstream(ret._odata[ss],tmp);
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //  avoid copy back routines for mult, mac, sub, add
  //////////////////////////////////////////////////////////////////////////////////////////////////////
    template<class obj1,class obj2,class obj3>
    void mult(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
    conformable(ret,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<rhs._grid->oSites();ss++){
      obj1 tmp;
      mult(&tmp,&lhs,&rhs._odata[ss]);
      vstream(ret._odata[ss],tmp);
    }
  }
  
  template<class obj1,class obj2,class obj3>
    void mac(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
    conformable(ret,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<rhs._grid->oSites();ss++){
      obj1 tmp;
      mac(&tmp,&lhs,&rhs._odata[ss]);
      vstream(ret._odata[ss],tmp);
    }
  }
  
  template<class obj1,class obj2,class obj3>
    void sub(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
    conformable(ret,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<rhs._grid->oSites();ss++){
      obj1 tmp;
      sub(&tmp,&lhs,&rhs._odata[ss]);
      vstream(ret._odata[ss],tmp);
    }
  }
  template<class obj1,class obj2,class obj3>
    void add(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
    conformable(ret,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<rhs._grid->oSites();ss++){
      obj1 tmp;
      add(&tmp,&lhs,&rhs._odata[ss]);
      vstream(ret._odata[ss],tmp);
    }
  }
  
  template<class sobj,class vobj>
  inline void axpy(Lattice<vobj> &ret,sobj a,const Lattice<vobj> &lhs,const Lattice<vobj> &rhs){
    conformable(lhs,rhs);
#pragma omp parallel for
    for(int ss=0;ss<lhs._grid->oSites();ss++){
      vobj tmp = a*lhs._odata[ss];
      vstream(ret._odata[ss],tmp+rhs._odata[ss]);
    }
  }

}
#endif
