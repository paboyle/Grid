#ifndef GRID_LATTICE_ARITH_H
#define GRID_LATTICE_ARITH_H

namespace Grid {

  template<class vobj>
  inline Lattice<vobj> operator -(const Lattice<vobj> &r)
  {
    Lattice<vobj> ret(r._grid);
#pragma omp parallel for
    for(int ss=0;ss<r._grid->oSites();ss++){
      ret._odata[ss]= -r._odata[ss];
    }
    return ret;
  }
  
  template<class vobj>
  inline void axpy(Lattice<vobj> &ret,double a,const Lattice<vobj> &lhs,const Lattice<vobj> &rhs){
    conformable(lhs,rhs);
#pragma omp parallel for
    for(int ss=0;ss<lhs._grid->oSites();ss++){
      axpy(&ret._odata[ss],a,&lhs._odata[ss],&rhs._odata[ss]);
    }
  }
  template<class vobj>
  inline void axpy(Lattice<vobj> &ret,std::complex<double> a,const Lattice<vobj> &lhs,const Lattice<vobj> &rhs){
    conformable(lhs,rhs);
#pragma omp parallel for
    for(int ss=0;ss<lhs._grid->oSites();ss++){
      axpy(&ret._odata[ss],a,&lhs._odata[ss],&rhs._odata[ss]);
    }
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
}
#endif
