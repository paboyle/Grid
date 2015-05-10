#ifndef GRID_LATTICE_OVERLOAD_H
#define GRID_LATTICE_OVERLOAD_H

namespace Grid {

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // unary negation
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class vobj>
  inline Lattice<vobj> operator -(const Lattice<vobj> &r)
  {
    Lattice<vobj> ret(r._grid);
#pragma omp parallel for
    for(int ss=0;ss<r._grid->oSites();ss++){
      vstream(ret._odata[ss], -r._odata[ss]);
    }
    return ret;
  } 
  /////////////////////////////////////////////////////////////////////////////////////
  // Lattice BinOp Lattice,
  //NB mult performs conformable check. Do not reapply here for performance.
  /////////////////////////////////////////////////////////////////////////////////////
  template<class left,class right>
    inline auto operator * (const Lattice<left> &lhs,const Lattice<right> &rhs)-> Lattice<decltype(lhs._odata[0]*rhs._odata[0])>
  {
    Lattice<decltype(lhs._odata[0]*rhs._odata[0])> ret(rhs._grid);
    mult(ret,lhs,rhs);
    return ret;
  }
  template<class left,class right>
    inline auto operator + (const Lattice<left> &lhs,const Lattice<right> &rhs)-> Lattice<decltype(lhs._odata[0]+rhs._odata[0])>
  {
    Lattice<decltype(lhs._odata[0]+rhs._odata[0])> ret(rhs._grid);
    add(ret,lhs,rhs);
    return ret;
  }
  template<class left,class right>
    inline auto operator - (const Lattice<left> &lhs,const Lattice<right> &rhs)-> Lattice<decltype(lhs._odata[0]-rhs._odata[0])>
  {
    Lattice<decltype(lhs._odata[0]-rhs._odata[0])> ret(rhs._grid);
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
      decltype(lhs*rhs._odata[0]) tmp=lhs*rhs._odata[ss]; 
      vstream(ret._odata[ss],tmp);
	   //      ret._odata[ss]=lhs*rhs._odata[ss];
    }
    return ret;
  }
  template<class left,class right>
    inline auto operator + (const left &lhs,const Lattice<right> &rhs) -> Lattice<decltype(lhs+rhs._odata[0])>
    {
      Lattice<decltype(lhs+rhs._odata[0])> ret(rhs._grid);
#pragma omp parallel for
      for(int ss=0;ss<rhs._grid->oSites(); ss++){
	decltype(lhs+rhs._odata[0]) tmp =lhs-rhs._odata[ss];  
	vstream(ret._odata[ss],tmp);
	//	ret._odata[ss]=lhs+rhs._odata[ss];
      }
        return ret;
    }
  template<class left,class right>
    inline auto operator - (const left &lhs,const Lattice<right> &rhs) -> Lattice<decltype(lhs-rhs._odata[0])>
  {
    Lattice<decltype(lhs-rhs._odata[0])> ret(rhs._grid);
#pragma omp parallel for
    for(int ss=0;ss<rhs._grid->oSites(); ss++){
      decltype(lhs-rhs._odata[0]) tmp=lhs-rhs._odata[ss];  
      vstream(ret._odata[ss],tmp);
      //      ret._odata[ss]=lhs-rhs._odata[ss];
    }
    return ret;
  }
    template<class left,class right>
      inline auto operator * (const Lattice<left> &lhs,const right &rhs) -> Lattice<decltype(lhs._odata[0]*rhs)>
    {
      Lattice<decltype(lhs._odata[0]*rhs)> ret(lhs._grid);
#pragma omp parallel for
      for(int ss=0;ss<lhs._grid->oSites(); ss++){
	decltype(lhs._odata[0]*rhs) tmp =lhs._odata[ss]*rhs;
	vstream(ret._odata[ss],tmp);
	//            ret._odata[ss]=lhs._odata[ss]*rhs;
      }
      return ret;
    }
    template<class left,class right>
      inline auto operator + (const Lattice<left> &lhs,const right &rhs) -> Lattice<decltype(lhs._odata[0]+rhs)>
    {
        Lattice<decltype(lhs._odata[0]+rhs)> ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<rhs._grid->oSites(); ss++){
	  decltype(lhs._odata[0]+rhs) tmp=lhs._odata[ss]+rhs; 
	  vstream(ret._odata[ss],tmp);
	  //	  ret._odata[ss]=lhs._odata[ss]+rhs;
        }
        return ret;
    }
    template<class left,class right>
      inline auto operator - (const Lattice<left> &lhs,const right &rhs) -> Lattice<decltype(lhs._odata[0]-rhs)>
    {
      Lattice<decltype(lhs._odata[0]-rhs)> ret(lhs._grid);
#pragma omp parallel for
      for(int ss=0;ss<rhs._grid->oSites(); ss++){
	  decltype(lhs._odata[0]-rhs) tmp=lhs._odata[ss]-rhs;
	  vstream(ret._odata[ss],tmp);
	  //	ret._odata[ss]=lhs._odata[ss]-rhs;
      }
      return ret;
    }


}
