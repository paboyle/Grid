    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/spin/Dirac.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef GRID_QCD_DIRAC_H
#define GRID_QCD_DIRAC_H
namespace Grid{

namespace QCD {


  class Gamma {

  public:

    const int Ns=4;
    
    enum GammaMatrix {
      Identity,                        
      GammaX, 
      GammaY, 
      GammaZ, 
      GammaT,  
      Gamma5,
      MinusIdentity,                        
      MinusGammaX, 
      MinusGammaY, 
      MinusGammaZ, 
      MinusGammaT,  
      MinusGamma5
      //      GammaXGamma5, // Rest are composite (willing to take hit for two calls sequentially)
      //      GammaYGamma5, // as they are less commonly used.
      //      GammaZGamma5, 
      //      GammaTGamma5,  
      //      SigmaXY,
      //      SigmaXZ,
      //      SigmaYZ,
      //      SigmaXT,
      //      SigmaYT,
      //      SigmaZT,
      //      MinusGammaXGamma5, easiest to form by composition
      //      MinusGammaYGamma5, as performance is not critical for these
      //      MinusGammaZGamma5, 
      //      MinusGammaTGamma5,  
      //      MinusSigmaXY,
      //      MinusSigmaXZ,
      //      MinusSigmaYZ,
      //      MinusSigmaXT,
      //      MinusSigmaYT,
      //      MinusSigmaZT
    };
    
    static GammaMatrix GammaMatrices[];
    static const char *GammaMatrixNames[];

    Gamma (GammaMatrix g) { _g=g; }

    GammaMatrix _g;

  };
  
    // Make gamma products (Chroma convention)
    SpinMatrix makeGammaProd(const unsigned int i);
    
    /* Gx
     *  0 0  0  i    
     *  0 0  i  0    
     *  0 -i 0  0
     * -i 0  0  0
     */
  // right multiplication makes sense for matrix args, not for vector since there is 
  // no concept of row versus columnar indices
    template<class vtype> inline void rmultMinusGammaX(iMatrix<vtype,Ns> &ret,const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(i,0) = timesI(rhs(i,3));
	ret(i,1) = timesI(rhs(i,2));
	ret(i,2) = timesMinusI(rhs(i,1));
	ret(i,3) = timesMinusI(rhs(i,0));
      }
    };
    template<class vtype> inline void rmultGammaX(iMatrix<vtype,Ns> &ret, const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(i,0) = timesMinusI(rhs(i,3));
	ret(i,1) = timesMinusI(rhs(i,2));
	ret(i,2) = timesI(rhs(i,1));
	ret(i,3) = timesI(rhs(i,0));
      }
    };
    template<class vtype> inline void multGammaX(iMatrix<vtype,Ns> &ret, const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(0,i) = timesI(rhs(3,i));
	ret(1,i) = timesI(rhs(2,i));
	ret(2,i) = timesMinusI(rhs(1,i));
	ret(3,i) = timesMinusI(rhs(0,i));
      }
    };
    template<class vtype> inline void multMinusGammaX(iMatrix<vtype,Ns> &ret, const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(0,i) = timesMinusI(rhs(3,i));
	ret(1,i) = timesMinusI(rhs(2,i));
	ret(2,i) = timesI(rhs(1,i));
	ret(3,i) = timesI(rhs(0,i));
      }
    };

    template<class vtype> inline void multGammaX(iVector<vtype,Ns> &ret, const iVector<vtype,Ns> &rhs){
      ret._internal[0] = timesI(rhs._internal[3]);
      ret._internal[1] = timesI(rhs._internal[2]);
      ret._internal[2] = timesMinusI(rhs._internal[1]);
      ret._internal[3] = timesMinusI(rhs._internal[0]);
    };
    template<class vtype> inline void multMinusGammaX(iVector<vtype,Ns> &ret, const iVector<vtype,Ns> &rhs){
	ret(0) = timesMinusI(rhs(3));
	ret(1) = timesMinusI(rhs(2));
	ret(2) = timesI(rhs(1));
	ret(3) = timesI(rhs(0));
    };


    /*Gy
     *  0 0  0  -1  [0] -+ [3]
     *  0 0  1  0   [1] +- [2]
     *  0 1  0  0
     * -1 0  0  0
     */
    template<class vtype> inline void rmultGammaY(iMatrix<vtype,Ns> &ret, const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(i,0) = -rhs(i,3);
	ret(i,1) =  rhs(i,2);
	ret(i,2) =  rhs(i,1);
	ret(i,3) = -rhs(i,0);
      }
    };
    template<class vtype> inline void rmultMinusGammaY(iMatrix<vtype,Ns> &ret, const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(i,0) =  rhs(i,3);
	ret(i,1) = -rhs(i,2);
	ret(i,2) = -rhs(i,1);
	ret(i,3) =  rhs(i,0);
      }
    };
    template<class vtype> inline void multGammaY(iMatrix<vtype,Ns> &ret, const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(0,i) = -rhs(3,i);
	ret(1,i) =  rhs(2,i);
	ret(2,i) =  rhs(1,i);
	ret(3,i) = -rhs(0,i);
      }
    };
    template<class vtype> inline void multMinusGammaY(iMatrix<vtype,Ns> &ret, const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(0,i) =  rhs(3,i);
	ret(1,i) = -rhs(2,i);
	ret(2,i) = -rhs(1,i);
	ret(3,i) =  rhs(0,i);
      }
    };
    template<class vtype> inline void multGammaY(iVector<vtype,Ns> &ret, const iVector<vtype,Ns> &rhs){
      ret(0) = -rhs(3);
      ret(1) =  rhs(2);
      ret(2) =  rhs(1);
      ret(3) = -rhs(0);
    };
    template<class vtype> inline void multMinusGammaY(iVector<vtype,Ns> &ret, const iVector<vtype,Ns> &rhs){
      ret(0) =  rhs(3);
      ret(1) = -rhs(2);
      ret(2) = -rhs(1);
      ret(3) =  rhs(0);
    };
    /*Gz
     *  0 0  i  0   [0]+-i[2]
     *  0 0  0 -i   [1]-+i[3]
     * -i 0  0  0
     *  0 i  0  0
     */
    template<class vtype> inline void rmultGammaZ(iMatrix<vtype,Ns> &ret, const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(i,0) = timesMinusI(rhs(i,2));
	ret(i,1) =      timesI(rhs(i,3));
	ret(i,2) =      timesI(rhs(i,0));
	ret(i,3) = timesMinusI(rhs(i,1));
      }
    };
    template<class vtype> inline void rmultMinusGammaZ(iMatrix<vtype,Ns> &ret, const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(i,0) =      timesI(rhs(i,2));
	ret(i,1) = timesMinusI(rhs(i,3));
	ret(i,2) = timesMinusI(rhs(i,0));
	ret(i,3) =      timesI(rhs(i,1));
      }
    };
    template<class vtype> inline void multGammaZ(iMatrix<vtype,Ns> &ret, const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(0,i) = timesI(rhs(2,i));
	ret(1,i) =timesMinusI(rhs(3,i));
	ret(2,i) =timesMinusI(rhs(0,i));
	ret(3,i) = timesI(rhs(1,i));
      }
    };
    template<class vtype> inline void multMinusGammaZ(iMatrix<vtype,Ns> &ret, const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(0,i) = timesMinusI(rhs(2,i));
	ret(1,i) = timesI(rhs(3,i));
	ret(2,i) = timesI(rhs(0,i));
	ret(3,i) = timesMinusI(rhs(1,i));
      }
    };
    template<class vtype> inline void multGammaZ(iVector<vtype,Ns> &ret, const iVector<vtype,Ns> &rhs){
      ret(0) = timesI(rhs(2));
      ret(1) =timesMinusI(rhs(3));
      ret(2) =timesMinusI(rhs(0));
      ret(3) = timesI(rhs(1));
    };
    template<class vtype> inline void multMinusGammaZ(iVector<vtype,Ns> &ret, const iVector<vtype,Ns> &rhs){
      ret(0) = timesMinusI(rhs(2));
      ret(1) = timesI(rhs(3));
      ret(2) = timesI(rhs(0));
      ret(3) = timesMinusI(rhs(1));
    };
    /*Gt
     *  0 0  1  0 [0]+-[2]
     *  0 0  0  1 [1]+-[3]
     *  1 0  0  0
     *  0 1  0  0
     */
    template<class vtype> inline void rmultGammaT(iMatrix<vtype,Ns> &ret, const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(i,0) = rhs(i,2);
	ret(i,1) = rhs(i,3);
	ret(i,2) = rhs(i,0);
	ret(i,3) = rhs(i,1);
      }
    };
    template<class vtype> inline void rmultMinusGammaT(iMatrix<vtype,Ns> &ret, const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(i,0) =- rhs(i,2);
	ret(i,1) =- rhs(i,3);
	ret(i,2) =- rhs(i,0);
	ret(i,3) =- rhs(i,1);
      }
    };
    template<class vtype> inline void multGammaT(iMatrix<vtype,Ns> &ret, const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(0,i) = rhs(2,i);
	ret(1,i) = rhs(3,i);
	ret(2,i) = rhs(0,i);
	ret(3,i) = rhs(1,i);
      }
    };
    template<class vtype> inline void multMinusGammaT(iMatrix<vtype,Ns> &ret, const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(0,i) =-rhs(2,i);
	ret(1,i) =-rhs(3,i);
	ret(2,i) =-rhs(0,i);
	ret(3,i) =-rhs(1,i);
      }
    };
    template<class vtype> inline void multGammaT(iVector<vtype,Ns> &ret, const iVector<vtype,Ns> &rhs){
      ret(0) = rhs(2);
      ret(1) = rhs(3);
      ret(2) = rhs(0);
      ret(3) = rhs(1);
    };
    template<class vtype> inline void multMinusGammaT(iVector<vtype,Ns> &ret, const iVector<vtype,Ns> &rhs){
      ret(0) =-rhs(2);
      ret(1) =-rhs(3);
      ret(2) =-rhs(0);
      ret(3) =-rhs(1);
    };
    /*G5
     *  1 0  0  0 [0]+-[2]
     *  0 1  0  0 [1]+-[3]
     *  0 0 -1  0
     *  0 0  0 -1
     */
    template<class vtype> inline void rmultGamma5(iMatrix<vtype,Ns> &ret, const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(i,0) = rhs(i,0);
	ret(i,1) = rhs(i,1);
	ret(i,2) =-rhs(i,2);
	ret(i,3) =-rhs(i,3);
      }
    };
    template<class vtype> inline void rmultMinusGamma5(iMatrix<vtype,Ns> &ret, const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(i,0) =-rhs(i,0);
	ret(i,1) =-rhs(i,1);
	ret(i,2) = rhs(i,2);
	ret(i,3) = rhs(i,3);
      }
    };

    template<class vtype> inline void multGamma5(iMatrix<vtype,Ns> &ret, const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(0,i) = rhs(0,i);
	ret(1,i) = rhs(1,i);
	ret(2,i) =-rhs(2,i);
	ret(3,i) =-rhs(3,i);
      }
    };
    template<class vtype> inline void multMinusGamma5(iMatrix<vtype,Ns> &ret, const iMatrix<vtype,Ns> &rhs){
      for(int i=0;i<Ns;i++){
	ret(0,i) =-rhs(0,i);
	ret(1,i) =-rhs(1,i);
	ret(2,i) = rhs(2,i);
	ret(3,i) = rhs(3,i);
      }
    };

    template<class vtype> inline void multGamma5(iVector<vtype,Ns> &ret, const iVector<vtype,Ns> &rhs){
      ret(0) = rhs(0);
      ret(1) = rhs(1);
      ret(2) =-rhs(2);
      ret(3) =-rhs(3);
    };
    template<class vtype> inline void multMinusGamma5(iVector<vtype,Ns> &ret, const iVector<vtype,Ns> &rhs){
      ret(0) =-rhs(0);
      ret(1) =-rhs(1);
      ret(2) = rhs(2);
      ret(3) = rhs(3);
    };



#ifdef GRID_WARN_SUBOPTIMAL
#warning "Optimisation alert switch over to multGammaX early "
#endif     

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Operator * : first case this is not a spin index, so recurse
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    
    // FIXME
    //
    // Optimisation; switch over to a "multGammaX(ret._internal,arg._internal)" style early and
    // note that doing so from the lattice operator will avoid copy back and case switch overhead, as
    // was done for the tensor math operator to remove operator * notation early
    //

    //left multiply
    template<class vtype> inline auto operator * ( const Gamma &G,const iScalar<vtype> &arg) ->
      typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type 

	{
	  iScalar<vtype> ret;
	  ret._internal=G*arg._internal;
	  return ret;
	}
    template<class vtype,int N> inline auto operator * ( const Gamma &G,const iVector<vtype,N> &arg) ->
      typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type 
	{
	  iVector<vtype,N> ret;
	  for(int i=0;i<N;i++){
	    ret._internal[i]=G*arg._internal[i];
	  }
	  return ret;
	}
    template<class vtype,int N> inline auto operator * ( const Gamma &G,const iMatrix<vtype,N> &arg) ->
      typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type 
	{
	  iMatrix<vtype,N> ret;
	  for(int i=0;i<N;i++){
	  for(int j=0;j<N;j++){
	    ret._internal[i][j]=G*arg._internal[i][j];
	  }}
	  return ret;
	}


    //right multiply
    template<class vtype> inline auto operator * (const iScalar<vtype> &arg, const Gamma &G) ->
      typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type 

	{
	  iScalar<vtype> ret;
	  ret._internal=arg._internal*G;
	  return ret;
	}
    template<class vtype,int N> inline auto operator * (const iVector<vtype,N> &arg, const Gamma &G) ->
      typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type 
	{
	  iVector<vtype,N> ret;
	  for(int i=0;i<N;i++){
	    ret._internal=arg._internal[i]*G;
	  }
	  return ret;
	}
    template<class vtype,int N> inline auto operator * (const iMatrix<vtype,N> &arg, const Gamma &G) ->
      typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type 
	{
	  iMatrix<vtype,N> ret;
	  for(int i=0;i<N;i++){
	  for(int j=0;j<N;j++){
	    ret._internal[i][j]=arg._internal[i][j]*G;
	  }}
	  return ret;
	}

    ////////////////////////////////////////////////////////
    // When we hit the spin index this matches and we stop
    ////////////////////////////////////////////////////////
    template<class vtype> inline auto operator * ( const Gamma &G,const iMatrix<vtype,Ns> &arg) ->
      typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,Ns>,SpinorIndex>::value,iMatrix<vtype,Ns> >::type 
      {
	iMatrix<vtype,Ns> ret;
	switch (G._g) {
	case Gamma::Identity:
	  ret = arg;
	  break;
	case Gamma::MinusIdentity:
	  ret = -arg;
	  break;
	case Gamma::GammaX:
	  multGammaX(ret,arg);
	  break;
	case Gamma::MinusGammaX:
	  multMinusGammaX(ret,arg);
	  break;
	case Gamma::GammaY:
	  multGammaY(ret,arg);
	  break;
	case Gamma::MinusGammaY:
	  multMinusGammaY(ret,arg);
	  break;
	case Gamma::GammaZ:
	  multGammaZ(ret,arg);
	  break;
	case Gamma::MinusGammaZ:
	  multMinusGammaZ(ret,arg);
	  break;
	case Gamma::GammaT:
	  multGammaT(ret,arg);
	  break;
	case Gamma::MinusGammaT:
	  multMinusGammaT(ret,arg);
	  break;
	case Gamma::Gamma5:
	  multGamma5(ret,arg);
	  break;
	case Gamma::MinusGamma5:
	  multMinusGamma5(ret,arg);
	  break;
	default:
	  assert(0);
	  break;
	}
	return ret;
      }
    // Could have used type trait for Matrix/vector and then an enable if to share code
    template<class vtype> inline auto operator * ( const Gamma &G,const iVector<vtype,Ns> &arg) ->
      typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type 
      {
	iVector<vtype,Ns> ret;
	switch (G._g) {
	case Gamma::Identity:
	  ret = arg;
	  break;
	case Gamma::MinusIdentity:
	  ret = -arg;
	  break;
	case Gamma::GammaX:
	  multGammaX(ret,arg);
	  break;
	case Gamma::MinusGammaX:
	  multMinusGammaX(ret,arg);
	  break;
	case Gamma::GammaY:
	  multGammaY(ret,arg);
	  break;
	case Gamma::MinusGammaY:
	  multMinusGammaY(ret,arg);
	  break;
	case Gamma::GammaZ:
	  multGammaZ(ret,arg);
	  break;
	case Gamma::MinusGammaZ:
	  multMinusGammaZ(ret,arg);
	  break;
	case Gamma::GammaT:
	  multGammaT(ret,arg);
	  break;
	case Gamma::MinusGammaT:
	  multMinusGammaT(ret,arg);
	  break;
	case Gamma::Gamma5:
	  multGamma5(ret,arg);
	  break;
	case Gamma::MinusGamma5:
	  multMinusGamma5(ret,arg);
	  break;
	default:
	  assert(0);
	  break;
	}
	return ret;
      }

    template<class vtype> inline auto operator * (const iMatrix<vtype,Ns> &arg, const Gamma &G) ->
      typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,Ns>,SpinorIndex>::value,iMatrix<vtype,Ns> >::type 
      {
	iMatrix<vtype,Ns> ret;
	switch (G._g) {
	case Gamma::Identity:
	  ret = arg;
	  break;
	case Gamma::MinusIdentity:
	  ret = -arg;
	  break;
	case Gamma::GammaX:
	  rmultGammaX(ret,arg);
	  break;
	case Gamma::MinusGammaX:
	  rmultMinusGammaX(ret,arg);
	  break;
	case Gamma::GammaY:
	  rmultGammaY(ret,arg);
	  break;
	case Gamma::MinusGammaY:
	  rmultMinusGammaY(ret,arg);
	  break;
	case Gamma::GammaZ:
	  rmultGammaZ(ret,arg);
	  break;
	case Gamma::MinusGammaZ:
	  rmultMinusGammaZ(ret,arg);
	  break;
	case Gamma::GammaT:
	  rmultGammaT(ret,arg);
	  break;
	case Gamma::MinusGammaT:
	  rmultMinusGammaT(ret,arg);
	  break;
	case Gamma::Gamma5:
	  rmultGamma5(ret,arg);
	  break;
	case Gamma::MinusGamma5:
	  rmultMinusGamma5(ret,arg);
	  break;
	default:
	  assert(0);
	  break;
	}
	return ret;
      }


    /* Output from test
./Grid_gamma 
Identity((1,0),(0,0),(0,0),(0,0))
        ((0,0),(1,0),(0,0),(0,0))
        ((0,0),(0,0),(1,0),(0,0))
        ((0,0),(0,0),(0,0),(1,0))   OK

GammaX  ((0,0),(0,0),(0,0),(0,1))
        ((0,0),(0,0),(0,1),(0,0))
        ((0,0),(0,-1),(0,0),(0,0))
        ((0,-1),(0,0),(0,0),(0,0))  OK
    * Gx
    *  0 0  0  i    
    *  0 0  i  0    
    *  0 -i 0  0
    * -i 0  0  0

GammaY  ((-0,-0),(-0,-0),(-0,-0),(-1,-0))
        ((0,0),(0,0),(1,0),(0,0))
        ((0,0),(1,0),(0,0),(0,0))          OK
        ((-1,-0),(-0,-0),(-0,-0),(-0,-0))
     *Gy
     *  0 0  0  -1  [0] -+ [3]
     *  0 0  1  0   [1] +- [2]
     *  0 1  0  0
     * -1 0  0  0

GammaZ  ((0,0),(0,0),(0,1),(0,0))
        ((0,0),(0,0),(0,0),(0,-1))
        ((0,-1),(0,0),(0,0),(0,0))
        ((0,0),(0,1),(0,0),(0,0))   OK
     *  0 0  i  0   [0]+-i[2]
     *  0 0  0 -i   [1]-+i[3]
     * -i 0  0  0
     *  0 i  0  0

GammaT  ((0,0),(0,0),(1,0),(0,0))
        ((0,0),(0,0),(0,0),(1,0))  OK
        ((1,0),(0,0),(0,0),(0,0))
        ((0,0),(1,0),(0,0),(0,0))
     *  0 0  1  0 [0]+-[2]
     *  0 0  0  1 [1]+-[3]
     *  1 0  0  0
     *  0 1  0  0

Gamma5  ((1,0),(0,0),(0,0),(0,0))
        ((0,0),(1,0),(0,0),(0,0))
        ((-0,-0),(-0,-0),(-1,-0),(-0,-0))
        ((-0,-0),(-0,-0),(-0,-0),(-1,-0))
     *  1 0  0  0 [0]+-[2]
     *  0 1  0  0 [1]+-[3]  OK
     *  0 0 -1  0
     *  0 0  0 -1
     */

}   //namespace QCD
} // Grid
#endif
