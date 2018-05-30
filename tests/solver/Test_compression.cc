    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_cg_prec.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#include<bitset>
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

//Preconditioning:   M psi = chi
//               =  M P^-1 P psi = chi
//               =  M P^-1 psi' = chi

//Solve for psi' using M P^-1 as operator, then apply  P^-1 psi' = psi 

//Inexact preconditioned CG requires slight modification because we want to avoid computing P^-1 exactly


/////////////////////////////////////////////////////////////
// Base classes for iterative processes based on operators
// single input vec, single output vec.
/////////////////////////////////////////////////////////////

//The compressor
#if 0


//Basic copy of WilsonCompressor for demonstration
template<class _Hspinor,class _Spinor, class projector>
class WilsonTestCompressorTemplate
{
 public:
  
  int mu,dag;  

  void Point(int p) { mu=p; };

  WilsonTestCompressorTemplate(int _dag=0){
    //printf("WilsonTestCompressorTemplate constructor\n");
    dag = _dag;
  }

  typedef _Spinor         SiteSpinor;
  typedef _Hspinor     SiteHalfSpinor;
  typedef _Hspinor SiteHalfCommSpinor;
  typedef typename SiteHalfSpinor::vector_type     vComplexIn;
  constexpr static int Nw=sizeof(SiteHalfSpinor)/sizeof(vComplexIn); //number of complex numbers in SiteHalfSpinor

  inline int CommDatumSize(void) {
    //printf("WilsonTestCompressorTemplate CommDatumSize\n");
    return sizeof(SiteHalfCommSpinor);
  }

  /*****************************************************/
  /* Compress includes precision change if mpi data is not same */
  /*****************************************************/
  inline void Compress(SiteHalfSpinor *buf,Integer o,const SiteSpinor &in) {
    //printf("WilsonTestCompressorTemplate Compress\n");
    projector::Proj(buf[o],in,mu,dag);
  }

  /*****************************************************/
  /* Exchange includes precision change if mpi data is not same */
  /*****************************************************/
  inline void Exchange(SiteHalfSpinor *mp,
                       SiteHalfSpinor *vp0,
                       SiteHalfSpinor *vp1,
		       Integer type,Integer o){
    //printf("WilsonTestCompressorTemplate Exchange\n");
    exchange(mp[2*o],mp[2*o+1],vp0[o],vp1[o],type);
  }

  /*****************************************************/
  /* Have a decompression step if mpi data is not same */
  /*****************************************************/
  inline void Decompress(SiteHalfSpinor *out,
			 SiteHalfSpinor *in, Integer o) {    
    //printf("WilsonTestCompressorTemplate Decompress\n");
    assert(0);
  }

  /*****************************************************/
  /* Compress Exchange                                 */
  /*****************************************************/
  inline void CompressExchange(SiteHalfSpinor *out0,
			       SiteHalfSpinor *out1,
			       const SiteSpinor *in,
			       Integer j,Integer k, Integer m,Integer type){
    //printf("WilsonTestCompressorTemplate CompressExchange\n");
    SiteHalfSpinor temp1, temp2,temp3,temp4;
    projector::Proj(temp1,in[k],mu,dag);
    projector::Proj(temp2,in[m],mu,dag);
    exchange(out0[j],out1[j],temp1,temp2,type);
  }

  /*****************************************************/
  /* Pass the info to the stencil */
  /*****************************************************/
  inline bool DecompressionStep(void) { return false; }

};

#elif 0


//Compressor that unpacks vectorized data to scalar
template<class _Hspinor,class _Spinor, class projector>
class WilsonTestCompressorTemplate
{
 public:
  
  int mu,dag;  

  void Point(int p) { mu=p; };

  WilsonTestCompressorTemplate(int _dag=0){
    dag = _dag;
  }

  typedef _Spinor         SiteSpinor;
  typedef _Hspinor     SiteHalfSpinor;
  typedef _Hspinor SiteHalfCommSpinor;
  typedef typename SiteHalfSpinor::vector_type     vComplexIn;
  constexpr static int Nw=sizeof(SiteHalfSpinor)/sizeof(vComplexIn); //number of complex numbers in SiteHalfSpinor

  typedef typename SiteHalfSpinor::scalar_object ScalarSiteHalfSpinor;

  constexpr static int Nsimd = vComplexIn::Nsimd();

  inline int CommDatumSize(void) {
    return Nsimd*sizeof(ScalarSiteHalfSpinor);
  }

  /*****************************************************/
  /* Compress includes precision change if mpi data is not same */
  /*****************************************************/
  inline void Compress(SiteHalfSpinor *buf,Integer o,const SiteSpinor &in) {
    SiteHalfSpinor hsp;
    projector::Proj(hsp,in,mu,dag);
    
    ScalarSiteHalfSpinor* to = (ScalarSiteHalfSpinor*)buf + o*Nsimd;
    
    std::vector<ScalarSiteHalfSpinor*> extract_args(Nsimd);
    for(int i=0;i<Nsimd;i++) extract_args[i] = to+i;
    extract1(hsp,extract_args,0);
  }

  /*****************************************************/
  /* Exchange includes precision change if mpi data is not same */
  /*****************************************************/
  inline void Exchange(SiteHalfSpinor *mp,
                       SiteHalfSpinor *vp0,
                       SiteHalfSpinor *vp1,
		       Integer type,Integer o){
    ScalarSiteHalfSpinor* vpp0 = (ScalarSiteHalfSpinor*)vp0 + o*Nsimd;
    ScalarSiteHalfSpinor* vpp1 = (ScalarSiteHalfSpinor*)vp1 + o*Nsimd;
    
    std::vector<ScalarSiteHalfSpinor*> merge_args0(Nsimd), merge_args1(Nsimd);
    for(int i=0;i<Nsimd;i++){
      merge_args0[i] = vpp0+i;
      merge_args1[i] = vpp1+i;
    }

    SiteHalfSpinor vt0,vt1;
    merge1(vt0,merge_args0,0);
    merge1(vt1,merge_args1,0);

    exchange(mp[2*o],mp[2*o+1],vt0,vt1,type);
  }

  /*****************************************************/
  /* Have a decompression step if mpi data is not same */
  /*****************************************************/
  inline void Decompress(SiteHalfSpinor *out,
			 SiteHalfSpinor *in, Integer o) {    
    ScalarSiteHalfSpinor* hin = (ScalarSiteHalfSpinor*)in + o*Nsimd;
    std::vector<ScalarSiteHalfSpinor*> merge_args(Nsimd);
    for(int i=0;i<Nsimd;i++) merge_args[i] = hin+i;
    merge1(out[o],merge_args,0);
  }

  /*****************************************************/
  /* Compress Exchange                                 */
  /*****************************************************/
  inline void CompressExchange(SiteHalfSpinor *out0,
			       SiteHalfSpinor *out1,
			       const SiteSpinor *in,
			       Integer j,Integer k, Integer m,Integer type){
    SiteHalfSpinor temp1, temp2,temp3,temp4;
    projector::Proj(temp1,in[k],mu,dag);
    projector::Proj(temp2,in[m],mu,dag);
    exchange(temp3,temp4,temp1,temp2,type);

    ScalarSiteHalfSpinor* hout0 = (ScalarSiteHalfSpinor*)out0 + j*Nsimd;
    ScalarSiteHalfSpinor* hout1 = (ScalarSiteHalfSpinor*)out1 + j*Nsimd;

    std::vector<ScalarSiteHalfSpinor*> extract_args0(Nsimd), extract_args1(Nsimd);
    for(int i=0;i<Nsimd;i++){
      extract_args0[i] = hout0+i;
      extract_args1[i] = hout1+i;
    }
    extract1(temp3,extract_args0,0);
    extract1(temp4,extract_args1,0);
  }

  /*****************************************************/
  /* Pass the info to the stencil */
  /*****************************************************/
  inline bool DecompressionStep(void) { return true; }

};

#else

//Access elements of std::complex
template<typename T>
inline T & cmplx_reim(std::complex<T> &c, const int reim){
  return reinterpret_cast<T(&)[2]>(c)[reim];
}

template<typename T>
inline const T & cmplx_reim(const std::complex<T> &c, const int reim){
  return reinterpret_cast<const T(&)[2]>(c)[reim];
}


//Pack and unpack float/double to fixed point representation of SZ bits
template<int SZ>
struct signedIntMap{};

template<>
struct signedIntMap<8>{ typedef int8_t type; };
template<>
struct signedIntMap<16>{ typedef int16_t type; };


template<typename T, int SZ>
inline typename signedIntMap<SZ>::type packN(T val){
  return typename signedIntMap<SZ>::type( (1<<(SZ-2) ) * val );
}
template<typename T, int SZ>
inline T unpackN(typename signedIntMap<SZ>::type val){
  return T(val)/(1<<(SZ-2));
}

template<typename T>
struct getHalfSpinorColors{
  //template <typename vtype> using iImplHalfSpinor        = iScalar<iVector<iVector<vtype, Dimension>, Nhs> >;
  enum { value = sizeof(typename T::element::element)/sizeof(typename T::element::element::element) };
};

//Compressor that compresses to a single magnitude and Nhs*Dimension fixed point integers of size packSize bits
template<class _Hspinor,class _Spinor, class projector, int packSize = 16>
class WilsonTestCompressorTemplate
{
 public:
  
  int mu,dag;  

  void Point(int p) { mu=p; };

  WilsonTestCompressorTemplate(int _dag=0){
    dag = _dag;
  }

  typedef _Spinor         SiteSpinor;
  typedef _Hspinor     SiteHalfSpinor;
  typedef _Hspinor SiteHalfCommSpinor;
  typedef typename SiteHalfSpinor::vector_type     vComplexIn;
  constexpr static int Nw=sizeof(SiteHalfSpinor)/sizeof(vComplexIn); //number of complex numbers in SiteHalfSpinor

  typedef typename SiteHalfSpinor::scalar_object ScalarSiteHalfSpinor;

  constexpr static int Nsimd = vComplexIn::Nsimd();
  constexpr static int Dimension = getHalfSpinorColors<SiteHalfSpinor>::value;
 
  typedef typename ScalarSiteHalfSpinor::scalar_type stype; //std::complex
  typedef typename stype::value_type srtype; //float/double

  //Pack and unpack *scalar* SiteHalfSpinor objects
  void packSpinor(void* tov, const ScalarSiteHalfSpinor &from){
    uint8_t* to = (uint8_t*)tov;
    typedef typename signedIntMap<packSize>::type packedType;

    srtype max = 0;
    for(int s=0;s<Nhs;s++)
      for(int c=0;c<Dimension;c++)
	for(int reim=0;reim<2;reim++)
	  if(fabs(cmplx_reim( from()(s)(c), reim )) > max )
	    max =  fabs(cmplx_reim( from()(s)(c), reim )) ;
  
    *( (srtype*)to ) = max; //copy the normalization to the buffer
    to += sizeof(srtype);
  
    packedType *top = (packedType*)to;
    packedType p;
    srtype q;
    for(int s=0;s<Nhs;s++)
      for(int c=0;c<Dimension;c++)
	for(int reim=0;reim<2;reim++){
	  q = cmplx_reim( from()(s)(c), reim );
	  if(max != 0.) q /= max;
	  *(top++) = packN<srtype,packSize>(q);
	}
  }

  void packSpinor(void* tov, const SiteHalfSpinor &from){
    uint8_t* to = (uint8_t*)tov;
    std::vector<ScalarSiteHalfSpinor> extracted(Nsimd);
    extract(from,extracted);

    static const int incr = sizeof(srtype) + Nhs*Dimension*2*sizeof(typename signedIntMap<packSize>::type);

    for(int i=0;i<Nsimd;i++){
      packSpinor((void*)to, extracted[i]);
      to += incr;
    }
  }


  void unpackSpinor(ScalarSiteHalfSpinor &to, void* fromv){
    uint8_t* from = (uint8_t*)fromv;
    typedef typename signedIntMap<packSize>::type packedType;

    srtype norm = *( (srtype*)from ); 
    from += sizeof(srtype);

    packedType *fromp = (packedType*)from;
    srtype q;
    for(int s=0;s<Nhs;s++)
      for(int c=0;c<Dimension;c++)
	for(int reim=0;reim<2;reim++){
	  q = unpackN<srtype,packSize>(*(fromp++) );
	  if(norm != 0.) q *= norm;
	  cmplx_reim( to()(s)(c), reim ) = q;
	}
  }

  void unpackSpinor(SiteHalfSpinor &to, void* fromv){
    uint8_t* from = (uint8_t*)fromv;
    std::vector<ScalarSiteHalfSpinor> unpacked(Nsimd);

    static const int incr = sizeof(srtype) + Nhs*Dimension*2*sizeof(typename signedIntMap<packSize>::type);

    for(int i=0;i<Nsimd;i++){
      unpackSpinor(unpacked[i],(void*)from);
      from += incr;
    }

    merge(to,unpacked);
  }

  inline int CommDatumSize(void) {
    return Nsimd*(  sizeof(srtype) + Nhs*Dimension*2*sizeof(typename signedIntMap<packSize>::type) );
  }

  /*****************************************************/
  /* Compress includes precision change if mpi data is not same */
  /*****************************************************/
  void Compress(SiteHalfSpinor *buf,Integer o,const SiteSpinor &in) {
    SiteHalfSpinor hsp;
    projector::Proj(hsp,in,mu,dag);
    
    uint8_t* to = (uint8_t*)buf + o*CommDatumSize();
    packSpinor(to, hsp);
  }

  /*****************************************************/
  /* Exchange includes precision change if mpi data is not same */
  /*****************************************************/
  void Exchange(SiteHalfSpinor *mp,
                       SiteHalfSpinor *vp0,
                       SiteHalfSpinor *vp1,
		       Integer type,Integer o){
    uint8_t* vpp0 = (uint8_t*)vp0 + o*CommDatumSize();
    uint8_t* vpp1 = (uint8_t*)vp1 + o*CommDatumSize();

    SiteHalfSpinor vt0, vt1;
    unpackSpinor(vt0, vpp0);
    unpackSpinor(vt1, vpp1);

    exchange(mp[2*o],mp[2*o+1],vt0,vt1,type);
  }

  /*****************************************************/
  /* Have a decompression step if mpi data is not same */
  /*****************************************************/
  void Decompress(SiteHalfSpinor *out,
			 SiteHalfSpinor *in, Integer o) {    
    uint8_t* hin = (uint8_t*)in + o*CommDatumSize();
    unpackSpinor(out[o],hin);
  }

  /*****************************************************/
  /* Compress Exchange                                 */
  /*****************************************************/
  void CompressExchange(SiteHalfSpinor *out0,
			       SiteHalfSpinor *out1,
			       const SiteSpinor *in,
			       Integer j,Integer k, Integer m,Integer type){
    SiteHalfSpinor temp1, temp2,temp3,temp4;
    projector::Proj(temp1,in[k],mu,dag);
    projector::Proj(temp2,in[m],mu,dag);
    exchange(temp3,temp4,temp1,temp2,type);

    uint8_t* hout0 = (uint8_t*)out0 + j*CommDatumSize();
    uint8_t* hout1 = (uint8_t*)out1 + j*CommDatumSize();
    packSpinor(hout0, temp3);
    packSpinor(hout1, temp4);
  }

  /*****************************************************/
  /* Pass the info to the stencil */
  /*****************************************************/
  inline bool DecompressionStep(void) { return true; }

};


#endif



template<typename HS,typename S, int packSize> using WilsonTestCompressor = WilsonTestCompressorTemplate<HS,S,WilsonProjector,packSize>;


template<class vobj,class cobj>
class WilsonStencilBasic : public CartesianStencil<vobj,cobj> {
public:
  double timer0;
  double timer1;
  double timer2;
  double timer3;
  double timer4;
  double timer5;
  double timer6;
  uint64_t callsi;
  void ZeroCountersi(void)
  {
    timer0=0;
    timer1=0;
    timer2=0;
    timer3=0;
    timer4=0;
    timer5=0;
    timer6=0;
    callsi=0;
  }
  void Reporti(int calls)
  {
    if ( timer0 ) std::cout << GridLogMessage << " timer0 (HaloGatherOpt) " <<timer0/calls <<std::endl;
    if ( timer1 ) std::cout << GridLogMessage << " timer1 (Communicate)   " <<timer1/calls <<std::endl;
    if ( timer2 ) std::cout << GridLogMessage << " timer2 (CommsMerge )   " <<timer2/calls <<std::endl;
    if ( timer3 ) std::cout << GridLogMessage << " timer3 (commsMergeShm) " <<timer3/calls <<std::endl;
    if ( timer4 ) std::cout << GridLogMessage << " timer4 " <<timer4 <<std::endl;
  }


  std::vector<int> same_node;
  std::vector<int> surface_list;

  WilsonStencilBasic(GridBase *grid,
		int npoints,
		int checkerboard,
		const std::vector<int> &directions,
		const std::vector<int> &distances)  
    : CartesianStencil<vobj,cobj> (grid,npoints,checkerboard,directions,distances) ,
    same_node(npoints)
  { 
    ZeroCountersi();
    surface_list.resize(0);
  };

  void BuildSurfaceList(int Ls,int vol4){

    // find same node for SHM
    // Here we know the distance is 1 for WilsonStencil
    for(int point=0;point<this->_npoints;point++){
      same_node[point] = this->SameNode(point);
    }
    
    for(int site = 0 ;site< vol4;site++){
      int local = 1;
      for(int point=0;point<this->_npoints;point++){
	if( (!this->GetNodeLocal(site*Ls,point)) && (!same_node[point]) ){ 
	  local = 0;
	}
      }
      if(local == 0) { 
	surface_list.push_back(site);
      }
    }
  }


  template < class compressor>
  void HaloExchangeOpt(const Lattice<vobj> &source,compressor &compress) 
  {
    std::vector<std::vector<CommsRequest_t> > reqs;
    this->HaloExchangeOptGather(source,compress);
    double t1=usecond();
    this->Communicate();
    double t2=usecond(); timer1 += t2-t1;
    this->CommsMerge(compress);
    double t3=usecond(); timer2 += t3-t2;
    this->CommsMergeSHM(compress);
    double t4=usecond(); timer3 += t4-t3;
  }
  
  template <class compressor>
  void HaloExchangeOptGather(const Lattice<vobj> &source,compressor &compress){
    this->Prepare();
    double t0=usecond();
    this->HaloGatherOpt(source,compress);
    double t1=usecond();
    timer0 += t1-t0;
    callsi++;
  }

  template <class compressor>
  void HaloGatherOpt(const Lattice<vobj> &source,compressor &compress)
  {
    this->halogtime-=usecond();
    this->HaloGather(source,compress);
    this->halogtime+=usecond();
  }
};




//This is hideous
template<class S, int packSize = 16>
class WilsonCompressedCommsImpl: public WilsonImpl<S,FundamentalRepresentation,CoeffReal>{
public:
  typedef WilsonImpl<S,FundamentalRepresentation,CoeffReal> WilsonBase;

#define INHERIT_BASE(TYPE) typedef typename WilsonBase::TYPE TYPE

  INHERIT_BASE(Gimpl);
  INHERIT_GIMPL_TYPES(Gimpl);

  INHERIT_BASE(Coeff_t);
  
  INHERIT_BASE(SiteSpinor);
  INHERIT_BASE(SitePropagator);
  INHERIT_BASE(SiteHalfSpinor);
  INHERIT_BASE(SiteHalfCommSpinor);    
  INHERIT_BASE(SiteDoubledGaugeField);

  INHERIT_BASE(FermionField);
  INHERIT_BASE(PropagatorField);
  INHERIT_BASE(DoubledGaugeField);

  //typedef WilsonCompressor<SiteHalfCommSpinor,SiteHalfSpinor, SiteSpinor> Compressor;
  typedef WilsonTestCompressor<SiteHalfSpinor, SiteSpinor, packSize> Compressor;  

  INHERIT_BASE(ImplParams);
  //INHERIT_BASE(StencilImpl);
  typedef WilsonStencilBasic<SiteSpinor, SiteHalfSpinor> StencilImpl;

  WilsonCompressedCommsImpl(const ImplParams &p = ImplParams()) : WilsonBase(p){}
  
  inline void multLink(SiteHalfSpinor &phi,
		       const SiteDoubledGaugeField &U,
		       const SiteHalfSpinor &chi,
		       int mu,
		       StencilEntry *SE,
		       StencilImpl &St) {
    mult(&phi(), &U(mu), &chi());
  }

#undef INHERIT_BASE
};

typedef WilsonCompressedCommsImpl<vComplexF,8> WilsonCompressedComms8ImplF;
typedef WilsonCompressedCommsImpl<vComplexD,8> WilsonCompressedComms8ImplD;
typedef WilsonCompressedCommsImpl<vComplexF,16> WilsonCompressedComms16ImplF;
typedef WilsonCompressedCommsImpl<vComplexD,16> WilsonCompressedComms16ImplD;


#define TO_INSTANTIATE \
  DOIT(WilsonCompressedComms8ImplF)\
  DOIT(WilsonCompressedComms8ImplD)\
  DOIT(WilsonCompressedComms16ImplF)\
  DOIT(WilsonCompressedComms16ImplD)

#include "InstantiateImpl.impl"

#undef TO_INSTANTIATE

typedef DomainWallFermion<WilsonCompressedComms8ImplD> DomainWallFermionCompressedComms8D;
typedef DomainWallFermion<WilsonCompressedComms8ImplF> DomainWallFermionCompressedComms8F;
typedef DomainWallFermion<WilsonCompressedComms16ImplD> DomainWallFermionCompressedComms16D;
typedef DomainWallFermion<WilsonCompressedComms16ImplF> DomainWallFermionCompressedComms16F;

template<typename T>
T parse(const std::string &name, std::istream &in){
  std::string p;
  in >> p;
  assert(p==name);
  char eq;
  in >> eq;
  assert(eq == '=');
  T out;
  in >> out;
  return out;
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int Ls=8;
  RealD mass=0.1;
  RealD outer_tol = 1e-8;
  RealD inner_tol_full = 1e-5;
  RealD inner_tol_half = 1e-5;
  RealD inner_tol_16c = 1e-5;
  RealD inner_tol_8c = 1e-5;

  RealD relup_delta_full = 0.1;
  RealD relup_delta_half = 0.1;
  RealD relup_delta_16c = 0.1;
  RealD relup_delta_8c = 0.1;

  std::string config_file = "";

  for(int i=1;i<argc;i++){
    if(std::string(argv[i]) == "--params"){
      std::ifstream f(argv[i+1]);
      f.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
      Ls = parse<int>("Ls", f);
#define PARSEIT(NM) NM = parse<RealD>(#NM, f)
      PARSEIT(mass);
      PARSEIT(outer_tol);
      PARSEIT(inner_tol_full);
      PARSEIT(inner_tol_half);
      PARSEIT(inner_tol_16c);
      PARSEIT(inner_tol_8c);
      PARSEIT(relup_delta_full);
      PARSEIT(relup_delta_half);
      PARSEIT(relup_delta_16c);
      PARSEIT(relup_delta_8c);
#undef PARSEIT

      //f >> outer_tol >> inner_tol_full >> inner_tol_half >> inner_tol_16c >> inner_tol_8c;      
    }else if(std::string(argv[i]) == "--config"){
      config_file = argv[i+1];
    }
  }
  
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexD::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  GridCartesian         * UGrid_f   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid_f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_f);
  GridCartesian         * FGrid_f   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_f);
  GridRedBlackCartesian * FrbGrid_f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_f);
  
  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeFermionD    src(FGrid); random(RNG5,src);
  LatticeFermionD result(FGrid); result=zero;
  LatticeGaugeFieldD Umu(UGrid);
  LatticeGaugeFieldF Umu_f(UGrid_f); 

  if(config_file.size() > 0){
    FieldMetaData header;
    NerscIO::readConfiguration(Umu,header,config_file);
  }else{
    SU3::HotConfiguration(RNG4,Umu);
  }

  precisionChange(Umu_f,Umu);
  
  RealD M5=1.8;

  LatticeFermionD    src_o(FrbGrid);
  pickCheckerboard(Odd,src_o,src);

  if(0){ //Test preconditioned CG
    LatticeFermionD result_o(FrbGrid);
    LatticeFermionD result_o_2(FrbGrid);
    result_o.checkerboard = Odd;
    result_o = zero;
    result_o_2.checkerboard = Odd;
    result_o_2 = zero;

    DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
    SchurDiagMooeeOperator<DomainWallFermionD,LatticeFermionD> HermOpEO(Ddwf);
    //DoNothingLinearOperator<LatticeFermionD> Prec;
    //FixedIterConjugateGradientPreconditioner<LatticeFermionD> Prec(HermOpEO, 20);
    SloppyConjugateGradientPreconditioner<LatticeFermionD> Prec(HermOpEO, 1e-2, 1000);

    std::cout << "Preconditioned CG" << std::endl;
    InexactPreconditionedConjugateGradient<LatticeFermionD> pCG(Prec,1.0e-8,10000);
    pCG(HermOpEO,src_o,result_o);

    std::cout << "Starting regular CG" << std::endl;
    ConjugateGradient<LatticeFermionD> CG(1.0e-8,10000);
    CG(HermOpEO,src_o,result_o_2);

    LatticeFermionD diff_o(FrbGrid);
    RealD diff = axpy_norm(diff_o, -1.0, result_o, result_o_2);

    std::cout << "pCG HermOp applications " << pCG.IterationsToComplete << "(outer) + " << Prec.InnerIterations << "(inner) = " << pCG.IterationsToComplete + Prec.InnerIterations << std::endl;
    std::cout << "CG HermOp applications " << CG.IterationsToComplete << std::endl;
    std::cout << "Diff between results: " << diff << std::endl;
  }

  if(0){ //Test compressor
    LatticeFermionD result_o(FrbGrid);
    LatticeFermionD result_o_2(FrbGrid);
    result_o.checkerboard = Odd;
    result_o = zero;
    result_o_2.checkerboard = Odd;
    result_o_2 = zero;

    DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
    SchurDiagMooeeOperator<DomainWallFermionD,LatticeFermionD> HermOpEO(Ddwf);

    DomainWallFermionCompressedComms16D DdwfC(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
    SchurDiagMooeeOperator<DomainWallFermionCompressedComms16D,LatticeFermionD> HermOpEOC(DdwfC);

    std::cout << "Starting regular CG with compressed operator" << std::endl;
    Integer iter1;
    {
      ConjugateGradient<LatticeFermionD> CG(1.0e-8,10000);
      CG.ErrorOnNoConverge = false;
      CG(HermOpEOC,src_o,result_o);
      iter1 = CG.IterationsToComplete;
    }
    Integer iter2;
    {
      std::cout << "Starting regular CG" << std::endl;
      ConjugateGradient<LatticeFermionD> CG(1.0e-8,10000);
      CG(HermOpEO,src_o,result_o_2);
      iter2 = CG.IterationsToComplete;
    }

    LatticeFermionD diff_o(FrbGrid);
    RealD diff = axpy_norm(diff_o, -1.0, result_o, result_o_2);

    std::cout << "CG HermOp CC applications " << iter1 << std::endl;
    std::cout << "CG HermOp applications " << iter2 << std::endl;
    std::cout << "Diff between results: " << diff << std::endl;
  }
  
  if(1){ //Compare mixed prec restarted single/single internal with same but with single/compressed
    LatticeFermionD result_o_full(FrbGrid);
    LatticeFermionD result_o_half(FrbGrid);
    LatticeFermionD result_o_16(FrbGrid);
    LatticeFermionD result_o_8(FrbGrid);
    result_o_full.checkerboard = Odd;
    result_o_full = zero;
    result_o_16 = result_o_8 = result_o_half = result_o_full;

    //Std
    DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
    SchurDiagMooeeOperator<DomainWallFermionD,LatticeFermionD> HermOpEO(Ddwf);

    DomainWallFermionF Ddwf_f(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5);
    SchurDiagMooeeOperator<DomainWallFermionF,LatticeFermionF> HermOpEO_f(Ddwf_f);

    //1/2 prec
    DomainWallFermionFH Ddwfhalf_f(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5);
    SchurDiagMooeeOperator<DomainWallFermionFH,LatticeFermionF> HermOpEOhalf_f(Ddwfhalf_f);

    //16
    DomainWallFermionCompressedComms16F DdwfC16_f(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5);
    SchurDiagMooeeOperator<DomainWallFermionCompressedComms16F,LatticeFermionF> HermOpEOC16_f(DdwfC16_f);

    //8
    DomainWallFermionCompressedComms8F DdwfC8_f(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5);
    SchurDiagMooeeOperator<DomainWallFermionCompressedComms8F,LatticeFermionF> HermOpEOC8_f(DdwfC8_f);

    #define ALGORITHM_MIXEDCG
    //#define ALGORITHM_RELUP
    //#define ALGORITHM_SLOPPY_PREC_CG

#ifdef ALGORITHM_MIXEDCG
    std::cout << "Starting mixed CG with single/compressed-16 inner\n";

    Integer inner_16, outer_16, patchup_16;
    {
      MixedPrecisionConjugateGradient<LatticeFermionD,LatticeFermionF> mCG(outer_tol, 10000, 50, FrbGrid_f, HermOpEOC16_f, HermOpEO);
      mCG.InnerTolerance = inner_tol_16c;
      mCG(src_o,result_o_16);
      inner_16 = mCG.TotalInnerIterations; outer_16 = mCG.TotalOuterIterations; patchup_16 = mCG.TotalFinalStepIterations;
    }

    std::cout << "Starting mixed CG with single/compressed-8 inner\n";
    Integer inner_8, outer_8, patchup_8;
    {
      MixedPrecisionConjugateGradient<LatticeFermionD,LatticeFermionF> mCG(outer_tol, 10000, 50, FrbGrid_f, HermOpEOC8_f, HermOpEO);
      mCG.InnerTolerance = inner_tol_8c;
      mCG(src_o,result_o_8);
      inner_8 = mCG.TotalInnerIterations; outer_8 = mCG.TotalOuterIterations; patchup_8 = mCG.TotalFinalStepIterations;
    }

    std::cout << "Starting mixed CG with single/half inner\n";
    Integer inner_half, outer_half, patchup_half;
    {
      MixedPrecisionConjugateGradient<LatticeFermionD,LatticeFermionF> mCG(outer_tol, 10000, 50, FrbGrid_f, HermOpEOhalf_f, HermOpEO);
      mCG.InnerTolerance = inner_tol_half;
      mCG(src_o,result_o_half);
      inner_half = mCG.TotalInnerIterations; outer_half = mCG.TotalOuterIterations; patchup_half = mCG.TotalFinalStepIterations;
    }

    std::cout << "Starting mixed CG with single/single inner\n";
    Integer inner_full, outer_full, patchup_full;
    {
      MixedPrecisionConjugateGradient<LatticeFermionD,LatticeFermionF> mCG(outer_tol, 10000, 50, FrbGrid_f, HermOpEO_f, HermOpEO);
      mCG.InnerTolerance = inner_tol_full;
      mCG(src_o,result_o_full);
      inner_full = mCG.TotalInnerIterations; outer_full = mCG.TotalOuterIterations; patchup_full = mCG.TotalFinalStepIterations;
    }
#elif defined(ALGORITHM_RELUP)
    std::cout << "Starting relup CG with single/compressed-16 inner\n";    
    Integer inner_16, outer_16, patchup_16;
    {
      ConjugateGradientReliableUpdate<LatticeFermionD,LatticeFermionF> relup(outer_tol, 2000, relup_delta_16c, FrbGrid_f, HermOpEOC16_f, HermOpEO);
      relup(src_o,result_o_16);
      inner_16 = relup.IterationsToComplete; outer_16 = relup.ReliableUpdatesPerformed; patchup_16 = relup.IterationsToCleanup;
    }

    std::cout << "Starting relup CG with single/compressed-8 inner\n";
    Integer inner_8, outer_8, patchup_8;
    {
      ConjugateGradientReliableUpdate<LatticeFermionD,LatticeFermionF> relup(outer_tol, 2000, relup_delta_8c, FrbGrid_f, HermOpEOC8_f, HermOpEO);
      relup.ErrorOnNoConverge = false;
      relup(src_o,result_o_8);
      inner_8 = relup.IterationsToComplete; outer_8 = relup.ReliableUpdatesPerformed; patchup_8 = relup.IterationsToCleanup;     
    }

    std::cout << "Starting relup CG with single/half inner\n";
    Integer inner_half, outer_half, patchup_half;
    {
      ConjugateGradientReliableUpdate<LatticeFermionD,LatticeFermionF> relup(outer_tol, 2000, relup_delta_half, FrbGrid_f, HermOpEOhalf_f, HermOpEO);
      relup(src_o,result_o_half);
      inner_half = relup.IterationsToComplete; outer_half = relup.ReliableUpdatesPerformed; patchup_half = relup.IterationsToCleanup;
    }

    std::cout << "Starting relup CG with single/single inner\n";
    Integer inner_full, outer_full, patchup_full;
    {
      ConjugateGradientReliableUpdate<LatticeFermionD,LatticeFermionF> relup(outer_tol, 2000, relup_delta_full, FrbGrid_f, HermOpEO_f, HermOpEO);
      relup(src_o,result_o_full);
      inner_full = relup.IterationsToComplete; outer_full = relup.ReliableUpdatesPerformed; patchup_full = relup.IterationsToCleanup;
    }
#elif defined(ALGORITHM_SLOPPY_PREC_CG)

    std::cout << "Starting sloppy pCG with single/compressed-16 inner\n";    
    Integer inner_16, outer_16;
    {
      SloppyConjugateGradientLowerPrecPreconditioner<LatticeFermionD,LatticeFermionF> prec(HermOpEOC16_f, FrbGrid_f, inner_tol_16c, 1000);
      InexactPreconditionedConjugateGradient<LatticeFermionD> CG(prec, outer_tol, 100);
      CG(HermOpEO,src_o,result_o_16);
      inner_16 = prec.InnerIterations; outer_16 = CG.IterationsToComplete;
    }

    std::cout << "Starting sloppy pCG with single/compressed-8 inner\n";    
    Integer inner_8, outer_8;
    {
      SloppyConjugateGradientLowerPrecPreconditioner<LatticeFermionD,LatticeFermionF> prec(HermOpEOC8_f, FrbGrid_f, inner_tol_8c, 1000);
      InexactPreconditionedConjugateGradient<LatticeFermionD> CG(prec, outer_tol, 100);
      CG(HermOpEO,src_o,result_o_8);
      inner_8 = prec.InnerIterations; outer_8 = CG.IterationsToComplete;
    }

    std::cout << "Starting sloppy pCG with single/half inner\n";    
    Integer inner_half, outer_half;
    {
      SloppyConjugateGradientLowerPrecPreconditioner<LatticeFermionD,LatticeFermionF> prec(HermOpEOhalf_f, FrbGrid_f, inner_tol_half, 1000);
      InexactPreconditionedConjugateGradient<LatticeFermionD> CG(prec, outer_tol, 100);
      CG(HermOpEO,src_o,result_o_half);
      inner_half = prec.InnerIterations; outer_half = CG.IterationsToComplete;
    }

    std::cout << "Starting sloppy pCG with single/single inner\n";    
    Integer inner_full, outer_full;
    {
      SloppyConjugateGradientLowerPrecPreconditioner<LatticeFermionD,LatticeFermionF> prec(HermOpEO_f, FrbGrid_f, inner_tol_full, 1000);
      InexactPreconditionedConjugateGradient<LatticeFermionD> CG(prec, outer_tol, 100);
      CG(HermOpEO,src_o,result_o_full);
      inner_full = prec.InnerIterations; outer_full = CG.IterationsToComplete;
    }

#endif

    std::cout << "Ls " << Ls << std::endl;
    std::cout << "Mass " << mass << std::endl;
    std::cout << "Outer tolerance " << outer_tol << std::endl;

#if defined(ALGORITHM_MIXEDCG) || defined(ALGORITHM_SLOPPY_PREC_CG)
    std::cout << "Inner tol full " << inner_tol_full << std::endl;
    std::cout << "Inner tol 1/2 prec " << inner_tol_half << std::endl;
    std::cout << "Inner tol compressed-16 " << inner_tol_16c << std::endl;
    std::cout << "Inner tol compressed-8 " << inner_tol_8c << std::endl;
#elif defined(ALGORITHM_RELUP)
    std::cout << "Relup delta full " << relup_delta_full << std::endl;
    std::cout << "Relup delta 1/2 prec " << relup_delta_half << std::endl;
    std::cout << "Relup delta compressed-16 " << relup_delta_16c << std::endl;
    std::cout << "Relup delta compressed-8 " << relup_delta_8c << std::endl;
#endif
    LatticeFermionD diff_o(FrbGrid);
    RealD diff = axpy_norm(diff_o, -1.0, result_o_16, result_o_full);
    std::cout << "Diff between results (s/c16): " << diff << std::endl;

    diff = axpy_norm(diff_o, -1.0, result_o_8, result_o_full);
    std::cout << "Diff between results (s/c8): " << diff << std::endl;

    diff = axpy_norm(diff_o, -1.0, result_o_half, result_o_full);
    std::cout << "Diff between results (s/h): " << diff << std::endl;

#if defined(ALGORITHM_MIXEDCG) || defined(ALGORITHM_RELUP)
    std::cout << "Iterations (s/c16) inner: " << inner_16 << " outer: " << outer_16 << " patchup: " << patchup_16 << std::endl;
    std::cout << "Iterations (s/c8) inner: " << inner_8 << " outer: " << outer_8 << " patchup: " << patchup_8 << std::endl;
    std::cout << "Iterations (s/h) inner: " << inner_half << " outer: " << outer_half << " patchup: " << patchup_half << std::endl;
    std::cout << "Iterations (s/s) inner: " << inner_full << " outer: " << outer_full << " patchup: " << patchup_full << std::endl;
#else
    std::cout << "Iterations (s/c16) inner: " << inner_16 << " outer: " << outer_16 << std::endl;
    std::cout << "Iterations (s/c8) inner: " << inner_8 << " outer: " << outer_8 << std::endl;
    std::cout << "Iterations (s/h) inner: " << inner_half << " outer: " << outer_half << std::endl;
    std::cout << "Iterations (s/s) inner: " << inner_full << " outer: " << outer_full << std::endl;
#endif

  }


  Grid_finalize();
}
