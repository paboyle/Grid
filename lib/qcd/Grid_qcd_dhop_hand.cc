#include <Grid.h>

#define REGISTER

#define LOAD_CHIMU \
  const vSpinColourVector & ref (in._odata[offset]);	\
    Chimu_00=ref()(0)(0);\
    Chimu_01=ref()(0)(1);\
    Chimu_02=ref()(0)(2);\
    Chimu_10=ref()(1)(0);\
    Chimu_11=ref()(1)(1);\
    Chimu_12=ref()(1)(2);\
    Chimu_20=ref()(2)(0);\
    Chimu_21=ref()(2)(1);\
    Chimu_22=ref()(2)(2);\
    Chimu_30=ref()(3)(0);\
    Chimu_31=ref()(3)(1);\
    Chimu_32=ref()(3)(2);

#define LOAD_CHI\
  const vHalfSpinColourVector &ref(buf[offset]);	\
    Chi_00 = ref()(0)(0);\
    Chi_01 = ref()(0)(1);\
    Chi_02 = ref()(0)(2);\
    Chi_10 = ref()(1)(0);\
    Chi_11 = ref()(1)(1);\
    Chi_12 = ref()(1)(2);

#define MULT_2SPIN(A)\
   auto & ref(U._odata[ss](A));	\
    U_00 = ref()(0,0);\
    U_10 = ref()(1,0);\
    U_20 = ref()(2,0);\
    U_01 = ref()(0,1);\
    U_11 = ref()(1,1);				\
    U_21 = ref()(2,1);\
    UChi_00 = U_00*Chi_00;\
    UChi_10 = U_00*Chi_10;\
    UChi_01 = U_10*Chi_00;\
    UChi_11 = U_10*Chi_10;\
    UChi_02 = U_20*Chi_00;\
    UChi_12 = U_20*Chi_10;\
    UChi_00+= U_01*Chi_01;\
    UChi_10+= U_01*Chi_11;\
    UChi_01+= U_11*Chi_01;\
    UChi_11+= U_11*Chi_11;\
    UChi_02+= U_21*Chi_01;\
    UChi_12+= U_21*Chi_11;\
    U_00 = ref()(0,2);\
    U_10 = ref()(1,2);\
    U_20 = ref()(2,2);\
    UChi_00+= U_00*Chi_02;\
    UChi_10+= U_00*Chi_12;\
    UChi_01+= U_10*Chi_02;\
    UChi_11+= U_10*Chi_12;\
    UChi_02+= U_20*Chi_02;\
    UChi_12+= U_20*Chi_12;

#define PERMUTE\
      permute(Chi_00,Chi_00,ptype);\
      permute(Chi_01,Chi_01,ptype);\
      permute(Chi_02,Chi_02,ptype);\
      permute(Chi_10,Chi_10,ptype);\
      permute(Chi_11,Chi_11,ptype);\
      permute(Chi_12,Chi_12,ptype);

//      hspin(0)=fspin(0)+timesI(fspin(3));
//      hspin(1)=fspin(1)+timesI(fspin(2));
#define XP_PROJ \
    Chi_00 = Chimu_00+timesI(Chimu_30);\
    Chi_01 = Chimu_01+timesI(Chimu_31);\
    Chi_02 = Chimu_02+timesI(Chimu_32);\
    Chi_10 = Chimu_10+timesI(Chimu_20);\
    Chi_11 = Chimu_11+timesI(Chimu_21);\
    Chi_12 = Chimu_12+timesI(Chimu_22);

#define YP_PROJ \
    Chi_00 = Chimu_00-Chimu_30;\
    Chi_01 = Chimu_01-Chimu_31;\
    Chi_02 = Chimu_02-Chimu_32;\
    Chi_10 = Chimu_10+Chimu_20;\
    Chi_11 = Chimu_11+Chimu_21;\
    Chi_12 = Chimu_12+Chimu_22;

#define ZP_PROJ \
  Chi_00 = Chimu_00+timesI(Chimu_20);		\
  Chi_01 = Chimu_01+timesI(Chimu_21);		\
  Chi_02 = Chimu_02+timesI(Chimu_22);		\
  Chi_10 = Chimu_10-timesI(Chimu_30);		\
  Chi_11 = Chimu_11-timesI(Chimu_31);		\
  Chi_12 = Chimu_12-timesI(Chimu_32);

#define TP_PROJ \
  Chi_00 = Chimu_00+Chimu_20;		\
  Chi_01 = Chimu_01+Chimu_21;		\
  Chi_02 = Chimu_02+Chimu_22;		\
  Chi_10 = Chimu_10+Chimu_30;		\
  Chi_11 = Chimu_11+Chimu_31;		\
  Chi_12 = Chimu_12+Chimu_32;


//      hspin(0)=fspin(0)-timesI(fspin(3));
//      hspin(1)=fspin(1)-timesI(fspin(2));
#define XM_PROJ \
    Chi_00 = Chimu_00-timesI(Chimu_30);\
    Chi_01 = Chimu_01-timesI(Chimu_31);\
    Chi_02 = Chimu_02-timesI(Chimu_32);\
    Chi_10 = Chimu_10-timesI(Chimu_20);\
    Chi_11 = Chimu_11-timesI(Chimu_21);\
    Chi_12 = Chimu_12-timesI(Chimu_22);

#define YM_PROJ \
    Chi_00 = Chimu_00+Chimu_30;\
    Chi_01 = Chimu_01+Chimu_31;\
    Chi_02 = Chimu_02+Chimu_32;\
    Chi_10 = Chimu_10-Chimu_20;\
    Chi_11 = Chimu_11-Chimu_21;\
    Chi_12 = Chimu_12-Chimu_22;

#define ZM_PROJ \
  Chi_00 = Chimu_00-timesI(Chimu_20);		\
  Chi_01 = Chimu_01-timesI(Chimu_21);		\
  Chi_02 = Chimu_02-timesI(Chimu_22);		\
  Chi_10 = Chimu_10+timesI(Chimu_30);		\
  Chi_11 = Chimu_11+timesI(Chimu_31);		\
  Chi_12 = Chimu_12+timesI(Chimu_32);

#define TM_PROJ \
  Chi_00 = Chimu_00-Chimu_20;		\
  Chi_01 = Chimu_01-Chimu_21;		\
  Chi_02 = Chimu_02-Chimu_22;		\
  Chi_10 = Chimu_10-Chimu_30;		\
  Chi_11 = Chimu_11-Chimu_31;		\
  Chi_12 = Chimu_12-Chimu_32;

//      fspin(0)=hspin(0);
//      fspin(1)=hspin(1);
//      fspin(2)=timesMinusI(hspin(1));
//      fspin(3)=timesMinusI(hspin(0));
#define XP_RECON\
  result_00 = UChi_00;\
  result_01 = UChi_01;\
  result_02 = UChi_02;\
  result_10 = UChi_10;\
  result_11 = UChi_11;\
  result_12 = UChi_12;\
  result_20 = timesMinusI(UChi_10);\
  result_21 = timesMinusI(UChi_11);\
  result_22 = timesMinusI(UChi_12);\
  result_30 = timesMinusI(UChi_00);\
  result_31 = timesMinusI(UChi_01);\
  result_32 = timesMinusI(UChi_02);

#define XP_RECON_ACCUM\
  result_00+=UChi_00;\
  result_01+=UChi_01;\
  result_02+=UChi_02;\
  result_10+=UChi_10;\
  result_11+=UChi_11;\
  result_12+=UChi_12;\
  result_20-=timesI(UChi_10);\
  result_21-=timesI(UChi_11);\
  result_22-=timesI(UChi_12);\
  result_30-=timesI(UChi_00);\
  result_31-=timesI(UChi_01);\
  result_32-=timesI(UChi_02);

#define XM_RECON\
  result_00 = UChi_00;\
  result_01 = UChi_01;\
  result_02 = UChi_02;\
  result_10 = UChi_10;\
  result_11 = UChi_11;\
  result_12 = UChi_12;\
  result_20 = timesI(UChi_10);\
  result_21 = timesI(UChi_11);\
  result_22 = timesI(UChi_12);\
  result_30 = timesI(UChi_00);\
  result_31 = timesI(UChi_01);\
  result_32 = timesI(UChi_02);

#define XM_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20+= timesI(UChi_10);\
  result_21+= timesI(UChi_11);\
  result_22+= timesI(UChi_12);\
  result_30+= timesI(UChi_00);\
  result_31+= timesI(UChi_01);\
  result_32+= timesI(UChi_02);

#define YP_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20+= UChi_10;\
  result_21+= UChi_11;\
  result_22+= UChi_12;\
  result_30-= UChi_00;\
  result_31-= UChi_01;\
  result_32-= UChi_02;

#define YM_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20-= UChi_10;\
  result_21-= UChi_11;\
  result_22-= UChi_12;\
  result_30+= UChi_00;\
  result_31+= UChi_01;\
  result_32+= UChi_02;

#define ZP_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20-= timesI(UChi_00);			\
  result_21-= timesI(UChi_01);			\
  result_22-= timesI(UChi_02);			\
  result_30+= timesI(UChi_10);			\
  result_31+= timesI(UChi_11);			\
  result_32+= timesI(UChi_12);

#define ZM_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20+= timesI(UChi_00);			\
  result_21+= timesI(UChi_01);			\
  result_22+= timesI(UChi_02);			\
  result_30-= timesI(UChi_10);			\
  result_31-= timesI(UChi_11);			\
  result_32-= timesI(UChi_12);

#define TP_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20+= UChi_00;			\
  result_21+= UChi_01;			\
  result_22+= UChi_02;			\
  result_30+= UChi_10;			\
  result_31+= UChi_11;			\
  result_32+= UChi_12;

#define TM_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20-= UChi_00;	\
  result_21-= UChi_01;	\
  result_22-= UChi_02;	\
  result_30-= UChi_10;	\
  result_31-= UChi_11;	\
  result_32-= UChi_12;

namespace Grid {
namespace QCD {

void DiracOptHand::DhopSite(CartesianStencil &st,LatticeDoubledGaugeField &U,
			    std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  &buf,
			    int ss,const LatticeFermion &in, LatticeFermion &out)
{
  REGISTER vComplex result_00; // 12 regs on knc
  REGISTER vComplex result_01;
  REGISTER vComplex result_02;
  
  REGISTER vComplex result_10;
  REGISTER vComplex result_11;
  REGISTER vComplex result_12;

  REGISTER vComplex result_20;
  REGISTER vComplex result_21;
  REGISTER vComplex result_22;

  REGISTER vComplex result_30;
  REGISTER vComplex result_31;
  REGISTER vComplex result_32; // 20 left

  REGISTER vComplex Chi_00;    // two spinor; 6 regs
  REGISTER vComplex Chi_01;
  REGISTER vComplex Chi_02;

  REGISTER vComplex Chi_10;
  REGISTER vComplex Chi_11;
  REGISTER vComplex Chi_12;   // 14 left

  REGISTER vComplex UChi_00;  // two spinor; 6 regs
  REGISTER vComplex UChi_01;
  REGISTER vComplex UChi_02;

  REGISTER vComplex UChi_10;
  REGISTER vComplex UChi_11;
  REGISTER vComplex UChi_12;  // 8 left

  REGISTER vComplex U_00;  // two rows of U matrix
  REGISTER vComplex U_10;
  REGISTER vComplex U_20;  
  REGISTER vComplex U_01;
  REGISTER vComplex U_11;
  REGISTER vComplex U_21;  // 2 reg left.

#define Chimu_00 Chi_00
#define Chimu_01 Chi_01
#define Chimu_02 Chi_02
#define Chimu_10 Chi_10
#define Chimu_11 Chi_11
#define Chimu_12 Chi_12
#define Chimu_20 UChi_00
#define Chimu_21 UChi_01
#define Chimu_22 UChi_02
#define Chimu_30 UChi_10
#define Chimu_31 UChi_11
#define Chimu_32 UChi_12


  int offset,local,perm, ptype;

  // Xp
  offset = st._offsets [Xp][ss];
  local  = st._is_local[Xp][ss];
  perm   = st._permute[Xp][ss];
  ptype  = st._permute_type[Xp];
  
  if ( local ) {
    LOAD_CHIMU;
    XP_PROJ;
    if ( perm) {
      PERMUTE;
    }
  } else { 
    LOAD_CHI;
  }

  {
    MULT_2SPIN(Xp);
  }
  XP_RECON;
  //  std::cout << "XP_RECON"<<std::endl;
  //  std::cout << result_00 <<" "<<result_01 <<" "<<result_02 <<std::endl;
  //  std::cout << result_10 <<" "<<result_11 <<" "<<result_12 <<std::endl;
  //  std::cout << result_20 <<" "<<result_21 <<" "<<result_22 <<std::endl;
  //  std::cout << result_30 <<" "<<result_31 <<" "<<result_32 <<std::endl;

  // Yp
  offset = st._offsets [Yp][ss];
  local  = st._is_local[Yp][ss];
  perm   = st._permute[Yp][ss];
  ptype  = st._permute_type[Yp];
  
  if ( local ) {
    LOAD_CHIMU;
    YP_PROJ;
    if ( perm) {
      PERMUTE;
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Yp);
  }
  YP_RECON_ACCUM;


  // Zp
  offset = st._offsets [Zp][ss];
  local  = st._is_local[Zp][ss];
  perm   = st._permute[Zp][ss];
  ptype  = st._permute_type[Zp];
  
  if ( local ) {
    LOAD_CHIMU;
    ZP_PROJ;
    if ( perm) {
      PERMUTE;
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Zp);
  }
  ZP_RECON_ACCUM;

  // Tp
  offset = st._offsets [Tp][ss];
  local  = st._is_local[Tp][ss];
  perm   = st._permute[Tp][ss];
  ptype  = st._permute_type[Tp];
  
  if ( local ) {
    LOAD_CHIMU;
    TP_PROJ;
    if ( perm) {
      PERMUTE;
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Tp);
  }
  TP_RECON_ACCUM;
  
  // Xm
  offset = st._offsets [Xm][ss];
  local  = st._is_local[Xm][ss];
  perm   = st._permute[Xm][ss];
  ptype  = st._permute_type[Xm];
  
  if ( local ) {
    LOAD_CHIMU;
    XM_PROJ;
    if ( perm) {
      PERMUTE;
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Xm);
  }
  XM_RECON_ACCUM;
  //  std::cout << "XM_RECON_ACCUM"<<std::endl;
  //  std::cout << result_00 <<" "<<result_01 <<" "<<result_02 <<std::endl;
  //  std::cout << result_10 <<" "<<result_11 <<" "<<result_12 <<std::endl;
  //  std::cout << result_20 <<" "<<result_21 <<" "<<result_22 <<std::endl;
  //  std::cout << result_30 <<" "<<result_31 <<" "<<result_32 <<std::endl;
  
  
  // Ym
  offset = st._offsets [Ym][ss];
  local  = st._is_local[Ym][ss];
  perm   = st._permute[Ym][ss];
  ptype  = st._permute_type[Ym];
  
  if ( local ) {
    LOAD_CHIMU;
    YM_PROJ;
    if ( perm) {
      PERMUTE;
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Ym);
  }
  YM_RECON_ACCUM;

  // Zm
  offset = st._offsets [Zm][ss];
  local  = st._is_local[Zm][ss];
  perm   = st._permute[Zm][ss];
  ptype  = st._permute_type[Zm];

  if ( local ) {
    LOAD_CHIMU;
    ZM_PROJ;
    if ( perm) {
      PERMUTE;
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Zm);
  }
  ZM_RECON_ACCUM;

  // Tm
  offset = st._offsets [Tm][ss];
  local  = st._is_local[Tm][ss];
  perm   = st._permute[Tm][ss];
  ptype  = st._permute_type[Tm];

  if ( local ) {
    LOAD_CHIMU;
    TM_PROJ;
    if ( perm) {
      PERMUTE;
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Tm);
  }
  TM_RECON_ACCUM;

  {
    vSpinColourVector & ref (out._odata[ss]);
    vstream(ref()(0)(0),result_00);
    vstream(ref()(0)(1),result_01);
    vstream(ref()(0)(2),result_02);
    vstream(ref()(1)(0),result_10);
    vstream(ref()(1)(1),result_11);
    vstream(ref()(1)(2),result_12);
    vstream(ref()(2)(0),result_20);
    vstream(ref()(2)(1),result_21);
    vstream(ref()(2)(2),result_22);
    vstream(ref()(3)(0),result_30);
    vstream(ref()(3)(1),result_31);
    vstream(ref()(3)(2),result_32);
  }
}

void DiracOptHand::DhopSiteDag(CartesianStencil &st,LatticeDoubledGaugeField &U,
			       std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  &buf,
			       int ss,const LatticeFermion &in, LatticeFermion &out)
{
  REGISTER vComplex result_00; // 12 regs on knc
  REGISTER vComplex result_01;
  REGISTER vComplex result_02;

  REGISTER vComplex result_10;
  REGISTER vComplex result_11;
  REGISTER vComplex result_12;

  REGISTER vComplex result_20;
  REGISTER vComplex result_21;
  REGISTER vComplex result_22;

  REGISTER vComplex result_30;
  REGISTER vComplex result_31;
  REGISTER vComplex result_32; // 20 left

  REGISTER vComplex Chi_00;    // two spinor; 6 regs
  REGISTER vComplex Chi_01;
  REGISTER vComplex Chi_02;

  REGISTER vComplex Chi_10;
  REGISTER vComplex Chi_11;
  REGISTER vComplex Chi_12;   // 14 left

  REGISTER vComplex UChi_00;  // two spinor; 6 regs
  REGISTER vComplex UChi_01;
  REGISTER vComplex UChi_02;

  REGISTER vComplex UChi_10;
  REGISTER vComplex UChi_11;
  REGISTER vComplex UChi_12;  // 8 left

  REGISTER vComplex U_00;  // two rows of U matrix
  REGISTER vComplex U_10;
  REGISTER vComplex U_20;  
  REGISTER vComplex U_01;
  REGISTER vComplex U_11;
  REGISTER vComplex U_21;  // 2 reg left.

#define Chimu_00 Chi_00
#define Chimu_01 Chi_01
#define Chimu_02 Chi_02
#define Chimu_10 Chi_10
#define Chimu_11 Chi_11
#define Chimu_12 Chi_12
#define Chimu_20 UChi_00
#define Chimu_21 UChi_01
#define Chimu_22 UChi_02
#define Chimu_30 UChi_10
#define Chimu_31 UChi_11
#define Chimu_32 UChi_12


  int offset,local,perm, ptype;

  // Xp
  offset = st._offsets [Xp][ss];
  local  = st._is_local[Xp][ss];
  perm   = st._permute[Xp][ss];
  ptype  = st._permute_type[Xp];
  
  if ( local ) {
    LOAD_CHIMU;
    XM_PROJ;
    if ( perm) {
      PERMUTE;
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Xp);
  }
  XM_RECON;
  
  // Yp
  offset = st._offsets [Yp][ss];
  local  = st._is_local[Yp][ss];
  perm   = st._permute[Yp][ss];
  ptype  = st._permute_type[Yp];
  
  if ( local ) {
    LOAD_CHIMU;
    YM_PROJ;
    if ( perm) {
      PERMUTE;
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Yp);
  }
  YM_RECON_ACCUM;


  // Zp
  offset = st._offsets [Zp][ss];
  local  = st._is_local[Zp][ss];
  perm   = st._permute[Zp][ss];
  ptype  = st._permute_type[Zp];
  
  if ( local ) {
    LOAD_CHIMU;
    ZM_PROJ;
    if ( perm) {
      PERMUTE;
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Zp);
  }
  ZM_RECON_ACCUM;

  // Tp
  offset = st._offsets [Tp][ss];
  local  = st._is_local[Tp][ss];
  perm   = st._permute[Tp][ss];
  ptype  = st._permute_type[Tp];
  
  if ( local ) {
    LOAD_CHIMU;
    TM_PROJ;
    if ( perm) {
      PERMUTE;
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Tp);
  }
  TM_RECON_ACCUM;
  
  // Xm
  offset = st._offsets [Xm][ss];
  local  = st._is_local[Xm][ss];
  perm   = st._permute[Xm][ss];
  ptype  = st._permute_type[Xm];
  
  if ( local ) {
    LOAD_CHIMU;
    XP_PROJ;
    if ( perm) {
      PERMUTE;
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Xm);
  }
  XP_RECON_ACCUM;
  
  
  // Ym
  offset = st._offsets [Ym][ss];
  local  = st._is_local[Ym][ss];
  perm   = st._permute[Ym][ss];
  ptype  = st._permute_type[Ym];
  
  if ( local ) {
    LOAD_CHIMU;
    YP_PROJ;
    if ( perm) {
      PERMUTE;
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Ym);
  }
  YP_RECON_ACCUM;

  // Zm
  offset = st._offsets [Zm][ss];
  local  = st._is_local[Zm][ss];
  perm   = st._permute[Zm][ss];
  ptype  = st._permute_type[Zm];

  if ( local ) {
    LOAD_CHIMU;
    ZP_PROJ;
    if ( perm) {
      PERMUTE;
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Zm);
  }
  ZP_RECON_ACCUM;

  // Tm
  offset = st._offsets [Tm][ss];
  local  = st._is_local[Tm][ss];
  perm   = st._permute[Tm][ss];
  ptype  = st._permute_type[Tm];

  if ( local ) {
    LOAD_CHIMU;
    TP_PROJ;
    if ( perm) {
      PERMUTE;
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Tm);
  }
  TP_RECON_ACCUM;

  {
    vSpinColourVector & ref (out._odata[ss]);
    vstream(ref()(0)(0),result_00);
    vstream(ref()(0)(1),result_01);
    vstream(ref()(0)(2),result_02);
    vstream(ref()(1)(0),result_10);
    vstream(ref()(1)(1),result_11);
    vstream(ref()(1)(2),result_12);
    vstream(ref()(2)(0),result_20);
    vstream(ref()(2)(1),result_21);
    vstream(ref()(2)(2),result_22);
    vstream(ref()(3)(0),result_30);
    vstream(ref()(3)(1),result_31);
    vstream(ref()(3)(2),result_32);
  }
}
}}
