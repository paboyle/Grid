/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: WilsonKernelsAsmBodyA64FX.h

    Copyright (C) 2020

Author:  Nils Meyer  <nils.meyer@ur.de>  Regensburg University

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

// GCC 10 messes up SVE instruction scheduling using -O3, but
// -O3 -fno-schedule-insns -fno-schedule-insns2 does wonders
// performance now is better than armclang 20.2

#ifdef KERNEL_DAG
#define DIR0_PROJ    XP_PROJ
#define DIR1_PROJ    YP_PROJ
#define DIR2_PROJ    ZP_PROJ
#define DIR3_PROJ    TP_PROJ
#define DIR4_PROJ    XM_PROJ
#define DIR5_PROJ    YM_PROJ
#define DIR6_PROJ    ZM_PROJ
#define DIR7_PROJ    TM_PROJ
#define DIR0_RECON   XP_RECON
#define DIR1_RECON   YP_RECON_ACCUM
#define DIR2_RECON   ZP_RECON_ACCUM
#define DIR3_RECON   TP_RECON_ACCUM
#define DIR4_RECON   XM_RECON_ACCUM
#define DIR5_RECON   YM_RECON_ACCUM
#define DIR6_RECON   ZM_RECON_ACCUM
#define DIR7_RECON   TM_RECON_ACCUM
#else
#define DIR0_PROJ    XM_PROJ
#define DIR1_PROJ    YM_PROJ
#define DIR2_PROJ    ZM_PROJ
#define DIR3_PROJ    TM_PROJ
#define DIR4_PROJ    XP_PROJ
#define DIR5_PROJ    YP_PROJ
#define DIR6_PROJ    ZP_PROJ
#define DIR7_PROJ    TP_PROJ
#define DIR0_RECON   XM_RECON
#define DIR1_RECON   YM_RECON_ACCUM
#define DIR2_RECON   ZM_RECON_ACCUM
#define DIR3_RECON   TM_RECON_ACCUM
#define DIR4_RECON   XP_RECON_ACCUM
#define DIR5_RECON   YP_RECON_ACCUM
#define DIR6_RECON   ZP_RECON_ACCUM
#define DIR7_RECON   TP_RECON_ACCUM
#endif

//using namespace std;

#undef SHOW
//#define SHOW

#undef WHERE

#ifdef INTERIOR_AND_EXTERIOR
#define WHERE "INT_AND_EXT"
#endif

#ifdef INTERIOR
#define WHERE "INT"
#endif

#ifdef EXTERIOR
#define WHERE "EXT"
#endif

//#pragma message("here")



////////////////////////////////////////////////////////////////////////////////
// Comms then compute kernel
////////////////////////////////////////////////////////////////////////////////
#ifdef INTERIOR_AND_EXTERIOR

#define ASM_LEG(Dir,NxtDir,PERMUTE_DIR,PROJ,RECON)			\
      basep = st.GetPFInfo(nent,plocal); nent++;			\
      if ( local ) {							            \
    LOAD_CHIMU(base);                                       \
    LOAD_TABLE(PERMUTE_DIR);                                \
    PROJ;							                        \
    MAYBEPERM(PERMUTE_DIR,perm);					        \
      } else {								                \
	  LOAD_CHI(base);							                \
      }									                    \
      base = st.GetInfo(ptype,local,perm,NxtDir,ent,plocal); ent++;	\
    MULT_2SPIN_1(Dir);					                    \
    PREFETCH_CHIMU(base);                                   \
    PREFETCH_CHIMU_L2(basep);                               \
    /* PREFETCH_GAUGE_L1(NxtDir); */                        \
    MULT_2SPIN_2;					                        \
    if (s == 0) {                                           \
      if ((Dir == 0) || (Dir == 4)) { PREFETCH_GAUGE_L2(Dir); } \
    }                                                       \
    RECON;								                    \

/*
NB: picking PREFETCH_GAUGE_L2(Dir+4); here results in performance penalty
    though I expected that it would improve on performance
*/

#define ASM_LEG_XP(Dir,NxtDir,PERMUTE_DIR,PROJ,RECON)	    \
  base = st.GetInfo(ptype,local,perm,Dir,ent,plocal); ent++; \
  PREFETCH1_CHIMU(base);						            \
  ASM_LEG(Dir,NxtDir,PERMUTE_DIR,PROJ,RECON)

#define RESULT(base,basep) SAVE_RESULT(base,basep);

#endif

////////////////////////////////////////////////////////////////////////////////
// Pre comms kernel -- prefetch like normal because it is mostly right
////////////////////////////////////////////////////////////////////////////////
#ifdef INTERIOR

#define ASM_LEG(Dir,NxtDir,PERMUTE_DIR,PROJ,RECON)			\
      basep = st.GetPFInfo(nent,plocal); nent++;			\
      if ( local ) {							\
  LOAD_CHIMU(base);                                       \
  LOAD_TABLE(PERMUTE_DIR);                                \
  PROJ;							                        \
  MAYBEPERM(PERMUTE_DIR,perm);					        \
      }else if ( st.same_node[Dir] ) {LOAD_CHI(base);}			\
      if ( local || st.same_node[Dir] ) {				\
  MULT_2SPIN_1(Dir);					                    \
  MULT_2SPIN_2;					                        \
  RECON;								\
      }									\
  base = st.GetInfo(ptype,local,perm,NxtDir,ent,plocal); ent++;	\
  PREFETCH_CHIMU(base);						\
  PREFETCH_CHIMU_L2(basep);                               \

#define ASM_LEG_XP(Dir,NxtDir,PERMUTE_DIR,PROJ,RECON)			\
  base = st.GetInfo(ptype,local,perm,Dir,ent,plocal); ent++;		\
  PREFETCH1_CHIMU(base);						\
  { ZERO_PSI; }								\
  ASM_LEG(Dir,NxtDir,PERMUTE_DIR,PROJ,RECON)

#define RESULT(base,basep) SAVE_RESULT(base,basep);

#endif

////////////////////////////////////////////////////////////////////////////////
// Post comms kernel
////////////////////////////////////////////////////////////////////////////////
#ifdef EXTERIOR

#define ASM_LEG(Dir,NxtDir,PERMUTE_DIR,PROJ,RECON)			\
  base = st.GetInfo(ptype,local,perm,Dir,ent,plocal); ent++;		\
  if((!local)&&(!st.same_node[Dir]) ) {					\
    LOAD_CHI(base);							\
    MULT_2SPIN_1(Dir);					                    \
    MULT_2SPIN_2;					                        \
    RECON;								\
    nmu++;								\
  }

#define ASM_LEG_XP(Dir,NxtDir,PERMUTE_DIR,PROJ,RECON)			\
  nmu=0;								\
  { ZERO_PSI;}								\
  base = st.GetInfo(ptype,local,perm,Dir,ent,plocal); ent++;		\
  if((!local)&&(!st.same_node[Dir]) ) {					\
    LOAD_CHI(base);							\
    MULT_2SPIN_1(Dir);					                    \
    MULT_2SPIN_2;					                        \
    RECON;								\
    nmu++;								\
  }

#define RESULT(base,basep) if (nmu){ ADD_RESULT(base,base);}

#endif


{
  int nmu;
  int local,perm, ptype;
  uint64_t base;
  uint64_t basep;
  const uint64_t plocal =(uint64_t) & in[0];

  MASK_REGS;
  int nmax=U.oSites();
  for(int site=0;site<Ns;site++) {
#ifndef EXTERIOR
    //    int sU =lo.Reorder(ssU);
    int sU =ssU;
    int ssn=ssU+1;     if(ssn>=nmax) ssn=0;
    //    int sUn=lo.Reorder(ssn);
    int sUn=ssn;
#else
    int sU =ssU;
    int ssn=ssU+1;     if(ssn>=nmax) ssn=0;
    int sUn=ssn;
#endif
    for(int s=0;s<Ls;s++) {
      ss =sU*Ls+s;
      ssn=sUn*Ls+s;
      int  ent=ss*8;// 2*Ndim
      int nent=ssn*8;

      uint64_t delta_base, delta_base_p;

   ASM_LEG_XP(Xp,Yp,PERMUTE_DIR3,DIR0_PROJ,DIR0_RECON);

#ifdef SHOW
      float rescale = 64. * 12.;
      std::cout << "=================================================================" << std::endl;
      std::cout << "ss = " << ss << "   ssn = " << ssn << std::endl;
      std::cout << "sU = " << sU << "   ssU = " << ssU << std::endl;
      std::cout << " " << std::endl;


      std::cout << "Dir = " << Xp << "        "  << WHERE<< std::endl;

      std::cout << "ent  nent  local  perm       = " << ent << "  " << nent << "  " << local << "  "  << perm << std::endl;
      std::cout << "st.same_node[Dir] = " << st.same_node[Xp] << std::endl;
      std::cout << "base              = " << (base - plocal)/rescale << std::endl;
      std::cout << "Basep             = " << (basep - plocal)/rescale << std::endl;
      //printf("U                 = %llu\n", (uint64_t)&[sU](Dir));
      std::cout << "----------------------------------------------------" << std::endl;
#endif

      ASM_LEG(Yp,Zp,PERMUTE_DIR2,DIR1_PROJ,DIR1_RECON);

#ifdef SHOW
      std::cout << "Dir = " << Yp << "        "  << WHERE<< std::endl;

      std::cout << "ent  nent  local  perm       = " << ent << "  " << nent << "  " << local << "  "  << perm << std::endl;
      std::cout << "st.same_node[Dir] = " << st.same_node[Yp] << std::endl;
      std::cout << "base              = " << (base - plocal)/rescale << std::endl;
      std::cout << "Basep             = " << (basep - plocal)/rescale << std::endl;
      //printf("U                 = %llu\n", (uint64_t)&[sU](Dir));
      std::cout << "----------------------------------------------------" << std::endl;
#endif

      ASM_LEG(Zp,Tp,PERMUTE_DIR1,DIR2_PROJ,DIR2_RECON);

#ifdef SHOW
      std::cout << "Dir = " << Zp << "        "  << WHERE<< std::endl;

      std::cout << "ent  nent  local  perm       = " << ent << "  " << nent << "  " << local << "  "  << perm << std::endl;
      std::cout << "st.same_node[Dir] = " << st.same_node[Zp] << std::endl;
      std::cout << "base              = " << (base - plocal)/rescale << std::endl;
      std::cout << "Basep             = " << (basep - plocal)/rescale << std::endl;
      //printf("U                 = %llu\n", (uint64_t)&[sU](Dir));
      std::cout << "----------------------------------------------------" << std::endl;
#endif

      ASM_LEG(Tp,Xm,PERMUTE_DIR0,DIR3_PROJ,DIR3_RECON);

#ifdef SHOW
      std::cout << "Dir = " << Tp << "        "  << WHERE<< std::endl;

      std::cout << "ent  nent  local  perm       = " << ent << "  " << nent << "  " << local << "  "  << perm << std::endl;
      std::cout << "st.same_node[Dir] = " << st.same_node[Tp] << std::endl;
      std::cout << "base              = " << (base - plocal)/rescale << std::endl;
      std::cout << "Basep             = " << (basep - plocal)/rescale << std::endl;
      //printf("U                 = %llu\n", (uint64_t)&[sU](Dir));
      std::cout << "----------------------------------------------------" << std::endl;
#endif

      ASM_LEG(Xm,Ym,PERMUTE_DIR3,DIR4_PROJ,DIR4_RECON);

#ifdef SHOW
      std::cout << "Dir = " << Xm << "        "  << WHERE<< std::endl;

      std::cout << "ent  nent  local  perm       = " << ent << "  " << nent << "  " << local << "  "  << perm << std::endl;
      std::cout << "st.same_node[Dir] = " << st.same_node[Xm] << std::endl;
      std::cout << "base              = " << (base - plocal)/rescale << std::endl;
      std::cout << "Basep             = " << (basep - plocal)/rescale << std::endl;
      //printf("U                 = %llu\n", (uint64_t)&[sU](Dir));
      std::cout << "----------------------------------------------------" << std::endl;
#endif

      // DC ZVA test
      // { uint64_t basestore = (uint64_t)&out[ss];
      //   PREFETCH_RESULT_L2_STORE(basestore); }


      ASM_LEG(Ym,Zm,PERMUTE_DIR2,DIR5_PROJ,DIR5_RECON);

#ifdef SHOW
      std::cout << "Dir = " << Ym << "        "  << WHERE<< std::endl;

      std::cout << "ent  nent  local  perm       = " << ent << "  " << nent << "  " << local << "  "  << perm << std::endl;
      std::cout << "st.same_node[Dir] = " << st.same_node[Ym] << std::endl;
      std::cout << "base              = " << (base - plocal)/rescale << std::endl;
      std::cout << "Basep             = " << (basep - plocal)/rescale << std::endl;
      //printf("U                 = %llu\n", (uint64_t)&[sU](Dir));
      std::cout << "----------------------------------------------------" << std::endl;
#endif

      // DC ZVA test
      //{ uint64_t basestore = (uint64_t)&out[ss];
      //  PREFETCH_RESULT_L2_STORE(basestore); }


      ASM_LEG(Zm,Tm,PERMUTE_DIR1,DIR6_PROJ,DIR6_RECON);

#ifdef SHOW
      std::cout << "Dir = " << Zm << "        "  << WHERE<< std::endl;

      std::cout << "ent  nent  local  perm       = " << ent << "  " << nent << "  " << local << "  "  << perm << std::endl;
      std::cout << "st.same_node[Dir] = " << st.same_node[Zm] << std::endl;
      std::cout << "base              = " << (base - plocal)/rescale << std::endl;
      std::cout << "Basep             = " << (basep - plocal)/rescale << std::endl;
      //printf("U                 = %llu\n", (uint64_t)&[sU](Dir));
      std::cout << "----------------------------------------------------" << std::endl;
#endif

      // DC ZVA test
      //{ uint64_t basestore = (uint64_t)&out[ss];
      //  PREFETCH_RESULT_L2_STORE(basestore); }


      ASM_LEG(Tm,Xp,PERMUTE_DIR0,DIR7_PROJ,DIR7_RECON);

#ifdef SHOW
      std::cout << "Dir = " << Tm << "        "  << WHERE<< std::endl;

      std::cout << "ent  nent  local  perm       = " << ent << "  " << nent << "  " << local << "  "  << perm << std::endl;
      std::cout << "st.same_node[Dir] = " << st.same_node[Tm] << std::endl;
      std::cout << "base              = " << (base - plocal)/rescale << std::endl;
      std::cout << "Basep             = " << (basep - plocal)/rescale << std::endl;
      //printf("U                 = %llu\n", (uint64_t)&[sU](Dir));
      std::cout << "----------------------------------------------------" << std::endl;
#endif

#ifdef EXTERIOR
      if (nmu==0) break;
      //      if (nmu!=0) std::cout << "EXT "<<sU<<std::endl;
#endif
      base = (uint64_t) &out[ss];
      basep= st.GetPFInfo(nent,plocal); ent++;
      basep = (uint64_t) &out[ssn];
      //PREFETCH_RESULT_L1_STORE(base);
      RESULT(base,basep);

#ifdef SHOW
      std::cout << "Dir = FINAL        " <<  WHERE<< std::endl;;

      base_ss = base;
      std::cout << "base              = " << (base - (uint64_t) &out[0])/rescale << std::endl;
      std::cout << "Basep             = " << (basep - plocal)/rescale << std::endl;
      //printf("U                 = %llu\n", (uint64_t)&[sU](Dir));
      std::cout << "----------------------------------------------------" << std::endl;
#endif

    }
    ssU++;
    UNLOCK_GAUGE(0);
  }
}

#undef DIR0_PROJ
#undef DIR1_PROJ
#undef DIR2_PROJ
#undef DIR3_PROJ
#undef DIR4_PROJ
#undef DIR5_PROJ
#undef DIR6_PROJ
#undef DIR7_PROJ
#undef DIR0_RECON
#undef DIR1_RECON
#undef DIR2_RECON
#undef DIR3_RECON
#undef DIR4_RECON
#undef DIR5_RECON
#undef DIR6_RECON
#undef DIR7_RECON
#undef ASM_LEG
#undef ASM_LEG_XP
#undef RESULT
