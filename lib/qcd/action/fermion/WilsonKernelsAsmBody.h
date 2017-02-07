#ifdef KERNEL_DAG
#define DIR0_PROJMEM(base) XP_PROJMEM(base);
#define DIR1_PROJMEM(base) YP_PROJMEM(base);
#define DIR2_PROJMEM(base) ZP_PROJMEM(base);
#define DIR3_PROJMEM(base) TP_PROJMEM(base);
#define DIR4_PROJMEM(base) XM_PROJMEM(base);
#define DIR5_PROJMEM(base) YM_PROJMEM(base);
#define DIR6_PROJMEM(base) ZM_PROJMEM(base);
#define DIR7_PROJMEM(base) TM_PROJMEM(base);
#define DIR0_RECON   XP_RECON
#define DIR1_RECON   YP_RECON_ACCUM
#define DIR2_RECON   ZP_RECON_ACCUM
#define DIR3_RECON   TP_RECON_ACCUM
#define DIR4_RECON   XM_RECON_ACCUM
#define DIR5_RECON   YM_RECON_ACCUM
#define DIR6_RECON   ZM_RECON_ACCUM
#define DIR7_RECON   TM_RECON_ACCUM
#else
#define DIR0_PROJMEM(base) XM_PROJMEM(base);
#define DIR1_PROJMEM(base) YM_PROJMEM(base);
#define DIR2_PROJMEM(base) ZM_PROJMEM(base);
#define DIR3_PROJMEM(base) TM_PROJMEM(base);
#define DIR4_PROJMEM(base) XP_PROJMEM(base);
#define DIR5_PROJMEM(base) YP_PROJMEM(base);
#define DIR6_PROJMEM(base) ZP_PROJMEM(base);
#define DIR7_PROJMEM(base) TP_PROJMEM(base);
#define DIR0_RECON   XM_RECON
#define DIR1_RECON   YM_RECON_ACCUM
#define DIR2_RECON   ZM_RECON_ACCUM
#define DIR3_RECON   TM_RECON_ACCUM
#define DIR4_RECON   XP_RECON_ACCUM
#define DIR5_RECON   YP_RECON_ACCUM
#define DIR6_RECON   ZP_RECON_ACCUM
#define DIR7_RECON   TP_RECON_ACCUM
#endif

////////////////////////////////////////////////////////////////////////////////
// Comms then compute kernel
////////////////////////////////////////////////////////////////////////////////
#ifdef INTERIOR_AND_EXTERIOR

#define ZERO_NMU(A) 
#define INTERIOR_BLOCK_XP(a,b,PERMUTE_DIR,PROJMEM,RECON) INTERIOR_BLOCK(a,b,PERMUTE_DIR,PROJMEM,RECON)
#define EXTERIOR_BLOCK_XP(a,b,RECON) EXTERIOR_BLOCK(a,b,RECON)

#define INTERIOR_BLOCK(a,b,PERMUTE_DIR,PROJMEM,RECON)	\
  LOAD64(%r10,isigns);                                  \
  PROJMEM(base);                                        \
  MAYBEPERM(PERMUTE_DIR,perm);                                  

#define EXTERIOR_BLOCK(a,b,RECON)             \
  LOAD_CHI(base);

#define COMMON_BLOCK(a,b,RECON)               \
  base = st.GetInfo(ptype,local,perm,b,ent,plocal); ent++;     \
  PREFETCH_CHIMU(base);                                         \
  MULT_2SPIN_DIR_PF(a,basep);					\
  LOAD64(%r10,isigns);                                          \
  RECON;                                                        

#define RESULT(base,basep) SAVE_RESULT(base,basep);

#endif

////////////////////////////////////////////////////////////////////////////////
// Pre comms kernel -- prefetch like normal because it is mostly right
////////////////////////////////////////////////////////////////////////////////
#ifdef INTERIOR

#define COMMON_BLOCK(a,b,RECON)       
#define ZERO_NMU(A) 

// No accumulate for DIR0
#define EXTERIOR_BLOCK_XP(a,b,RECON)				\
  ZERO_PSI;							\
  base = st.GetInfo(ptype,local,perm,b,ent,plocal); ent++;	

#define EXTERIOR_BLOCK(a,b,RECON)  \
  base = st.GetInfo(ptype,local,perm,b,ent,plocal); ent++;     

#define INTERIOR_BLOCK_XP(a,b,PERMUTE_DIR,PROJMEM,RECON) INTERIOR_BLOCK(a,b,PERMUTE_DIR,PROJMEM,RECON)

#define INTERIOR_BLOCK(a,b,PERMUTE_DIR,PROJMEM,RECON)		\
  LOAD64(%r10,isigns);						\
  PROJMEM(base);                                                \
  MAYBEPERM(PERMUTE_DIR,perm);                                  \
  base = st.GetInfo(ptype,local,perm,b,ent,plocal); ent++;	\
  PREFETCH_CHIMU(base);						\
  MULT_2SPIN_DIR_PF(a,basep);					\
  LOAD64(%r10,isigns);                                          \
  RECON;                                                        

#define RESULT(base,basep) SAVE_RESULT(base,basep);

#endif

////////////////////////////////////////////////////////////////////////////////
// Post comms kernel
////////////////////////////////////////////////////////////////////////////////
#ifdef EXTERIOR

#define ZERO_NMU(A) nmu=0;

#define INTERIOR_BLOCK_XP(a,b,PERMUTE_DIR,PROJMEM,RECON) \
  ZERO_PSI;   base = st.GetInfo(ptype,local,perm,b,ent,plocal); ent++;		

#define EXTERIOR_BLOCK_XP(a,b,RECON) EXTERIOR_BLOCK(a,b,RECON)

#define INTERIOR_BLOCK(a,b,PERMUTE_DIR,PROJMEM,RECON)			\
  base = st.GetInfo(ptype,local,perm,b,ent,plocal); ent++;		

#define EXTERIOR_BLOCK(a,b,RECON)				\
    nmu++;							\
    LOAD_CHI(base);						\
    MULT_2SPIN_DIR_PF(a,base);					\
    base = st.GetInfo(ptype,local,perm,b,ent,plocal); ent++;	\
    LOAD64(%r10,isigns);					\
    RECON;                                                        

#define COMMON_BLOCK(a,b,RECON)			

#define RESULT(base,basep) if (nmu){  ADD_RESULT(base,base);}

#endif

{
  int nmu;
  int local,perm, ptype;
  uint64_t base;
  uint64_t basep;
  const uint64_t plocal =(uint64_t) & in._odata[0];

  COMPLEX_SIGNS(isigns);
  MASK_REGS;
  int nmax=U._grid->oSites();
  for(int site=0;site<Ns;site++) {
    int sU =lo.Reorder(ssU);
    int ssn=ssU+1;     if(ssn>=nmax) ssn=0;
    int sUn=lo.Reorder(ssn);
#ifndef EXTERIOR
    LOCK_GAUGE(0);
#endif
    for(int s=0;s<Ls;s++) {
      ss =sU*Ls+s;
      ssn=sUn*Ls+s; 
      int  ent=ss*8;// 2*Ndim
      int nent=ssn*8;

      ZERO_NMU(0);
      base  = st.GetInfo(ptype,local,perm,Xp,ent,plocal); ent++;
#ifndef EXTERIOR
      PF_GAUGE(Xp); 
      PREFETCH1_CHIMU(base);
#endif
      ////////////////////////////////
      // Xp
      ////////////////////////////////
      basep = st.GetPFInfo(nent,plocal); nent++;
      if ( local ) {
	INTERIOR_BLOCK_XP(Xp,Yp,PERMUTE_DIR3,DIR0_PROJMEM,DIR0_RECON);
      } else { 
	EXTERIOR_BLOCK_XP(Xp,Yp,DIR0_RECON);
      }
      COMMON_BLOCK(Xp,Yp,DIR0_RECON);
      ////////////////////////////////
      // Yp
      ////////////////////////////////
      basep = st.GetPFInfo(nent,plocal); nent++;
      if ( local ) {
	INTERIOR_BLOCK(Yp,Zp,PERMUTE_DIR2,DIR1_PROJMEM,DIR1_RECON);
      } else { 
	EXTERIOR_BLOCK(Yp,Zp,DIR1_RECON);
      }
      COMMON_BLOCK(Yp,Zp,DIR1_RECON);
      ////////////////////////////////
      // Zp
      ////////////////////////////////
      basep = st.GetPFInfo(nent,plocal); nent++;
      if ( local ) {
	INTERIOR_BLOCK(Zp,Tp,PERMUTE_DIR1,DIR2_PROJMEM,DIR2_RECON);
      } else { 
	EXTERIOR_BLOCK(Zp,Tp,DIR2_RECON);
      }
      COMMON_BLOCK(Zp,Tp,DIR2_RECON);
      ////////////////////////////////
      // Tp
      ////////////////////////////////
      basep = st.GetPFInfo(nent,plocal); nent++;
      if ( local ) {
	INTERIOR_BLOCK(Tp,Xm,PERMUTE_DIR0,DIR3_PROJMEM,DIR3_RECON);
      } else { 
	EXTERIOR_BLOCK(Tp,Xm,DIR3_RECON);
      }
      COMMON_BLOCK(Tp,Xm,DIR3_RECON);
      ////////////////////////////////
      // Xm
      ////////////////////////////////
      //  basep= st.GetPFInfo(nent,plocal); nent++;
      if ( local ) {
	INTERIOR_BLOCK(Xm,Ym,PERMUTE_DIR3,DIR4_PROJMEM,DIR4_RECON);
      } else { 
	EXTERIOR_BLOCK(Xm,Ym,DIR4_RECON);
      }
      COMMON_BLOCK(Xm,Ym,DIR4_RECON);
      ////////////////////////////////
      // Ym
      ////////////////////////////////
      basep= st.GetPFInfo(nent,plocal); nent++;
      if ( local ) {
	INTERIOR_BLOCK(Ym,Zm,PERMUTE_DIR2,DIR5_PROJMEM,DIR5_RECON);
      } else { 
	EXTERIOR_BLOCK(Ym,Zm,DIR5_RECON);
      }
      COMMON_BLOCK(Ym,Zm,DIR5_RECON);
      ////////////////////////////////
      // Zm
      ////////////////////////////////
      basep= st.GetPFInfo(nent,plocal); nent++;
      if ( local ) {
	INTERIOR_BLOCK(Zm,Tm,PERMUTE_DIR1,DIR6_PROJMEM,DIR6_RECON);
      } else { 
	EXTERIOR_BLOCK(Zm,Tm,DIR6_RECON);
      }
      COMMON_BLOCK(Zm,Tm,DIR6_RECON);
      ////////////////////////////////
      // Tm
      ////////////////////////////////
      basep= st.GetPFInfo(nent,plocal); nent++;
      if ( local ) {
	INTERIOR_BLOCK(Tm,Xp,PERMUTE_DIR0,DIR7_PROJMEM,DIR7_RECON);
      } else { 
	EXTERIOR_BLOCK(Tm,Xp,DIR7_RECON);
      }
      COMMON_BLOCK(Tm,Xp,DIR7_RECON);

      base = (uint64_t) &out._odata[ss];
      basep= st.GetPFInfo(nent,plocal); nent++;
      RESULT(base,basep);
    }
    ssU++;
    UNLOCK_GAUGE(0);
  }
}

#undef DIR0_PROJMEM
#undef DIR1_PROJMEM
#undef DIR2_PROJMEM
#undef DIR3_PROJMEM
#undef DIR4_PROJMEM
#undef DIR5_PROJMEM
#undef DIR6_PROJMEM
#undef DIR7_PROJMEM
#undef DIR0_RECON
#undef DIR1_RECON
#undef DIR2_RECON
#undef DIR3_RECON
#undef DIR4_RECON
#undef DIR5_RECON
#undef DIR6_RECON
#undef DIR7_RECON
#undef EXTERIOR_BLOCK
#undef INTERIOR_BLOCK
#undef EXTERIOR_BLOCK_XP
#undef INTERIOR_BLOCK_XP
#undef COMMON_BLOCK
#undef ZERO_NMU
#undef RESULT
