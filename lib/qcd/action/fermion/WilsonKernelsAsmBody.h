{
  int local,perm, ptype;
  uint64_t base;
  uint64_t basep;
  const uint64_t plocal =(uint64_t) & in._odata[0];

  //  vComplexF isigns[2] = { signs[0], signs[1] };
  //COMPLEX_TYPE is vComplexF of vComplexD depending 
  //on the chosen precision
  COMPLEX_TYPE *isigns = &signs[0];

  MASK_REGS;
  int nmax=U._grid->oSites();
  for(int site=0;site<Ns;site++) {
  int sU =lo.Reorder(ssU);
  int ssn=ssU+1; 
  if(ssn>=nmax) ssn=0;
  int sUn=lo.Reorder(ssn);
  for(int s=0;s<Ls;s++) {
  ss =sU*Ls+s;
  ssn=sUn*Ls+s; 
  ////////////////////////////////
  // Xp
  ////////////////////////////////
  int  ent=ss*8;// 2*Ndim
  int nent=ssn*8;

  PF_GAUGE(Xp); 
  base  = st.GetInfo(ptype,local,perm,Xp,ent,plocal); ent++;
  PREFETCH1_CHIMU(base);

  basep = st.GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);
#ifdef KERNEL_DAG
    XP_PROJMEM(base);
#else 
    XM_PROJMEM(base);
#endif
    MAYBEPERM(PERMUTE_DIR3,perm);
  } else { 
    LOAD_CHI(base);
  }
  base = st.GetInfo(ptype,local,perm,Yp,ent,plocal); ent++;
  PREFETCH_CHIMU(base);
  {
    MULT_2SPIN_DIR_PFXP(Xp,basep);
  }
  LOAD64(%r10,isigns);
#ifdef KERNEL_DAG
  XP_RECON;
#else
  XM_RECON;
#endif
  ////////////////////////////////
  // Yp
  ////////////////////////////////
  basep = st.GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
    YP_PROJMEM(base);
#else
    YM_PROJMEM(base);
#endif
    MAYBEPERM(PERMUTE_DIR2,perm);
  } else { 
    LOAD_CHI(base);
  }
  base  = st.GetInfo(ptype,local,perm,Zp,ent,plocal); ent++;
  PREFETCH_CHIMU(base);
  {
    MULT_2SPIN_DIR_PFYP(Yp,basep);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
  YP_RECON_ACCUM;
#else
  YM_RECON_ACCUM;
#endif

  ////////////////////////////////
  // Zp
  ////////////////////////////////
  basep = st.GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
    ZP_PROJMEM(base);
#else
    ZM_PROJMEM(base);
#endif
    MAYBEPERM(PERMUTE_DIR1,perm);
  } else { 
    LOAD_CHI(base);
  }
  base  = st.GetInfo(ptype,local,perm,Tp,ent,plocal); ent++;
  PREFETCH_CHIMU(base);
  {
    MULT_2SPIN_DIR_PFZP(Zp,basep);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
  ZP_RECON_ACCUM;
#else
  ZM_RECON_ACCUM;
#endif

  ////////////////////////////////
  // Tp
  ////////////////////////////////
  basep = st.GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
    TP_PROJMEM(base);
#else
    TM_PROJMEM(base);
#endif
    MAYBEPERM(PERMUTE_DIR0,perm);
  } else { 
    LOAD_CHI(base);
  }
  base = st.GetInfo(ptype,local,perm,Xm,ent,plocal); ent++;
  PREFETCH_CHIMU(base);
  {
    MULT_2SPIN_DIR_PFTP(Tp,basep);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
  TP_RECON_ACCUM;
#else
  TM_RECON_ACCUM;
#endif

  ////////////////////////////////
  // Xm
  ////////////////////////////////
#ifndef STREAM_STORE
  basep= (uint64_t) &out._odata[ss];
#endif
  //  basep= st.GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
    XM_PROJMEM(base);
#else
    XP_PROJMEM(base);
#endif
    MAYBEPERM(PERMUTE_DIR3,perm);
  } else { 
    LOAD_CHI(base);
  }
  base = st.GetInfo(ptype,local,perm,Ym,ent,plocal); ent++;
  PREFETCH_CHIMU(base);
  {
    MULT_2SPIN_DIR_PFXM(Xm,basep);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
  XM_RECON_ACCUM;
#else
  XP_RECON_ACCUM;
#endif

  ////////////////////////////////
  // Ym
  ////////////////////////////////
  basep= st.GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
    YM_PROJMEM(base);
#else
    YP_PROJMEM(base);
#endif
    MAYBEPERM(PERMUTE_DIR2,perm);
  } else { 
    LOAD_CHI(base);
  }
  base = st.GetInfo(ptype,local,perm,Zm,ent,plocal); ent++;
  PREFETCH_CHIMU(base);
  {
    MULT_2SPIN_DIR_PFYM(Ym,basep);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
  YM_RECON_ACCUM;
#else
  YP_RECON_ACCUM;
#endif

  ////////////////////////////////
  // Zm
  ////////////////////////////////
  basep= st.GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
    ZM_PROJMEM(base);
#else
    ZP_PROJMEM(base);
#endif
    MAYBEPERM(PERMUTE_DIR1,perm);
  } else { 
    LOAD_CHI(base);
  }
  base = st.GetInfo(ptype,local,perm,Tm,ent,plocal); ent++;
  PREFETCH_CHIMU(base);
  {
    MULT_2SPIN_DIR_PFZM(Zm,basep);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
  ZM_RECON_ACCUM;
#else
  ZP_RECON_ACCUM;
#endif

  ////////////////////////////////
  // Tm
  ////////////////////////////////
  basep= st.GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
    TM_PROJMEM(base);
#else
    TP_PROJMEM(base);
#endif
    MAYBEPERM(PERMUTE_DIR0,perm);
  } else { 
    LOAD_CHI(base);
  }
  base= (uint64_t) &out._odata[ss];
#ifndef STREAM_STORE
  PREFETCH_CHIMU(base);
#endif
  {
    MULT_2SPIN_DIR_PFTM(Tm,basep);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
  TM_RECON_ACCUM;
#else
  TP_RECON_ACCUM;
#endif

  basep= st.GetPFInfo(nent,plocal); nent++;
  SAVE_RESULT(base,basep);
  
  }
  ssU++;
  }
}
