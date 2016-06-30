{
  int local,perm, ptype;
  uint64_t base;
  uint64_t basep;
  const uint64_t plocal =(uint64_t) & in._odata[0];

  //  vComplexF isigns[2] = { signs[0], signs[1] };
  vComplexF *isigns = &signs[0];

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
    XM_PROJMEM(base);
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
  XM_RECON;

  ////////////////////////////////
  // Yp
  ////////////////////////////////
  basep = st.GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    YM_PROJMEM(base);
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
  YM_RECON_ACCUM;

  ////////////////////////////////
  // Zp
  ////////////////////////////////
  basep = st.GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    ZM_PROJMEM(base);
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
  ZM_RECON_ACCUM;

  ////////////////////////////////
  // Tp
  ////////////////////////////////
  basep = st.GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    TM_PROJMEM(base);
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
  TM_RECON_ACCUM;

  ////////////////////////////////
  // Xm
  ////////////////////////////////
  basep= (uint64_t) &out._odata[ss];
  //  basep= st.GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    XP_PROJMEM(base);
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
  XP_RECON_ACCUM;

  ////////////////////////////////
  // Ym
  ////////////////////////////////
  basep= st.GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    YP_PROJMEM(base);
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
  YP_RECON_ACCUM;

  ////////////////////////////////
  // Zm
  ////////////////////////////////
  basep= st.GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    ZP_PROJMEM(base);
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
  ZP_RECON_ACCUM;

  ////////////////////////////////
  // Tm
  ////////////////////////////////
  basep= st.GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    TP_PROJMEM(base);
    MAYBEPERM(PERMUTE_DIR0,perm);
  } else { 
    LOAD_CHI(base);
  }
  base= (uint64_t) &out._odata[ss];
  PREFETCH_CHIMU(base);
  {
    MULT_2SPIN_DIR_PFTM(Tm,basep);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
  TP_RECON_ACCUM;

  basep= st.GetPFInfo(nent,plocal); nent++;
  SAVE_RESULT(base,basep);
  
  }
  ssU++;
  }
}
