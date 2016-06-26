{
  int locala,perma, ptypea;
  int localb,permb, ptypeb;
  int localc,permc, ptypec;
  uint64_t basea, baseb, basec;

  const uint64_t plocal =(uint64_t) & in._odata[0];

  //  vComplexF isigns[2] = { signs[0], signs[1] };
  vComplexF *isigns = &signs[0];

  MASK_REGS;

  for(int site=0;site<Ns;site++) {
  int sU =lo.Reorder(ssU);
  for(int s=0;s<Ls;s++) {
  ss     =sU*Ls+s;

  ////////////////////////////////
  // Xp
  ////////////////////////////////
  int ent=ss*8;// 2*Ndim
  basea = st.GetInfo(ptypea,locala,perma,Xp,ent,plocal); ent++;
  PREFETCH1_CHIMU(basea);
  PF_GAUGE(Xp); 

  baseb = st.GetInfo(ptypeb,localb,permb,Yp,ent,plocal); ent++;
  PREFETCH_CHIMU(baseb);
  basec = st.GetInfo(ptypec,localc,permc,Zp,ent,plocal); ent++;
  PREFETCH_CHIMU(basec);

  if ( locala ) {
    LOAD64(%r10,isigns);
    XM_PROJMEM(basea);
    MAYBEPERM(PERMUTE_DIR3,perma);
  } else { 
    LOAD_CHI(basea);
  }
  {
    MULT_2SPIN_DIR_PFXP(Xp,baseb);
  }
  LOAD64(%r10,isigns);
  XM_RECON;

  ////////////////////////////////
  // Yp
  ////////////////////////////////
  basea = st.GetInfo(ptypea,locala,perma,Tp,ent,plocal); ent++;
  PREFETCH_CHIMU(basea);
  if ( localb ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    YM_PROJMEM(baseb);
    MAYBEPERM(PERMUTE_DIR2,permb);
  } else { 
    LOAD_CHI(baseb);
  }
  {
    MULT_2SPIN_DIR_PFYP(Yp,basec);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
  YM_RECON_ACCUM;

  ////////////////////////////////
  // Zp
  ////////////////////////////////
  baseb = st.GetInfo(ptypeb,localb,permb,Xm,ent,plocal); ent++;
  PREFETCH_CHIMU(baseb);
  if ( localc ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    ZM_PROJMEM(basec);
    MAYBEPERM(PERMUTE_DIR1,permc);
  } else { 
    LOAD_CHI(basec);
  }
  {
    MULT_2SPIN_DIR_PFZP(Zp,basea);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
  ZM_RECON_ACCUM;

  ////////////////////////////////
  // Tp
  ////////////////////////////////
  basec = st.GetInfo(ptypec,localc,permc,Ym,ent,plocal); ent++;
  PREFETCH_CHIMU(basec);
  if ( locala ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    TM_PROJMEM(basea);
    MAYBEPERM(PERMUTE_DIR0,perma);
  } else { 
    LOAD_CHI(basea);
  }
  {
    MULT_2SPIN_DIR_PFTP(Tp,baseb);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
  TM_RECON_ACCUM;

  ////////////////////////////////
  // Xm
  ////////////////////////////////
  basea = st.GetInfo(ptypea,locala,perma,Zm,ent,plocal); ent++;
  PREFETCH_CHIMU(basea);
  if ( localb ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    XP_PROJMEM(baseb);
    MAYBEPERM(PERMUTE_DIR3,permb);
  } else { 
    LOAD_CHI(baseb);
  }
  {
    MULT_2SPIN_DIR_PFXM(Xm,basec);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
  XP_RECON_ACCUM;

  ////////////////////////////////
  // Ym
  ////////////////////////////////
  baseb = st.GetInfo(ptypeb,localb,permb,Tm,ent,plocal); ent++;
  PREFETCH_CHIMU(baseb);
  if ( localc ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    YP_PROJMEM(basec);
    MAYBEPERM(PERMUTE_DIR2,permc);
  } else { 
    LOAD_CHI(basec);
  }
  {
    MULT_2SPIN_DIR_PFYM(Ym,basea);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
  YP_RECON_ACCUM;

  ////////////////////////////////
  // Zm
  ////////////////////////////////
  basec = (uint64_t)&out._odata[ss];
  PREFETCH_CHIMU(basec);
  if ( locala ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    ZP_PROJMEM(basea);
    MAYBEPERM(PERMUTE_DIR1,perma);
  } else { 
    LOAD_CHI(basea);
  }
  {
    MULT_2SPIN_DIR_PFZM(Zm,baseb);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
  ZP_RECON_ACCUM;

  ////////////////////////////////
  // Tm
  ////////////////////////////////
  //  basea = st.GetInfo(ptypea,locala,perma,Xp,ent,plocal); ent++;
  //  PREFETCH_CHIMU(basea);
  if ( localb ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    TP_PROJMEM(baseb);
    MAYBEPERM(PERMUTE_DIR0,permb);
  } else { 
    LOAD_CHI(baseb);
  }
  {
    MULT_2SPIN_DIR_PFTM(Tm,basec);
  }
  //  baseb = st.GetInfo(ptypeb,localb,permb,Yp,ent,plocal); ent++;
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
  TP_RECON_ACCUM;

  SAVE_RESULT(&out._odata[ss],basec);
  
  }
  ssU++;
  }
}
