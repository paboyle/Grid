{
  int locala,perma, ptypea;
  int localb,permb, ptypeb;
  uint64_t basea, baseb;
  const uint64_t plocal =(uint64_t) & in._odata[0];

  //  vComplexF isigns[2] = { signs[0], signs[1] };
  vComplexF *isigns = &signs[0];

  MASK_REGS;

  for(int site=0;site<Ns;site++) {
  int sU=lo.Reorder(ssU);  
  for(int s=0;s<Ls;s++) {
  ss=sU*Ls+s;
  ////////////////////////////////
  // Xp
  ////////////////////////////////
  int ent=ss*8;// 2*Ndim
  basea = st.GetInfo(ptypea,locala,perma,Xp,ent,plocal); ent++;
  baseb = st.GetInfo(ptypeb,localb,permb,Yp,ent,plocal); ent++;

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
  basea = st.GetInfo(ptypea,locala,perma,Zp,ent,plocal); ent++;
  if ( localb ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    YM_PROJMEM(baseb);
    MAYBEPERM(PERMUTE_DIR2,permb);
  } else { 
    LOAD_CHI(baseb);
  }
  {
    MULT_2SPIN_DIR_PFYP(Yp,basea);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
  YM_RECON_ACCUM;

  ////////////////////////////////
  // Zp
  ////////////////////////////////
  baseb = st.GetInfo(ptypeb,localb,permb,Tp,ent,plocal); ent++;
  if ( locala ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    ZM_PROJMEM(basea);
    MAYBEPERM(PERMUTE_DIR1,perma);
  } else { 
    LOAD_CHI(basea);
  }
  {
    MULT_2SPIN_DIR_PFZP(Zp,baseb);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
  ZM_RECON_ACCUM;

  ////////////////////////////////
  // Tp
  ////////////////////////////////
  basea = st.GetInfo(ptypea,locala,perma,Xm,ent,plocal); ent++;
  if ( localb ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    TM_PROJMEM(baseb);
    MAYBEPERM(PERMUTE_DIR0,permb);
  } else { 
    LOAD_CHI(baseb);
  }
  {
    MULT_2SPIN_DIR_PFTP(Tp,basea);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
  TM_RECON_ACCUM;

  ////////////////////////////////
  // Xm
  ////////////////////////////////
  baseb = st.GetInfo(ptypeb,localb,permb,Ym,ent,plocal); ent++;
  if ( locala ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    XP_PROJMEM(basea);
    MAYBEPERM(PERMUTE_DIR3,perma);
  } else { 
    LOAD_CHI(basea);
  }
  {
    MULT_2SPIN_DIR_PFXM(Xm,baseb);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
  XP_RECON_ACCUM;

  ////////////////////////////////
  // Ym
  ////////////////////////////////
  basea = st.GetInfo(ptypea,locala,perma,Zm,ent,plocal); ent++;
  if ( localb ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    YP_PROJMEM(baseb);
    MAYBEPERM(PERMUTE_DIR2,permb);
  } else { 
    LOAD_CHI(baseb);
  }
  {
    MULT_2SPIN_DIR_PFYM(Ym,basea);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
  YP_RECON_ACCUM;

  ////////////////////////////////
  // Zm
  ////////////////////////////////
  baseb = st.GetInfo(ptypeb,localb,permb,Tm,ent,plocal); ent++;
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
  basea = (uint64_t)&out._odata[ss];
  if ( localb ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
    TP_PROJMEM(baseb);
    MAYBEPERM(PERMUTE_DIR0,permb);
  } else { 
    LOAD_CHI(baseb);
  }
  baseb = st.GetInfo(ptypeb,localb,permb,Xp,ent,plocal);
  {
    MULT_2SPIN_DIR_PFTM(Tm,basea);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
  TP_RECON_ACCUM;

  SAVE_RESULT(&out._odata[ss],baseb);

  } 
  ssU++;
  }
}
