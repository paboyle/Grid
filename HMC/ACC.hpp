namespace Grid{
  struct WFParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(WFParameters,
	    int, tau,
	    std::string, data_name,
	    std::string, path);
       

    template <class ReaderClass >
    WFParameters(Reader<ReaderClass>& Reader){
      read(Reader, "WilsonFlow", *this);
    }

  };

  struct ConfParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(ConfParameters,
           std::string, conf_prefix,
	   int, StartConfiguration,
	   int, EndConfiguration);
  
    template <class ReaderClass >
    ConfParameters(Reader<ReaderClass>& Reader){
      read(Reader, "Configurations", *this);
    }

  };

  struct ACFParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(ACFParameters,
	   int, MDtime_div_fac,				    
	   int, MScut,	          // Cutoff time for Madras-Sokal approx. to error in ACC
	   std::vector<int>, space_block_sizes,
	   int, R,                // Summation radius of Master field tecnnique
	   int, isFullTimeAvg);

    template <class ReaderClass >
    ACFParameters(Reader<ReaderClass>& Reader){
      read(Reader, "Autocorrelations", *this);
    }

  };
}

template <class T> void readFile(T& out, std::string const fname){
  Grid::emptyUserRecord record;
  Grid::ScidacReader RD;
  RD.open(fname);
  RD.readScidacFieldRecord(out,record);
  RD.close();
}

template <class T> void writeFile(T& in, std::string const fname){  
#ifdef HAVE_LIME
  // Ref: https://github.com/paboyle/Grid/blob/feature/scidac-wp1/tests/debug/Test_general_coarse_hdcg_phys48.cc#L111
  // the same function also in ActionBase.h
  std::cout << Grid::GridLogMessage << "Writes to: " << fname << std::endl;
  Grid::emptyUserRecord record;
  Grid::ScidacWriter WR(in.Grid()->IsBoss());
  WR.open(fname);
  WR.writeScidacFieldRecord(in,record,0);
  WR.close();
#endif
}

template <class L, typename A> void MS_approx(Grid::GridBase *Coarse, std::vector<L,A>  const& in, int T, int W, int block_size, int tau,
					      std::string const e_name, std::string ACC_Type ){
  using namespace Grid;
  // applicable only if blocking the src fields
  // Otherwise, for block averaged G_B(t), its error is hard to estimate as we need to consider correlation of G_x(t) at diff. sites wihting a block

  RealD sumG[in.size()], ACC;
  for (int t=0; t<T; t++) {
    sumG[t] = TensorRemove(sum(in[t])).real();
    std::cout << GridLogMessage << ACC_Type + " MFACC " + e_name + " (Madras-Sokal Approx): " << tau << " " << block_size << " " << 1 << " " << t << " "
              << sumG[t]/sumG[0] << std::endl;
    std::cout << GridLogMessage << ACC_Type + " MFCOV "  + e_name + " (Madras-Sokal Approx): " << tau << " " << block_size << " " << 1 << " " << t << " "
	      << sumG[t] << std::endl;
  }
  
  L var(Coarse), tmp(Coarse), Gt(Coarse), G0(Coarse);
  for (int t=0; 2*t<T-W; t++){ // \sigma_{\rho}(t) uses \rho(t') up to 2t+W
    ACC = sumG[t]/sumG[0];
    var = (1/sumG[t])*in[t] - (1/sumG[0])*in[0];
    var = (2.0/RealD(T))*var*var;
    for (int k=1; k<W+t; k++) {
      tmp = (1.0/sumG[t]/sumG[t]+2.0/sumG[0]/sumG[0])*in[k]*in[k] +
	in[std::abs(k-t)]*( (1.0/sumG[t]/sumG[t])*in[k+t] - (2.0/sumG[0]/sumG[t])*in[k]) -(2.0/sumG[0]/sumG[t])*in[k]*in[k+t];
      var = var + (2.0/RealD(T))*tmp;
    }
    std::cout << GridLogMessage << ACC_Type + " Variance " + e_name + " (Madras-Sokal Approx): " << tau << " " << block_size << " " << 1 << " " << t << " "
	      << TensorRemove(sum(var)).real()*ACC*ACC << std::endl;
  }
}

template <class L, typename A> void binning(Grid::GridBase *Coarse, std::vector<L,A>  const& in, int T, int n_bin, int block_size, int tau,
					    std::string const e_name, std::string ACC_Type){
  using namespace Grid;
  // error of G_b(t) is first computed via binning for each block b
  // G_b(t) with diff. b is considered as indep. measurement of G(t) => sigma_{G_(t)}^2 is reduced by V_b for the avg over blocks
  // Then, proceed to compute error of G(t)/G(0)
  // Computing ACC_b(t) = G_b(t)/G_b(0), 
  
  L avg(Coarse), avg0(Coarse), var(Coarse), var0(Coarse), cov(Coarse), tmp(Coarse);
  for (int t=0; t<T; t++){
    
    avg = Zero(); var = Zero(); cov = Zero();
    
    // Find avg field over binns
    for (int b=0; b<n_bin; b++)
      avg = avg + (1/RealD(n_bin))*in[b*T + t];
    if ( t== 0 ) avg0 = avg;
    RealD ACC = TensorRemove(sum(avg)).real()/TensorRemove(sum(avg0)).real();
    std::cout << GridLogMessage << ACC_Type + " MFACC " + e_name + " (Binning): " << tau << " " << block_size << " " << n_bin << " " << t << " "
              << ACC << std::endl;

    // Find variance of Gt/G0
    for	(int b=0; b<n_bin; b++){
      tmp = in[b*T + t] - avg;
      var = var + (1/RealD((n_bin-1)*n_bin))*(tmp*tmp);
      cov = cov + (1/RealD((n_bin-1)*n_bin))*tmp*(in[b*T] - avg0);
    }
    RealD G0 = real(sum(avg0)), Gt = real(sum(avg)); // the factor of 1/N_B is cancel out by 1/V^2 from Cov[G0,Gt] in the expression of var of Gt/G0
    if ( t == 0 ) var0 = var;
    std::cout << GridLogMessage << ACC_Type + " Variance " + e_name + " (Binning): " << tau << " " << block_size << " " << n_bin << " " << t << " "
              << (real(sum(var))/Gt/Gt + real(sum(var0))/G0/G0 - 2.0*real(sum(cov))/G0/Gt)*ACC*ACC 
	      << " " << (real(sum(var))/Gt/Gt + real(sum(var0))/G0/G0)*ACC*ACC << " " << - 2.0*real(sum(cov))/G0/Gt*ACC*ACC << std::endl;

  }
}

template <class L, typename A> void binning2(Grid::GridBase *Coarse, std::vector<L,A>  const& in, int T, int n_bin, int block_size, int tau,
					     std::string const e_name, std::string ACC_Type){
  using namespace Grid;
  // if bin_size > 30, the mean over binns approx dist. like Gaussian; the validity of this  estimation rests on Central Limit Theorem
  // we avg G_x(t) over blocks first and then take the ratio G(t)/G(0) before avg over bins
  // block avg might give a better sig as it incorpolates more stats; both should be unbiased
  
  RealD avg, avg0, var, tmp;
  std::vector<RealD> in_sum(in.size());

  for (int i=0; i<in.size(); i++) {
    in_sum[i] = TensorRemove(sum(in[i])).real();
    std::cout << GridLogMessage << ACC_Type + " MFCOV " + e_name + " (Binning2): " << tau << " " << block_size << " " << n_bin << " " << i << " "
	      << in_sum[i] << std::endl;
  }
  
  // Find avg field over binns
  for (int t=0; t<T; t++){
    avg = 0.0; var = 0.0;
    for (int b=0; b<n_bin; b++) 
      avg += (1/RealD(n_bin))*in_sum[b*T + t]/in_sum[b*T];
    std::cout << GridLogMessage << ACC_Type + " MFACC " + e_name + " (Binning2): " << tau << " " << block_size << " " << n_bin << " " << t << " "
              << avg << std::endl;

    // Find variance of Gt/G0 
    for (int b=0; b<n_bin; b++){
      tmp = in_sum[b*T + t]/in_sum[b*T] - avg;
      var = var + (1/RealD((n_bin-1)*n_bin))*(tmp*tmp);
    }
    std::cout << GridLogMessage << ACC_Type + " Variance " + e_name + " (Binning2): " << tau << " " << block_size << " " << n_bin << " " << t << " "
              << var << " " << (1-avg*avg)*(1-avg*avg)/RealD(Coarse->gSites()-3)*RealD(n_bin) << std::endl;
  }
}

template <class L, typename A> void MF_approx(Grid::GridBase *Coarse, std::vector<L,A>  const& in, int T, int R, int block_size, int tau,
                                              std::string const e_name, std::string ACC_Type ){
  using namespace Grid;

  // Assume: block averaging
  
  int R_b = R/block_size;
  RealD sumG[in.size()], ACC;
  for (int t=0; t<T; t++) {
    sumG[t] = TensorRemove(sum(in[t])).real();
    std::cout << GridLogMessage << ACC_Type + " MFACC " + e_name + " (Master-Field Approx): " << tau << " " << block_size << " " << R << " " << t << " "
              << sumG[t]/sumG[0] << std::endl;
  }

  RealD var, cov;
  L tmp(Coarse), Gt(Coarse), G0(Coarse), one(Coarse); one = typename L::scalar_type(1.0,0.0);
  G0 = in[0];
  G0 = G0 - (sumG[0]/RealD(Coarse->gSites())) * one;
  //std::cout << GridLogMessage << "in MF" << T <<  " " << block_size << " "<< R_b << " " << R << std::endl;
  for(int t=0; t<T; t++) {
    Gt = in[t] - (sumG[t]/RealD(Coarse->gSites())) * one;
    tmp = Gt;
    var = TensorRemove(sum(Gt*Gt)).real(); // C_tt
    cov = TensorRemove(sum(Gt*G0)).real(); // C_t0
    
    for (int r=1; r<=R_b; r++){
      // Shift Gt by y s.t. |y| <= r
      for(int x=-r; x<=r; x++)
	for(int y=-r; y<=r; y++)
	  for(int z=-r; z<=r; z++)
	    for(int s=-r; s<=r; s++) {
	      int r2_c = x*x+y*y+z*z+s*s;
	      if ( (r-1)*(r-1) < r2_c && r2_c <= r*r) {
		
		int d[4] = {x,y,z,s};
		Gt = tmp;
		for(int mu=0; mu<Nd; mu++) 
		  Gt = Cshift(Gt, mu, d[mu]);
		var += TensorRemove(sum(Gt*tmp)).real();
		cov += TensorRemove(sum(Gt*G0)).real();
	      }
	    }
      
      // Note: No division by the volume of the coarse lattice
      std::cout << GridLogMessage << ACC_Type + " Variance " + e_name + " (Master-Field Approx): " << tau << " " << block_size << " " << r << " " << t << " "
		<< sumG[t] << " " << var << " " << cov << std::endl;
    }
  }
}
