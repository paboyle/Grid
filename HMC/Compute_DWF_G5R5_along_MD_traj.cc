/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: HMC/analyze_snapshots.cc

Copyright (C) 2017

Author: Guido Cossu <guido.cossu@ed.ac.uk>
Author: Shuhei Yamamoto <syamamoto@bnl.gov>

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>
#include <string>

namespace Grid{
  struct SSParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(SSParameters,
				    std::string, fname_dt,
				    std::string, f_config,
				    int, i_config,
				    std::string, path,
				    std::string, save_path); 
       

    template <class ReaderClass >
    SSParameters(Reader<ReaderClass>& Reader){
      read(Reader, "SnapShots", *this);
    }

  };

  struct WFParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(WFParameters,
				    int, steps,
				    double, step_size,
				    int, meas_interval,
				    double, maxTau); // for the adaptive algorithm


    template <class ReaderClass >
    WFParameters(Reader<ReaderClass>& Reader){
      read(Reader, "WilsonFlow", *this);
    }

  };

  struct LanczosParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(LanczosParameters,
				    RealD, resid,
                                    RealD, mass ,
                                    RealD, M5 ,
                                    Integer, Ls,
                                    Integer, Nstop,
                                    Integer, Nk,
                                    Integer, Np,
				    Integer, Nsave,
                                    RealD, ChebyLow,
                                    RealD, ChebyHigh,
                                    Integer, ChebyOrder,
				    bool, recycle_src,
				    bool, use_prev_evec,
				    std::string, prev_t,
                                    std::string, fpath,
                                    std::string, fname,
                                    std::string, outpath)

    LanczosParameters() {
      // Default values
      resid             = 1e-5;
      mass              = 0.01;
      M5                = 1.8;
      Ls                = 16;
      Nk                = 20;
      Nstop             = Nk;
      Np                = 80;
      recycle_src       = true;
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

int main(int argc, char **argv) {
  using namespace Grid;
  
  Grid_init(&argc, &argv);
  GridLogLayout();

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  GridCartesian Grid(latt_size, simd_layout, mpi_layout);
  
  typedef Grid::XmlReader Serialiser;
  Serialiser        Reader("evecs.xml");//, false, "root");
  SSParameters      SSPar(Reader);
  WFParameters      WFPar(Reader);
  LanczosParameters LanPar;
  read(Reader,"LanczosParameters",LanPar);

  int Ls = LanPar.Ls;
  
  typedef DomainWallFermionD FermionOp;
  typedef typename DomainWallFermionD::FermionField FermionField;
 
  GridCartesian *         UGrid   = &Grid;
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);

  std::vector<int> seeds5({5, 6, 7, 8});
  GridParallelRNG RNG5(FGrid);     RNG5.SeedFixedIntegers(seeds5);
  
  
  std::vector<std::string> ts;
  if(Grid.IsBoss()){
    std::ifstream fin_t(SSPar.fname_dt);
    if (!fin_t) assert(false);
    {
      std::string t;
      while (std::getline(fin_t, t))
	ts.push_back(t);
    }
  }
  int vec_size = ts.size();
  MPI_Bcast(&vec_size, 1, MPI_INT, 0, Grid.communicator);
  for(int i=0; i<vec_size; i++){
    if(!Grid.IsBoss()) ts.push_back("");
    int size = ts[i].size();
    MPI_Bcast(&size, 1, MPI_INT, 0, Grid.communicator);
    if(!Grid.IsBoss()) ts[i].resize(size);
    MPI_Bcast((void *) ts[i].c_str(), size, MPI_CHAR, 0, Grid.communicator);
  }

  int num_meas = WFPar.steps/WFPar.meas_interval+1;
  LatticeGaugeField Ulat(&Grid), Usmr(&Grid), Uflow(&Grid);
  std::vector<FermionField> evec_prev_smr(num_meas,FGrid), evec_prev_lat(num_meas,FGrid); 
  FieldMetaData header;
  std::string flowed_field;
  if(LanPar.use_prev_evec){
    FermionField tmp(FGrid);
    std::string fname;
    std::string fname_pre  = LanPar.outpath + "/" + std::to_string(SSPar.i_config) + "/evec_";
    std::string fname_post = "_"+std::to_string(SSPar.i_config)+"_"+LanPar.prev_t;
    
    for(int i=0; i<num_meas; i++){
      std::string tau = (i==0)? "0" : "4";

      evec_prev_smr[i] = Zero();
      evec_prev_lat[i] = Zero();
      
      for(int j=0; j<LanPar.Nsave; j++){
	fname = fname_pre + std::to_string(j)+"_tau_"+tau+"_smr"+fname_post;
	readFile(tmp,fname);
	evec_prev_smr[i] = evec_prev_smr[i] + tmp;
	fname = fname_pre + std::to_string(j)+"_tau_"+tau+"_lat"+fname_post;
        readFile(tmp,fname);
	evec_prev_lat[i] = evec_prev_lat[i] + tmp;
      }
    }
  }
  else{
    for(int i=0; i<num_meas; i++){
      gaussian(RNG5, evec_prev_smr[i]);
      gaussian(RNG5, evec_prev_lat[i]);
    }
  }
  std::cout << std::setprecision(15);
  for (const std::string t : ts){
    /****************  Read Config   ************************/
    if(t=="0"){
      std::cout << GridLogMessage << SSPar.f_config+" t="+t <<std::endl;
      NerscIO::readConfiguration(Ulat,header,SSPar.f_config+"_lat."+std::to_string(SSPar.i_config));
      NerscIO::readConfiguration(Usmr,header,SSPar.f_config+"_lat_smr."+std::to_string(SSPar.i_config));
    }
    else {
      std::cout << GridLogMessage << SSPar.path+"U_lat_"+t <<std::endl;
      NerscIO::readConfiguration(Ulat,header,SSPar.path+"U_lat_"+t);
      NerscIO::readConfiguration(Usmr,header,SSPar.path+"U_smr_"+t);
    }
    
    /****************  Define Measurements on Flowed Field   ************************/
    WilsonFlow<PeriodicGimplR> WF(WFPar.step_size,WFPar.steps,WFPar.meas_interval);
    
    WF.addMeasurement(WFPar.meas_interval,
		      [&t,&flowed_field,&WFPar,&SSPar,&LanPar,&UGrid,&UrbGrid,&FGrid,&FrbGrid,&evec_prev_smr,&evec_prev_lat]
		      (int step, RealD tw, const typename PeriodicGimplR::GaugeField &constU){

      LatticeGaugeField U(UGrid); U = constU;
      std::string tau = std::to_string(int(std::round(tw)));

      std::string fpath  = SSPar.path;
      std::string spath  = SSPar.save_path;

      
      /************** Compute G5R5 Evecs on Flowed Fields    ****************************/

      RealD resid   = LanPar.resid;
      int   mass    = LanPar.mass;
      RealD M5      = LanPar.M5;
      int   Ls      = LanPar.Ls;
      int   Nk      = LanPar.Nk;
      int   Nstop   = LanPar.Nstop;
      int   Np      = LanPar.Np;
      int   Nsave   = LanPar.Nsave;
      int   Nm      = Nk + Np;
      int   MaxIt   = 10000;
    
      FermionOp                                                  Ddwf(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
      MdagMLinearOperator<FermionOp,FermionField>                HermOp(Ddwf);
      Gamma5R5HermitianLinearOperator<FermionOp, LatticeFermion> G5R5Herm(Ddwf);
      
      Chebyshev<FermionField>      Cheby   (LanPar.ChebyLow,LanPar.ChebyHigh,LanPar.ChebyOrder);
      FunctionHermOp<FermionField> OpCheby (Cheby,HermOp);
      PlainHermOp<FermionField>    Op      (HermOp);
      PlainHermOp<FermionField>    Op2     (G5R5Herm);
      
      ImplicitlyRestartedLanczos<FermionField> IRL(OpCheby, Op, Nstop, Nk, Nm, resid, MaxIt);
    
      /***********************************************************************/
      /*                    compute eigenvectors of D_H^2                    */
      /***********************************************************************/
      int Nconv;
      std::vector<RealD> eval(Nm);
      std::vector<FermionField> evec(Nm, FGrid);
      FermionField src(FGrid);
      int flw_meas_count = step/WFPar.meas_interval;
      std::cout << GridLogMessage << flowed_field << " t= " << t << " t_w= " << tau << " " << flw_meas_count << " " << step << " " << WFPar.meas_interval << std::endl;
      if(flowed_field == "lat")
	src = evec_prev_lat[flw_meas_count];
      else
	src = evec_prev_smr[flw_meas_count];
      IRL.calc(eval, evec, src, Nconv);

      std::cout << GridLogMessage << flowed_field << " t= " << t << " t_w= " << tau << " " <<  mass <<" : " << eval        << std::endl;
      std::cout << GridLogMessage << flowed_field << " t= " << t << " t_w= " << tau << " #evecs "   << evec.size() << std::endl;
      std::cout << GridLogMessage << flowed_field << " t= " << t << " t_w= " << tau << " Nconv  "   << Nconv       << std::endl;
      std::cout << GridLogMessage << flowed_field << " t= " << t << " t_w= " << tau << " Nm     "   << Nm          << std::endl;
      if ( Nconv > evec.size() ) Nconv = evec.size();
    
      /***********************************************************************/
      /*                       orthogonalization                             */
      /***********************************************************************/
      std::cout << GridLogMessage << flowed_field << " t= " << t << " t_w= " << tau << " Start orthogonalization "      << std::endl;
      // calculat the matrix
      std::cout << GridLogMessage << flowed_field << " t= " << t << " t_w= " << tau  << " calculate the matrix element" << std::endl;
      std::vector<LatticeFermion> G5R5Mevec(Nconv, FGrid);
      std::vector<LatticeFermion> finalevec(Nconv, FGrid);
      std::vector<RealD> eMe(Nconv), eMMe(Nconv);
      for(int i = 0; i < Nconv; i++){
	std::cout << GridLogMessage << flowed_field << " t= " << t << " t_w= " << tau << " calculate the matrix element["<<i<<"]" << std::endl;
	G5R5Herm.HermOpAndNorm(evec[i], G5R5Mevec[i], eMe[i], eMMe[i]);
      }
      std::cout << GridLogMessage << flowed_field << " t= " << t << " t_w= " << tau << " Re<evec, G5R5M(evec)>: "    << std::endl;
      std::cout << GridLogMessage << flowed_field << " t= " << t << " t_w= " << tau << " " << eMe                    << std::endl;
      std::cout << GridLogMessage << flowed_field << " t= " << t << " t_w= " << tau << " <G5R5M(evec), G5R5M(evec)>" << std::endl;
      std::cout << GridLogMessage << flowed_field << " t= " << t << " t_w= " << tau << " " << eMMe                   << std::endl;
      std::vector<std::vector<ComplexD>> VevecG5R5Mevec(Nconv);
      Eigen::MatrixXcd evecG5R5Mevec = Eigen::MatrixXcd::Zero(Nconv, Nconv);
      for(int i = 0; i < Nconv; i++){
	VevecG5R5Mevec[i].resize(Nconv); //can be just RealD tmp; or static conversion?
	for(int j = 0; j < Nconv; j++){
	  VevecG5R5Mevec[i][j] = innerProduct(evec[i], G5R5Mevec[j]);
	  evecG5R5Mevec(i, j) = VevecG5R5Mevec[i][j];
	}
      }
      // calculate eigenvector
      std::cout << GridLogMessage << flowed_field << " t= " << t << " t_w= " << tau << " Eigen solver" << std::endl;
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(evecG5R5Mevec);
      std::vector<RealD> eigeneval(Nconv);
      std::vector<std::vector<ComplexD>> eigenevec(Nconv);
      for(int i = 0; i < Nconv; i++){
	eigeneval[i] = eigensolver.eigenvalues()[i];
	eigenevec[i].resize(Nconv);
	for(int j = 0; j < Nconv; j++){
	  eigenevec[i][j] = eigensolver.eigenvectors()(i, j);
	}
      }
      //rotation
      std::cout << GridLogMessage << flowed_field << " t= " << t  << " t_w= " << tau << " Do rotation" << std::endl;
      for(int i = 0; i < Nconv; i++){
	finalevec[i] = finalevec[i] - finalevec[i];
	for(int j = 0; j < Nconv; j++){
	  finalevec[i] = eigenevec[j][i]*evec[j] + finalevec[i];
	}
      }
      //normalize again;
      for(int i = 0; i < Nconv; i++){
	RealD tmp_RealD = norm2(finalevec[i]);
	tmp_RealD = 1./pow(tmp_RealD, 0.5);
	finalevec[i] = finalevec[i]*tmp_RealD;
      }
      
      //check
      for(int i = 0; i < Nconv; i++){
	G5R5Herm.HermOpAndNorm(finalevec[i], G5R5Mevec[i], eMe[i], eMMe[i]);
      }
      
      /***********************************************************************/
      /*                    sort the eigenvectors                            */
      /***********************************************************************/
      std::vector<LatticeFermion> finalevec_copy(Nconv, FGrid);
      std::vector<RealD>          eMe_sort(eMe), eMe_sort2;
      for(int i = 0; i < Nconv; i++)
	finalevec_copy[i] = finalevec[i];

      /*******       Sort according to evals first          ******************/
      // Sort with sign first
      sort(eMe_sort.begin(), eMe_sort.end());
      // Fill eMe_sort2 with pos evals
      for(int i = 0; i < Nconv; i++)
	if(eMe_sort[i] >= 0)
	  eMe_sort2.push_back(eMe_sort[i]);
      // Sort in magnitude but + comes before - for pair of \pm evals of similar magnitude
      for(int i = 0; i < eMe_sort.size(); i++){
	int miss = 1;
        if(eMe_sort[i] < 0){
	  for(int j=0; j<eMe_sort2.size(); j++)
	    if(abs(eMe_sort2[j]) > abs(eMe_sort[i])){ // hits the wall
	      int pos = j;
	      if( j== 0 || abs(abs(eMe_sort2[j]) + eMe_sort[i]) < abs(abs(eMe_sort2[j-1]) + eMe_sort[i]) )
		pos += 1;
	      eMe_sort2.insert(eMe_sort2.begin()+pos,eMe_sort[i]);
	      miss = 0;
	      break;
	    }
	  // greater than any pos evals in magnitude
	  if(miss)eMe_sort2.push_back(eMe_sort[i]);
	}
      }

      for(int i = 0; i < Nconv; i++){
	for(int j = 0; j < Nconv; j++)
	  if(eMe[j] == eMe_sort2[i])
	    finalevec[i] = finalevec_copy[j];
      }
      if(LanPar.recycle_src){
	if(flowed_field == "lat"){
	  evec_prev_lat[flw_meas_count] = Zero();
	  for(auto evec:finalevec)
	    evec_prev_lat[flw_meas_count] += evec;
	}
	else{
	  evec_prev_smr[flw_meas_count] = Zero();
	  for(auto evec:finalevec)
	    evec_prev_smr[flw_meas_count] += evec;
	}
      }

      // check and save
      for(int i = 0; i < Nconv; i++){
	G5R5Herm.HermOpAndNorm(finalevec[i], G5R5Mevec[i], eMe[i], eMMe[i]);
	if(i < Nsave)
	  writeFile(finalevec[i],
		    LanPar.outpath + "/" + std::to_string(SSPar.i_config) + "/evec_" +
		    std::to_string(i)+"_tau_"+tau+"_"+flowed_field+"_"+std::to_string(SSPar.i_config)+"_"+t);
      }
      std::cout << GridLogMessage << flowed_field << " t= " << t << " t_w= " << tau << " Sorted Re<evec, G5R5M(evec)>: "    << std::endl;
      std::cout << GridLogMessage << flowed_field << " t= " << t << " t_w= " << tau << " Sorted Evals G5R5M " << eMe        << std::endl;
      std::cout << GridLogMessage << flowed_field << " t= " << t << " t_w= " << tau << " Sorted <G5R5M(evec), G5R5M(evec)>" << std::endl;
      std::cout << GridLogMessage << flowed_field << " t= " << t << " t_w= " << tau << " Sorted Evals G5R5M^2 "  << eMMe    << std::endl;
      
    
      /***********************************************************************/
      /*                   calculate chirality matrix                        */
      /***********************************************************************/

      std::vector<LatticeFermion>        G5evec(Nconv, FGrid);
      std::vector<std::vector<ComplexD>> chiral_matrix(Nconv);
      std::vector<std::vector<RealD>>    chiral_matrix_real(Nconv);
      for(int i = 0; i < Nconv; i++){
	G5evec[i] = Zero();
	for(int j = 0; j < Ls/2; j++){
	  axpby_ssp(G5evec[i], -1., finalevec[i], 0., G5evec[i], j, j);
	}
	for(int j = Ls/2; j < Ls; j++){
	  axpby_ssp(G5evec[i], 1., finalevec[i], 0., G5evec[i], j, j);
	}
      }
      // Compute spectral reconstruction of topological charge density
      std::string sp_file = LanPar.outpath + "/" + std::to_string(SSPar.i_config) + "/sp_sum_tau_"+tau+"_"+flowed_field+"_"+std::to_string(SSPar.i_config)+"_"+t;
      LatticeComplexD sp_sum(FGrid); sp_sum = Zero();
      for(int i = 0; i < Nconv; i++) {
	RealD sign = (eMe[i]>=0)? 1.0 : -1.0;
	RealD abs_lambda = sqrt(eMe[i]*eMe[i] - mass*mass);
	sp_sum = sp_sum - localInnerProduct(finalevec[i],G5evec[i]) + 0.5*sign*abs_lambda*localInnerProduct(finalevec[i],finalevec[i]);
      }
      LatticeComplexD sp_sum4D(UGrid), tmp_F(UGrid); sp_sum4D = Zero();
      for(int i=0; i<Ls;i++){
	ExtractSlice(tmp_F,sp_sum,i,0);
	sp_sum4D = sp_sum4D + tmp_F;
      }
      writeFile(sp_sum4D,sp_file);

      // Compute spectral reconstruction of noise part of topological charge density
      sp_file = LanPar.outpath + "/" + std::to_string(SSPar.i_config) + "/sp_noise_tau_"+tau+"_"+flowed_field+"_"+std::to_string(SSPar.i_config)+"_"+t;
      sp_sum = Zero();
      for(int i = 0; i < Nconv; i++) {
	if(eMe[i]<0){
	    RealD abs_lambda = sqrt(eMe[i]*eMe[i] - mass*mass);
	    // pos eval comes before neg eval
	    sp_sum = sp_sum - localInnerProduct(finalevec[i],G5evec[i]) - 0.5*abs_lambda*localInnerProduct(finalevec[i],finalevec[i]);
	    sp_sum = sp_sum - localInnerProduct(finalevec[i-1],G5evec[i-1]) + 0.5*abs_lambda*localInnerProduct(finalevec[i-1],finalevec[i-1]);
	  }
      }
      sp_sum4D = Zero();
      for(int i=0; i<Ls;i++){
        ExtractSlice(tmp_F,sp_sum,i,0);
        sp_sum4D = sp_sum4D + tmp_F;
      }
      writeFile(sp_sum4D,sp_file);
      
      for(int i = 0; i < Nconv; i++){
	chiral_matrix_real[i].resize(Nconv);
	chiral_matrix[i].resize(Nconv);
#if 0
	// Compute Evec Density and Save it
	auto evdensity = localInnerProduct(finalevec[i],finalevec[i] );
	writeFile(evdensity,
		  LanPar.outpath + "/" + std::to_string(SSPar.i_config) + "/evec_density" +
		  "_" + std::to_string(i)+"_tau_"+tau+"_"+flowed_field+"_"+std::to_string(SSPar.i_config)+"_"+t);
#endif
	// Compute Chiral Density and Save it
	for(int j = 0; j < Nconv; j++){
	  chiral_matrix[i][j] = innerProduct(finalevec[i], G5evec[j]);
	  chiral_matrix_real[i][j] = abs(chiral_matrix[i][j]);

	  std::cout << GridLogMessage << flowed_field << " t= " << t << " t_w= " << tau << " chiral_matrix_cplx "<<i<<" "<<j<<" "
		    << real(chiral_matrix[i][j]) << " " << imag(chiral_matrix[i][j]) << std::endl;
	  std::cout << GridLogMessage << flowed_field << " t= " << t << " t_w= " << tau << " chiral_matrix_real "<<i<<" "<<j<<" "<< chiral_matrix_real[i][j] << std::endl;
#if 0
	  if ( chiral_matrix_real[i][j] > 0.8 ) {
	    auto g5density = localInnerProduct(finalevec[i], G5evec[j]);
	    writeFile(g5density,
		      LanPar.outpath + "/" + std::to_string(SSPar.i_config) + "/chiral_density_" +
		      std::to_string(i)+"_"+std::to_string(j)+"_tau_"+tau+"_"+flowed_field+"_"+std::to_string(SSPar.i_config)+"_"+t);
	  }
#endif
	}
      }
      for(int i = 0; i < Nconv; i++){
	if(chiral_matrix[i][i].real() < 0.){
	  chiral_matrix_real[i][i] = -1. * chiral_matrix_real[i][i];
	}
      }
      
      // Save Chiral matrix for the config as a text file
      if( UGrid->IsBoss()){
	FILE *fp = fopen((LanPar.outpath + "/" + std::to_string(SSPar.i_config) +
			  "/chiral_matrix_real_"+"tau_"+tau+"_"+flowed_field+"_"+std::to_string(SSPar.i_config)+"_"+t).c_str(),"w");
	assert(fp!=NULL);
	for(int i = 0; i < Nconv; i++){
	  for(int j = 0; j < Nconv; j++){
	  fprintf(fp,"%lf ",chiral_matrix_real[i][j]);
	  }
	  fprintf(fp,"\n");
	}
	fclose(fp);
      }

    });
    
    /*********************   Perform Measurements  ***********************/
    flowed_field = "smr";
    WF.smear(Uflow, Usmr);
    flowed_field = "lat";
    WF.smear(Uflow, Ulat);
  }
  
  Grid_finalize();
}  // main


/*
Input file example


JSON

{
    "WilsonFlow":{
	"steps": 200,
	"step_size": 0.01,
	"meas_interval": 50,
  "maxTau": 2.0
    },
    "Configurations":{
	"conf_prefix": "ckpoint_lat",
	"rng_prefix": "ckpoint_rng",
	"StartConfiguration": 3000,
	"EndConfiguration": 3000,
	"Skip": 5
    }
}


*/
