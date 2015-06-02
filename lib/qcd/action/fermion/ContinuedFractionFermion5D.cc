#include <Grid.h>

namespace Grid {

  namespace QCD {

    RealD  ContinuedFractionFermion5D::M           (const LatticeFermion &psi, LatticeFermion &chi)
    {
      LatticeFermion D(psi._grid);

      DW(psi,D,DaggerNo); 

      int sign=1;
      for(int s=0;s<Ls;s++){
	if ( s==0 ) {
	  ag5xpby_ssp(chi,cc[0]*Beta[0]*sign*scale,D,sqrt_cc[0],psi,s,s+1); // Multiplies Dw by G5 so Hw
	} else if ( s==(Ls-1) ){
	  RealD R=(1.0+mass)/(1.0-mass);
	  ag5xpby_ssp(chi,Beta[s]*scale,D,sqrt_cc[s-1],psi,s,s-1);
	  ag5xpby_ssp(chi,R,psi,1.0,chi,s,s);
	} else {
	  ag5xpby_ssp(chi,cc[s]*Beta[s]*sign*scale,D,sqrt_cc[s],psi,s,s+1);
  	  axpby_ssp(chi,1.0,chi,sqrt_cc[s-1],psi,s,s-1);
	}
	sign=-sign; 
      }
      return norm2(chi);
    }
    RealD  ContinuedFractionFermion5D::Mdag        (const LatticeFermion &psi, LatticeFermion &chi)
    {
      // This matrix is already hermitian. (g5 Dw) = Dw dag g5 = (g5 Dw)dag
      // The rest of matrix is symmetric.
      // Can ignore "dag"
      return M(psi,chi);
    }
    void   ContinuedFractionFermion5D::Meooe       (const LatticeFermion &psi, LatticeFermion &chi)
    {
      Dhop(psi,chi,DaggerNo); // Dslash on diagonal. g5 Dslash is hermitian
      
      int sign=1;
      for(int s=0;s<Ls;s++){
	if ( s==(Ls-1) ){
	  ag5xpby_ssp(chi,Beta[s]*scale,chi,0.0,chi,s,s);
	} else {
	  ag5xpby_ssp(chi,cc[s]*Beta[s]*sign*scale,chi,0.0,chi,s,s);
	}
	sign=-sign; 
    }

    }
    void   ContinuedFractionFermion5D::MeooeDag    (const LatticeFermion &psi, LatticeFermion &chi)
    {
      Meooe(psi,chi);
    }
    void   ContinuedFractionFermion5D::Mooee       (const LatticeFermion &psi, LatticeFermion &chi)
    {
      double dw_diag = (4.0-this->M5)*scale;
    
      int sign=1;
      for(int s=0;s<Ls;s++){
	if ( s==0 ) {
	  ag5xpby_ssp(chi,cc[0]*Beta[0]*sign*dw_diag,psi,sqrt_cc[0],psi,s,s+1); // Multiplies Dw by G5 so Hw
	} else if ( s==(Ls-1) ){
	  // Drop the CC here.
	  double R=(1+this->mass)/(1-this->mass);
	  ag5xpby_ssp(chi,Beta[s]*dw_diag,psi,sqrt_cc[s-1],psi,s,s-1);
	  ag5xpby_ssp(chi,R,psi,1.0,chi,s,s);
	} else {
	  ag5xpby_ssp(chi,cc[s]*Beta[s]*sign*dw_diag,psi,sqrt_cc[s],psi,s,s+1);
	  axpby_ssp(chi,1.0,chi,sqrt_cc[s-1],psi,s,s-1);
	}
	sign=-sign; 
      }
    }

    void   ContinuedFractionFermion5D::MooeeDag    (const LatticeFermion &psi, LatticeFermion &chi)
    {
      Mooee(psi,chi);
    }
    void   ContinuedFractionFermion5D::MooeeInv    (const LatticeFermion &psi, LatticeFermion &chi)
    {
      // Apply Linv
      axpby_ssp(chi,1.0/cc_d[0],psi,0.0,psi,0,0);
      for(int s=1;s<Ls;s++){
	axpbg5y_ssp(chi,1.0/cc_d[s],psi,-1.0/See[s-1],chi,s,s-1);
      }
      // Apply Dinv
      for(int s=0;s<Ls;s++){
	ag5xpby_ssp(chi,1.0/See[s],chi,0.0,chi,s,s); //only appearance of See[0]
      }
      // Apply Uinv = (Linv)^T
      axpby_ssp(chi,1.0/cc_d[Ls-1],chi,0.0,chi,this->Ls-1,this->Ls-1);
      for(int s=Ls-2;s>=0;s--){
	axpbg5y_ssp(chi,1.0/cc_d[s],chi,-1.0*cc_d[s+1]/See[s]/cc_d[s],chi,s,s+1);
      }
    }
    void   ContinuedFractionFermion5D::MooeeInvDag (const LatticeFermion &psi, LatticeFermion &chi)
    {
      MooeeInv(psi,chi);
    }
    
    // Constructors
    ContinuedFractionFermion5D::ContinuedFractionFermion5D(
							   LatticeGaugeField &_Umu,
							   GridCartesian         &FiveDimGrid,
							   GridRedBlackCartesian &FiveDimRedBlackGrid,
							   GridCartesian         &FourDimGrid,
							   GridRedBlackCartesian &FourDimRedBlackGrid,
							   RealD _mass,RealD M5) :
      WilsonFermion5D(_Umu,
		      FiveDimGrid, FiveDimRedBlackGrid,
		      FourDimGrid, FourDimRedBlackGrid,M5),
      mass(_mass)
    {
    }

  }
}

