    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/PartialFractionFermion5D.cc

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#include <Grid.h>
namespace Grid {
  namespace QCD {


    template<class Impl>
    void  PartialFractionFermion5D<Impl>::Mdir (const FermionField &psi, FermionField &chi,int dir,int disp){
      // this does both dag and undag but is trivial; make a common helper routing

      int sign = 1;
      int Ls = this->Ls;

      this->DhopDir(psi,chi,dir,disp);

      int nblock=(Ls-1)/2;
      for(int b=0;b<nblock;b++){
	int s = 2*b;
	ag5xpby_ssp(chi,-scale,chi,0.0,chi,s,s); 
	ag5xpby_ssp(chi, scale,chi,0.0,chi,s+1,s+1); 
      }
      ag5xpby_ssp(chi,p[nblock]*scale/amax,chi,0.0,chi,Ls-1,Ls-1);

    }
    template<class Impl>
    void   PartialFractionFermion5D<Impl>::Meooe_internal(const FermionField &psi, FermionField &chi,int dag)
    {
      int Ls = this->Ls;
      int sign = dag ? (-1) : 1;

      if ( psi.checkerboard == Odd ) {
	this->DhopEO(psi,chi,DaggerNo);
      } else {
	this->DhopOE(psi,chi,DaggerNo);
      }

      int nblock=(Ls-1)/2;
      for(int b=0;b<nblock;b++){
	int s = 2*b;
	ag5xpby_ssp(chi,-scale,chi,0.0,chi,s,s); 
	ag5xpby_ssp(chi, scale,chi,0.0,chi,s+1,s+1); 
      }
      ag5xpby_ssp(chi,p[nblock]*scale/amax,chi,0.0,chi,Ls-1,Ls-1);
    }

    template<class Impl>
    void   PartialFractionFermion5D<Impl>::Mooee_internal(const FermionField &psi, FermionField &chi,int dag)
    {
      // again dag and undag are trivially related
      int sign = dag ? (-1) : 1;
      int Ls = this->Ls;
      
      int nblock=(Ls-1)/2;
      for(int b=0;b<nblock;b++){
	
	int s = 2*b;
	RealD pp = p[nblock-1-b];
	RealD qq = q[nblock-1-b];
	
	// Do each 2x2 block aligned at s and multiplies Dw site diagonal by G5 so Hw
	ag5xpby_ssp(chi,-dw_diag*scale,psi,amax*sqrt(qq)*scale,psi, s  ,s+1); 
	ag5xpby_ssp(chi, dw_diag*scale,psi,amax*sqrt(qq)*scale,psi, s+1,s);
	axpby_ssp  (chi, 1.0, chi,sqrt(amax*pp)*scale*sign,psi,s+1,Ls-1);
      }
      
      {
	RealD R=(1+mass)/(1-mass);
	//R g5 psi[Ls-1] + p[0] H
	ag5xpbg5y_ssp(chi,R*scale,psi,p[nblock]*scale*dw_diag/amax,psi,Ls-1,Ls-1);
	
	for(int b=0;b<nblock;b++){
	  int s = 2*b+1;
	  RealD pp = p[nblock-1-b];
	  axpby_ssp(chi,1.0,chi,-sqrt(amax*pp)*scale*sign,psi,Ls-1,s);
	}
      }
    }

    template<class Impl>
    void   PartialFractionFermion5D<Impl>::MooeeInv_internal(const FermionField &psi, FermionField &chi,int dag)
    {
      int sign = dag ? (-1) : 1;
      int Ls = this->Ls;

      FermionField tmp(psi._grid);
      
      ///////////////////////////////////////////////////////////////////////////////////////
      //Linv
      ///////////////////////////////////////////////////////////////////////////////////////
      int nblock=(Ls-1)/2;

      axpy(chi,0.0,psi,psi); // Identity piece
      
      for(int b=0;b<nblock;b++){
	int s = 2*b;
	RealD pp = p[nblock-1-b];
	RealD qq = q[nblock-1-b];
	RealD coeff1=sign*sqrt(amax*amax*amax*pp*qq) / ( dw_diag*dw_diag + amax*amax* qq);
	RealD coeff2=sign*sqrt(amax*pp)*dw_diag / ( dw_diag*dw_diag + amax*amax* qq); // Implicit g5 here
	axpby_ssp  (chi,1.0,chi,coeff1,psi,Ls-1,s);
	axpbg5y_ssp(chi,1.0,chi,coeff2,psi,Ls-1,s+1);
      }
      
      ///////////////////////////////////////////////////////////////////////////////////////
      //Dinv (note D isn't really diagonal -- just diagonal enough that we can still invert)
      // Compute Seeinv (coeff of gamma5)
      ///////////////////////////////////////////////////////////////////////////////////////
      RealD R=(1+mass)/(1-mass);
      RealD Seeinv = R + p[nblock]*dw_diag/amax;
      for(int b=0;b<nblock;b++){
	Seeinv += p[nblock-1-b]*dw_diag/amax / ( dw_diag*dw_diag/amax/amax + q[nblock-1-b]);
      }    
      Seeinv = 1.0/Seeinv;
      
      for(int b=0;b<nblock;b++){
	int s = 2*b;
	RealD pp = p[nblock-1-b];
	RealD qq = q[nblock-1-b];
	RealD coeff1=dw_diag / ( dw_diag*dw_diag + amax*amax* qq); // Implicit g5 here
	RealD coeff2=amax*sqrt(qq) / ( dw_diag*dw_diag + amax*amax* qq);
	ag5xpby_ssp  (tmp,-coeff1,chi,coeff2,chi,s,s+1);
	ag5xpby_ssp  (tmp, coeff1,chi,coeff2,chi,s+1,s);
      }
      ag5xpby_ssp  (tmp, Seeinv,chi,0.0,chi,Ls-1,Ls-1);
      
      ///////////////////////////////////////////////////////////////////////////////////////
      // Uinv
      ///////////////////////////////////////////////////////////////////////////////////////
      for(int b=0;b<nblock;b++){
	int s = 2*b;
	RealD pp = p[nblock-1-b];
	RealD qq = q[nblock-1-b];
	RealD coeff1=-sign*sqrt(amax*amax*amax*pp*qq) / ( dw_diag*dw_diag + amax*amax* qq);
	RealD coeff2=-sign*sqrt(amax*pp)*dw_diag / ( dw_diag*dw_diag + amax*amax* qq); // Implicit g5 here
	axpby_ssp  (chi,1.0/scale,tmp,coeff1/scale,tmp,s,Ls-1);
	axpbg5y_ssp(chi,1.0/scale,tmp,coeff2/scale,tmp,s+1,Ls-1);
      }
      axpby_ssp  (chi, 1.0/scale,tmp,0.0,tmp,Ls-1,Ls-1);
    }

    template<class Impl>
    void   PartialFractionFermion5D<Impl>::M_internal(const FermionField &psi, FermionField &chi,int dag)
    {
      FermionField D(psi._grid);
  
      int Ls = this->Ls;
      int sign = dag ? (-1) : 1;

      // For partial frac Hw case (b5=c5=1) chroma quirkily computes
      //
      // Conventions for partfrac appear to be a mess.
      // Tony's Nara lectures have
      //
      // BlockDiag(  H/p_i  1             | 1       )    
      //          (  1      p_i H / q_i^2 | 0       )  
      //           ---------------------------------
      //           ( -1      0                | R  +p0 H  )
      //
      //Chroma     ( -2H    2sqrt(q_i)    |   0         )
      //           (2 sqrt(q_i)   2H      |  2 sqrt(p_i) )
      //           ---------------------------------
      //           ( 0     -2 sqrt(p_i)   |  2 R gamma_5 + p0 2H
      //
      // Edwards/Joo/Kennedy/Wenger
      //
      // Here, the "beta's" selected by chroma to scale the unphysical bulk constraint fields
      // incorporate the approx scale factor. This is obtained by propagating the
      // scale on "H" out to the off diagonal elements as follows:
      //
      // BlockDiag(  H/p_i  1             | 1       ) 
      //          (  1      p_i H / q_i^2 | 0       )  
      //           ---------------------------------
      //          ( -1      0                | R  + p_0 H  )
      //
      // becomes:
      // BlockDiag(  H/ sp_i  1               | 1             ) 
      //          (  1      sp_i H / s^2q_i^2 | 0             )  
      //           ---------------------------------
      //           ( -1      0                | R + p_0/s H   )
      //
      //
      // This is implemented in Chroma by
      //           p0' = p0/approxMax
      //           p_i' = p_i*approxMax
      //           q_i' = q_i*approxMax*approxMax
      //
      // After the equivalence transform is applied the matrix becomes
      // 
      //Chroma     ( -2H    sqrt(q'_i)    |   0         )
      //           (sqrt(q'_i)   2H       |   sqrt(p'_i) )
      //           ---------------------------------
      //           ( 0     -sqrt(p'_i)    |  2 R gamma_5 + p'0 2H
      //
      //     =     ( -2H    sqrt(q_i)amax    |   0              )
      //           (sqrt(q_i)amax   2H       |   sqrt(p_i*amax) )
      //           ---------------------------------
      //           ( 0     -sqrt(p_i)*amax   |  2 R gamma_5 + p0/amax 2H
      //

      this->DW(psi,D,DaggerNo); 

      int nblock=(Ls-1)/2;
      for(int b=0;b<nblock;b++){
	
	int s = 2*b;
	double pp = p[nblock-1-b];
	double qq = q[nblock-1-b];
	
	// Do each 2x2 block aligned at s and
	ag5xpby_ssp(chi,-1.0*scale,D,amax*sqrt(qq)*scale,psi, s  ,s+1); // Multiplies Dw by G5 so Hw
	ag5xpby_ssp(chi, 1.0*scale,D,amax*sqrt(qq)*scale,psi, s+1,s);
	
	// Pick up last column
	axpby_ssp  (chi, 1.0, chi,sqrt(amax*pp)*scale*sign,psi,s+1,Ls-1);
      }
	
      {
	double R=(1+this->mass)/(1-this->mass);
	//R g5 psi[Ls] + p[0] H
	ag5xpbg5y_ssp(chi,R*scale,psi,p[nblock]*scale/amax,D,Ls-1,Ls-1);
	
	for(int b=0;b<nblock;b++){
	  int s = 2*b+1;
	  double pp = p[nblock-1-b];
	  axpby_ssp(chi,1.0,chi,-sqrt(amax*pp)*scale*sign,psi,Ls-1,s);
	}
      }

    }

    template<class Impl>
    RealD  PartialFractionFermion5D<Impl>::M    (const FermionField &in, FermionField &out)
    {
      M_internal(in,out,DaggerNo);
      return norm2(out);
    }
    template<class Impl>
    RealD  PartialFractionFermion5D<Impl>::Mdag (const FermionField &in, FermionField &out)
    {
      M_internal(in,out,DaggerYes);
      return norm2(out);
    }

    template<class Impl>
    void PartialFractionFermion5D<Impl>::Meooe       (const FermionField &in, FermionField &out)
    {
      Meooe_internal(in,out,DaggerNo);
    }
    template<class Impl>
    void PartialFractionFermion5D<Impl>::MeooeDag    (const FermionField &in, FermionField &out)
    {
      Meooe_internal(in,out,DaggerYes);
    }
    template<class Impl>
    void PartialFractionFermion5D<Impl>::Mooee       (const FermionField &in, FermionField &out)
    {
      Mooee_internal(in,out,DaggerNo);
    }
    template<class Impl>
    void PartialFractionFermion5D<Impl>::MooeeDag    (const FermionField &in, FermionField &out)
    {
      Mooee_internal(in,out,DaggerYes);
    }

    template<class Impl>
    void PartialFractionFermion5D<Impl>::MooeeInv    (const FermionField &in, FermionField &out)
    {
      MooeeInv_internal(in,out,DaggerNo);
    }
    template<class Impl>
    void PartialFractionFermion5D<Impl>::MooeeInvDag (const FermionField &in, FermionField &out)
    {
      MooeeInv_internal(in,out,DaggerYes);
    }


  // force terms; five routines; default to Dhop on diagonal
    template<class Impl>
   void PartialFractionFermion5D<Impl>::MDeriv  (GaugeField &mat,const FermionField &U,const FermionField &V,int dag)
  {
    int Ls = this->Ls;

    FermionField D(V._grid);

    int nblock=(Ls-1)/2;
    for(int b=0;b<nblock;b++){
      int s = 2*b;
      ag5xpby_ssp(D,-scale,U,0.0,U,s,s); 
      ag5xpby_ssp(D, scale,U,0.0,U,s+1,s+1); 
    }
    ag5xpby_ssp(D,p[nblock]*scale/amax,U,0.0,U,Ls-1,Ls-1);

    this->DhopDeriv(mat,D,V,DaggerNo); 
  };
    template<class Impl>
   void PartialFractionFermion5D<Impl>::MoeDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag)
  {
    int Ls = this->Ls;

    FermionField D(V._grid);

    int nblock=(Ls-1)/2;
    for(int b=0;b<nblock;b++){
      int s = 2*b;
      ag5xpby_ssp(D,-scale,U,0.0,U,s,s); 
      ag5xpby_ssp(D, scale,U,0.0,U,s+1,s+1); 
    }
    ag5xpby_ssp(D,p[nblock]*scale/amax,U,0.0,U,Ls-1,Ls-1);

    this->DhopDerivOE(mat,D,V,DaggerNo); 
  };
    template<class Impl>
   void PartialFractionFermion5D<Impl>::MeoDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag)
  {
    int Ls = this->Ls;

    FermionField D(V._grid);

    int nblock=(Ls-1)/2;
    for(int b=0;b<nblock;b++){
      int s = 2*b;
      ag5xpby_ssp(D,-scale,U,0.0,U,s,s); 
      ag5xpby_ssp(D, scale,U,0.0,U,s+1,s+1); 
    }
    ag5xpby_ssp(D,p[nblock]*scale/amax,U,0.0,U,Ls-1,Ls-1);

    this->DhopDerivEO(mat,D,V,DaggerNo); 
  };

    template<class Impl>
    void  PartialFractionFermion5D<Impl>::SetCoefficientsTanh(Approx::zolotarev_data *zdata,RealD scale){
      SetCoefficientsZolotarev(1.0/scale,zdata);
    }
    template<class Impl>
    void  PartialFractionFermion5D<Impl>::SetCoefficientsZolotarev(RealD zolo_hi,Approx::zolotarev_data *zdata){

      // check on degree matching
      //      std::cout<<GridLogMessage << Ls << " Ls"<<std::endl;
      //      std::cout<<GridLogMessage << zdata->n  << " - n"<<std::endl;
      //      std::cout<<GridLogMessage << zdata->da << " -da "<<std::endl;
      //      std::cout<<GridLogMessage << zdata->db << " -db"<<std::endl;
      //      std::cout<<GridLogMessage << zdata->dn << " -dn"<<std::endl;
      //      std::cout<<GridLogMessage << zdata->dd << " -dd"<<std::endl;
      int Ls = this->Ls;

      assert(Ls == (2*zdata->da -1) );

      // Part frac
      //      RealD R;
      R=(1+mass)/(1-mass);
      dw_diag = (4.0-this->M5);

      //      std::vector<RealD> p; 
      //      std::vector<RealD> q;
      p.resize(zdata->da);
      q.resize(zdata->dd);
	
      for(int n=0;n<zdata->da;n++){
	p[n] = zdata -> alpha[n];
      }
      for(int n=0;n<zdata->dd;n++){
	q[n] = -zdata -> ap[n];
      }
      
      scale= part_frac_chroma_convention ? 2.0 : 1.0; // Chroma conventions annoy me

      amax=zolo_hi;
    }

      // Constructors
    template<class Impl>
    PartialFractionFermion5D<Impl>::PartialFractionFermion5D(GaugeField &_Umu,
							     GridCartesian         &FiveDimGrid,
							     GridRedBlackCartesian &FiveDimRedBlackGrid,
							     GridCartesian         &FourDimGrid,
							     GridRedBlackCartesian &FourDimRedBlackGrid,
							     RealD _mass,RealD M5,
							     const ImplParams &p) :
      WilsonFermion5D<Impl>(_Umu,
			    FiveDimGrid, FiveDimRedBlackGrid,
			    FourDimGrid, FourDimRedBlackGrid,M5,p),
      mass(_mass)

    {
      int Ls = this->Ls;

      assert((Ls&0x1)==1); // Odd Ls required
      int nrational=Ls-1;


      Approx::zolotarev_data *zdata = Approx::higham(1.0,nrational);

      // NB: chroma uses a cast to "float" for the zolotarev range(!?).
      // this creates a real difference in the operator which I do not like but we can replicate here
      // to demonstrate compatibility
      //      RealD eps = (zolo_lo / zolo_hi);
      //      zdata = bfm_zolotarev(eps,nrational,0);
      
      SetCoefficientsTanh(zdata,1.0);

      Approx::zolotarev_free(zdata);

    }
 
    FermOpTemplateInstantiate(PartialFractionFermion5D);

 }
}

