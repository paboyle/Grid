#ifndef  GRID_QCD_DOMAIN_WALL_FERMION_H
#define  GRID_QCD_DOMAIN_WALL_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    class DomainWallFermion : public CayleyFermion5D
    {
    public:

      // Constructors
      DomainWallFermion(LatticeGaugeField &_Umu,
			GridCartesian         &FiveDimGrid,
			GridRedBlackCartesian &FiveDimRedBlackGrid,
			GridCartesian         &FourDimGrid,
			GridRedBlackCartesian &FourDimRedBlackGrid,
			RealD _mass,RealD _M5) : 

      CayleyFermion5D(_Umu,
		      FiveDimGrid,
		      FiveDimRedBlackGrid,
		      FourDimGrid,
		      FourDimRedBlackGrid,_mass,_M5)

      {
	RealD eps = 1.0;

	zdata = Approx::grid_higham(eps,this->Ls);// eps is ignored for higham
	assert(zdata->n==this->Ls);
 
	///////////////////////////////////////////////////////////
	// The Cayley coeffs (unprec)
	///////////////////////////////////////////////////////////
	this->omega.resize(this->Ls);
	this->bs.resize(this->Ls);
	this->cs.resize(this->Ls);
	this->as.resize(this->Ls);
	
	for(int i=0; i < this->Ls; i++){
	  this->as[i] = 1.0;
	  this->omega[i] = ((double)zdata -> gamma[i]);
	  double bb=1.0;
	  this->bs[i] = 0.5*(bb/(this->omega[i]) + 1.0);
	  this->cs[i] = 0.5*(bb/(this->omega[i]) - 1.0);
	}

	////////////////////////////////////////////////////////
	// Constants for the preconditioned matrix Cayley form
	////////////////////////////////////////////////////////
	this->bee.resize(this->Ls);
	this->cee.resize(this->Ls);
	this->beo.resize(this->Ls);
	this->ceo.resize(this->Ls);

	for(int i=0;i<this->Ls;i++){
	  this->bee[i]=as[i]*(bs[i]*(4.0-M5) +1.0);
	  this->cee[i]=as[i]*(1.0-cs[i]*(4.0-M5));
	  this->beo[i]=as[i]*bs[i];
	  this->ceo[i]=-as[i]*cs[i];
	}

	aee.resize(this->Ls);
	aeo.resize(this->Ls);
	for(int i=0;i<this->Ls;i++){
	  aee[i]=cee[i];
	  aeo[i]=ceo[i];
	}

	//////////////////////////////////////////
	// LDU decomposition of eeoo
	//////////////////////////////////////////
	dee.resize(this->Ls);
	lee.resize(this->Ls);
	leem.resize(this->Ls);
	uee.resize(this->Ls);
	ueem.resize(this->Ls);

	for(int i=0;i<this->Ls;i++){
	  
	  dee[i] = bee[i];
	  
	  if ( i < this->Ls-1 ) {

	    lee[i] =-cee[i+1]/bee[i]; // sub-diag entry on the ith column
	    
	    leem[i]=this->mass*cee[this->Ls-1]/bee[0];
	    for(int j=0;j<i;j++)  leem[i]*= aee[j]/bee[j+1];
	    
	    uee[i] =-aee[i]/bee[i];   // up-diag entry on the ith row
	    
	    ueem[i]=this->mass;
	    for(int j=1;j<=i;j++) ueem[i]*= cee[j]/bee[j];
	    ueem[i]*= aee[0]/bee[0];
	    
	  } else { 
	    lee[i] =0.0;
	    leem[i]=0.0;
	    uee[i] =0.0;
	    ueem[i]=0.0;
	  }
	}
	
	{ 
	  double delta_d=mass*cee[this->Ls-1];
	  for(int j=0;j<this->Ls-1;j++) delta_d *= cee[j]/bee[j];
	  dee[this->Ls-1] += delta_d;
	}
      }

    };

  }
}

#endif
