//--------------------------------------------------------------------
/*! @file Integrator_algorithm.h
 * @brief Declaration of classes for the Molecular Dynamics algorithms
 *
 * @author Guido Cossu
 */
//--------------------------------------------------------------------

#ifndef INTEGRATOR_ALG_INCLUDED
#define INTEGRATOR_ALG_INCLUDED

namespace Grid{
  namespace QCD{

   /* PAB:
    *
    * Recursive leapfrog; explanation of nested stepping
    *
    * Nested 1:4; units in dt for top level integrator
    *
    * CHROMA                           IroIro
    *   0        1                      0              
    *  P 1/2                           P 1/2
    *          P 1/16                                  P1/16
    *                 U 1/8                                   U1/8
    *          P 1/8                                   P1/8
    *                 U 1/8                                   U1/8
    *          P 1/8                                   P1/8
    *                 U 1/8                                 U1/8
    *          P 1/8                                   P1/8
    *                 U 1/8                                 U1/8
    *          P 1/16                                  P1/8
    *  P 1                             P 1
    *          P 1/16                    * skipped --- avoids revaluating force
    *                 U 1/8                                 U1/8
    *          P 1/8                                   P1/8
    *                 U 1/8                                 U1/8
    *          P 1/8                                   P1/8
    *                 U 1/8                                 U1/8
    *          P 1/8                                   P1/8
    *                 U 1/8                                 U1/8
    *          P 1/16                                  P1/8
    *  P 1                             P 1
    *          P 1/16                    * skipped
    *                 U 1/8                                 U1/8
    *          P 1/8                                   P1/8
    *                 U 1/8                                 U1/8
    *          P 1/8                                   P1/8
    *                 U 1/8                                 U1/8
    *          P 1/8                                   P1/8
    *                 U 1/8                                 U1/8
    *          P 1/16                    * skipped
    *  P 1                             P 1
    *          P 1/16                                  P1/8
    *                 U 1/8                                 U1/8
    *          P 1/8                                   P1/8
    *                 U 1/8                                 U1/8
    *          P 1/8                                   P1/8
    *                 U 1/8                                 U1/8
    *          P 1/8                                   P1/8
    *                 U 1/8                                 U1/8
    *          P 1/16                                  P1/16
    *  P 1/2                            P 1/2
    */    

    template<class GaugeField> class LeapFrog {
    public:

      typedef LeapFrog<GaugeField> Algorithm;

      void step (GaugeField& U, 
		 int level, std::vector<int>& clock,
		 Integrator<GaugeField,Algorithm> * Integ){

	// level  : current level
	// fl     : final level
	// eps    : current step size
	
	int fl = Integ->as.size() -1;
	
	// Get current level step size
	int fin = 2*Integ->Params.MDsteps;
	for(int l=0; l<=level; ++l) fin*= Integ->as[l].multiplier;
	fin = fin-1;

	double eps = Integ->Params.stepsize;
	for(int l=0; l<=level; ++l) eps/= Integ->as[l].multiplier;
	
	int multiplier = Integ->as[level].multiplier;
	for(int e=0; e<multiplier; ++e){

	  int first_step,last_step;

	  first_step = (clock[level]==0);

	  if(first_step){    // initial half step
	    Integ->update_P(U, level,eps/2.0);	    ++clock[level];
	  }

	  if(level == fl){          // lowest level
	    Integ->update_U(U, eps);
	  }else{                 // recursive function call
	    step(U, level+1,clock, Integ);
	  }

	  last_step  = (clock[level]==fin);
	  int mm = last_step ? 1 : 2;
	  Integ->update_P(U, level,mm*eps/2.0);	    
	  clock[level]+=mm;

	}
      }
    };

    template<class GaugeField> class MinimumNorm2 {
    public:
      typedef MinimumNorm2<GaugeField> Algorithm;

    private:
      const double lambda = 0.1931833275037836;
    public:

      void step (GaugeField& U, 
		 int level, std::vector<int>& clock,
		 Integrator<GaugeField,Algorithm>* Integ){

	// level  : current level
	// fl     : final level
	// eps    : current step size

	int fl = Integ->as.size() -1;

	double eps = Integ->Params.stepsize;
	
	for(int l=0; l<=level; ++l) eps/= 2.0*Integ->as[l].multiplier;
	
	// which is final half step
	int fin = Integ->as[0].multiplier;
	for(int l=1; l<=level; ++l) fin*= 2.0*Integ->as[l].multiplier;
	fin = 3*Integ->Params.MDsteps*fin -1;

	int multiplier = Integ->as[level].multiplier;
	for(int e=0; e<multiplier; ++e){       // steps per step

	  int first_step,last_step;

	  first_step = (clock[level]==0);

	  if(first_step){    // initial half step 
	    Integ->update_P(U,level,lambda*eps);   ++clock[level];
	  }
	  
	  if(level == fl){          // lowest level 
	    Integ->update_U(U,0.5*eps);
	  }else{                 // recursive function call 
	    step(U,level+1,clock, Integ);
	  }
	  
	  Integ->update_P(U,level,(1.0-2.0*lambda)*eps); ++clock[level];
	  
	  if(level == fl){          // lowest level 
	    Integ->update_U(U,0.5*eps);
	  }else{                 // recursive function call 
	    step(U,level+1,clock, Integ);
	  }    
	  
	  last_step  = (clock[level]==fin);
	  int mm = (last_step) ? 1 : 2;
	  Integ->update_P(U,level,lambda*eps*mm); clock[level]+=mm;

	}
      }
    };

  }
}

#endif//INTEGRATOR_INCLUDED
