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

    /*
Chroma: Recursive min norm
00132     Real dtau = traj_length / Real(n_steps);
00133     Real lambda_dt = dtau*lambda;
00134     Real dtauby2 = dtau / Real(2);
00135     Real one_minus_2lambda_dt = (Real(1)-Real(2)*lambda)*dtau;
00136     Real two_lambda_dt = lambda_dt*Real(2);
00137 
00138     // Its sts so:
00139     expSdt(s, lambda_dt); 
00140     for(int i=0; i < n_steps-1; i++) {  // N-1 full steps
00141       // Roll the exp(lambda_dt T) here and start
00142       // Next iter into one
00143       subIntegrator(s, dtauby2);         <--- either leapU or next integrator
00144       expSdt(s, one_minus_2lambda_dt);  
00145       subIntegrator(s, dtauby2);         <--- either leapU or next integrator
00146       expSdt(s, two_lambda_dt); 
00147     }
00148     // Last step, can't roll the first and last exp(lambda_dt T) 
00149     // together.
00150     subIntegrator(s, dtauby2);
00151     expSdt(s, one_minus_2lambda_dt);
00152     subIntegrator(s, dtauby2);
00153     expSdt(s, lambda_dt);

    * 
    */


    class MinimumNorm2{
      const double lambda = 0.1931833275037836;

    public:
      void step (LatticeLorentzColourMatrix& U, 
		 int level, std::vector<int>& clock,
		 Integrator<MinimumNorm2>* Integ){

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

	  if(clock[level] == 0){    // initial half step 
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
	  
	  // Handle the final half step
	  std::cout << GridLogMessage << " P[" <<level<<"] clock = "<<clock[level]<<"/"<<fin<<std::endl;
	  int mm = (clock[level]==fin) ? 1 : 2;
	  Integ->update_P(U,level,lambda*eps*mm); clock[level]+=mm;

	}
		
      }
      
    };



    /*

Chroma: Recursive leapfrog


00124     Real dtau = traj_length / n_steps;
00125     Real dtauby2 = dtau/2;
00126 
00127     // Its sts so:
00128     expSdt(s, dtauby2);  // First half step
00129     for(int i=0; i < n_steps-1; i++) {  // N-1 full steps
00130       subIntegrator(s, dtau);  <--- either leapU or next integrator
00131       expSdt(s, dtau);
00132     }
00133     subIntegrator(s, dtau);     // Last Full Step
00134     expSdt(s, dtauby2);  // Last Half Step

    * Nested 1:4; units in dt for top level integrator
    * CHROMA                           GUIDO
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
    *          P 1/16                    * skipped
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
    * Total
    */    
    class LeapFrog{
    public:
      void step (LatticeLorentzColourMatrix& U, 
		 int level, std::vector<int>& clock,
		 Integrator<LeapFrog>* Integ){

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

	  if ( level==0 ) {
	    first_step = (clock[level]==0);
	    last_step  = (clock[level]==fin);
	  } else {
	    first_step = (e==0);
	    last_step  = (e==multiplier-1);
	  }

	  if(first_step){    // initial half step
	    Integ->update_P(U, level,eps/2.0);	    ++clock[level];
	  }

	  if(level == fl){          // lowest level
	    Integ->update_U(U, eps);
	  }else{                 // recursive function call
	    step(U, level+1,clock, Integ);
	  }

	  // Handle the final half step
	  std::cout << GridLogMessage << " P[" <<level<<"] clock = "<<clock[level]<<"/"<<fin<<std::endl;
	  int mm = last_step ? 1 : 2;
	  Integ->update_P(U, level,mm*eps/2.0);	    
	  clock[level]+=mm;

	}
      }
    };


  }
}

#endif//INTEGRATOR_INCLUDED
