/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/PrecGeneralisedConjugateResidual.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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
#ifndef GRID_PREC_GCR_NON_HERM_H
#define GRID_PREC_GCR_NON_HERM_H

///////////////////////////////////////////////////////////////////////////////////////////////////////
//VPGCR Abe and Zhang, 2005.
//INTERNATIONAL JOURNAL OF NUMERICAL ANALYSIS AND MODELING
//Computing and Information Volume 2, Number 2, Pages 147-161
//NB. Likely not original reference since they are focussing on a preconditioner variant.
//    but VPGCR was nicely written up in their paper
///////////////////////////////////////////////////////////////////////////////////////////////////////
NAMESPACE_BEGIN(Grid);

#define GCRLogLevel std::cout << GridLogMessage <<std::string(level,'\t')<< " Level "<<level<<" " 

template<class Field>
class PrecGeneralisedConjugateResidualNonHermitian : public LinearFunction<Field> {
public:                                                
  using LinearFunction<Field>::operator();
  RealD   Tolerance;
  Integer MaxIterations;
  int verbose;
  int mmax;
  int nstep;
  int steps;
  int level;
  GridStopWatch PrecTimer;
  GridStopWatch MatTimer;
  GridStopWatch LinalgTimer;

  LinearFunction<Field>     &Preconditioner;
  LinearOperatorBase<Field> &Linop;

  void Level(int lv) { level=lv; };

  PrecGeneralisedConjugateResidualNonHermitian(RealD tol,Integer maxit,LinearOperatorBase<Field> &_Linop,LinearFunction<Field> &Prec,int _mmax, int _nstep) : 
    Tolerance(tol), 
    MaxIterations(maxit),
    Linop(_Linop),
    Preconditioner(Prec),
    mmax(_mmax),
    nstep(_nstep)         // what is nstep vs mmax? one is the number of inner iterations
  { 
    level=1;
    verbose=1;
  };

  // virtual method stubs for updating GCR polynomial
  virtual void LogBegin(void){
    std::cout << "GCR::LogBegin() "<<std::endl;
  };
  virtual void LogIteration(int k, ComplexD a, std::vector<ComplexD> betas){
    std::cout << "GCR::LogIteration() "<<std::endl;
  };
  virtual void LogComplete(std::vector<ComplexD>& alphas, std::vector<std::vector<ComplexD>>& betas) {
    std::cout << "GCR::LogComplete() "<<std::endl;
  };

  void operator() (const Field &src, Field &psi){

    //    psi=Zero();
    RealD cp, ssq,rsq;
    ssq=norm2(src);
    rsq=Tolerance*Tolerance*ssq;
      
    Field r(src.Grid());

    PrecTimer.Reset();
    MatTimer.Reset();
    LinalgTimer.Reset();

    GridStopWatch SolverTimer;
    SolverTimer.Start();

    steps=0;
    for(int k=0;k<MaxIterations;k++){

      cp=GCRnStep(src,psi,rsq);

      GCRLogLevel <<"PGCR("<<mmax<<","<<nstep<<") "<< steps <<" steps cp = "<<cp<<" target "<<rsq <<std::endl;

      if(cp<rsq) {
        SolverTimer.Stop();

        Linop.Op(psi,r);
        axpy(r,-1.0,src,r);
        RealD tr = norm2(r);
        GCRLogLevel<<"PGCR: Converged on iteration " <<steps
          << " computed residual "<<sqrt(cp/ssq)
          << " true residual "    <<sqrt(tr/ssq)
          << " target "           <<Tolerance <<std::endl;

        GCRLogLevel<<"PGCR Time elapsed: Total  "<< SolverTimer.Elapsed() <<std::endl;
        return;
      }

    }
    GCRLogLevel<<"Variable Preconditioned GCR did not converge"<<std::endl;
    //    assert(0);
  }

  RealD GCRnStep(const Field &src, Field &psi,RealD rsq){

    RealD cp;
    ComplexD a, b;
    //    ComplexD zAz;
    RealD zAAz;
    ComplexD rq;

    GridBase *grid = src.Grid();

    Field r(grid);
    Field z(grid);
    Field tmp(grid);
    Field ttmp(grid);
    Field Az(grid);

    ////////////////////////////////
    // history for flexible orthog
    ////////////////////////////////
    std::vector<Field> q(mmax,grid);        // q = Ap
    std::vector<Field> p(mmax,grid);        // store mmax conjugate momenta
    std::vector<RealD> qq(mmax);            // qq = (Ap)^2 = <p|A^\dagger A |p> (denom of \alpha)
      
    GCRLogLevel<< "PGCR nStep("<<nstep<<")"<<std::endl;

    //////////////////////////////////
    // initial guess x0 is taken as nonzero.
    // r0=src-A x0 = src
    //////////////////////////////////
    MatTimer.Start();
    Linop.Op(psi,Az);
    //    zAz = innerProduct(Az,psi);
    zAAz= norm2(Az);
    MatTimer.Stop();
    

    LinalgTimer.Start();
    r=src-Az;
    LinalgTimer.Stop();
    GCRLogLevel<< "PGCR true residual r = src - A psi   "<< norm2(r) <<std::endl;

    this->LogBegin();       // initialize polynomial GCR if needed (TODO think about placement of this)
    
    /////////////////////
    // p = Prec(r)
    /////////////////////

    PrecTimer.Start();
    Preconditioner(r,z);
    PrecTimer.Stop();

    MatTimer.Start();
    Linop.Op(z,Az);
    MatTimer.Stop();

    LinalgTimer.Start();

    //    zAz = innerProduct(Az,psi);
    zAAz= norm2(Az);

    //p[0],q[0],qq[0] 
    p[0]= z;
    q[0]= Az;
    qq[0]= zAAz;

    std::cout << "||init p - src||: " << norm2(p[0] - src) << std::endl;   // for debugging
    
    cp =norm2(r);
    LinalgTimer.Stop();

    std::vector<ComplexD> all_alphas;
    std::vector<std::vector<ComplexD>> all_betas;

    for(int k=0;k<nstep;k++){

      steps++;

      int kp     = k+1;
      int peri_k = k %mmax;     // only store mmax vectors; just roll around if needed
      int peri_kp= kp%mmax;

      // std::cout << "peri_kp = " << peri_kp << std::endl;

      LinalgTimer.Start();
      rq= innerProduct(q[peri_k],r); // what if rAr not real?
      a = rq/qq[peri_k];              // compute alpha_j

      all_alphas.push_back(a);

      axpy(psi,a,p[peri_k],psi);      // update psi --> psi + \alpha p

      cp = axpy_norm(r,-a,q[peri_k],r);       // update r --> r - \alpha D p. Note q = Dp
      LinalgTimer.Stop();

      // LogIterationA(k + 1, a);

      GCRLogLevel<< "GCR step["<<steps<<"]  resid " << cp << " target " <<rsq<<std::endl; 

      // moving this to end of loop so that it doesn't exit beforehand
      // TODO if I want to uncomment this, I have to split the LogIteration again and put LogIterationA() beforehand
      // if((k==nstep-1)||(cp<rsq)){
      //   return cp;
      // }


      PrecTimer.Start();
      Preconditioner(r,z);// solve Az = r
      PrecTimer.Stop();

      MatTimer.Start();
      Linop.Op(z,Az);
      MatTimer.Stop();
      //      zAz = innerProduct(Az,psi);
      zAAz= norm2(Az);

      LinalgTimer.Start();

      q[peri_kp]=Az;
      p[peri_kp]=z;

      // Field Dsrc (grid);
      // Linop.Op(src, Dsrc);
      // std::cout << "||q[peri_kp] - D(src)||: " << norm2(q[peri_kp] - Dsrc) << std::endl;   // for debugging

          // // delete after testing
          // std::cout << "Testing Dsq on one for GCR: " << std::endl;
          // Field myField (grid);
          // myField = 1.0;
          // Field out1 (grid); Field out2 (grid);
          // Linop.HermOp(myField, out1);
          // Linop.Op(myField, out2);
          // std::cout << "Dsq.Hermop(ones) has norm " << norm2(out1) << std::endl;
          // std::cout << "Dsq.Op(ones) has norm " << norm2(out2) << std::endl;

      // basically northog = k+1 if mmax is large
      int northog = ((kp)>(mmax-1))?(mmax-1):(kp);  // if more than mmax done, we orthog all mmax history.
      // std::cout << "northog: " << northog << std::endl;
      std::vector<ComplexD> betas (northog);
      // std::cout << "peri_kp: " << peri_kp << std::endl;
      // we iterate backwards counting down from the current k+1 index (peri_kp) because we 
      for(int back=0;back<northog;back++){

        int peri_back=(k-back)%mmax;   	  assert((k-back)>=0);
        // std::cout << "peri_back: " << peri_back << std::endl;

        // b=-real(innerProduct(q[peri_back],Az))/qq[peri_back];
        b=-(innerProduct(q[peri_back],Az))/qq[peri_back];     // TODO try complex beta
        p[peri_kp]=p[peri_kp]+b*p[peri_back];
        q[peri_kp]=q[peri_kp]+b*q[peri_back];

        // LogIterationB(peri_back, b);
        // betas[back] = b;    // may need to change the indexing if I ever do it with restarts
        // std::cout << "[DEBUG] pushing beta for back = " << back << ", peri_back = " << peri_back << std::endl;

        betas[peri_back] = b;    // may need to change the indexing if I ever do it with restarts

      }
      qq[peri_kp]=norm2(q[peri_kp]); // could use axpy_norm
      LinalgTimer.Stop();

      // log iteration and update GCR polynomial if necessary.
      all_betas.push_back(betas);
      LogIteration(k + 1, a, betas);

      // finish if necessary
      if((k==nstep-1)||(cp<rsq)){
        std::cout << "All alphas: " << std::endl << all_alphas << std::endl;
        std::cout << "All betas: " << std::endl << all_betas << std::endl;
        LogComplete(all_alphas, all_betas);
        std::cout << "Exiting GCR." << std::endl;
        return cp;
      }

    }
    assert(0); // never reached
    return cp;
  }
};

class PolynomialFile: Serializable {
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PolynomialFile, 
      std::vector<std::vector<std::complex<double>>>, data,
      std::vector<std::vector<std::complex<double>>>, betas,
      std::vector<std::complex<double>>,              alphas
    );
};

// Optionally record the GCR polynomial. [PO]: TODO
template <class Field>
class PGCRPolynomial : public PrecGeneralisedConjugateResidualNonHermitian<Field> {
public:
  std::vector<ComplexD> ak;
  std::vector<std::vector<ComplexD>> bk;
  // std::vector<ComplexD> poly_p;
  std::vector<std::vector<ComplexD>> poly_p;
  std::vector<ComplexD> poly_Ap;        // polynomial in Ap_j (only store it for last p)
  std::vector<ComplexD> poly_r;
  std::vector<ComplexD> polynomial;

  PolynomialFile& PF;

public:
  PGCRPolynomial(RealD tol, Integer maxit,LinearOperatorBase<Field> &_Linop, LinearFunction<Field> &Prec, int _mmax, int _nstep, PolynomialFile& _PF)
    : PrecGeneralisedConjugateResidualNonHermitian<Field>(tol, maxit, _Linop, Prec, _mmax, _nstep), PF(_PF)
  {};

  // think this applies the polynomial in A = Linop to a field src. The coeffs are 
  // stored in the vector `polynomial`.
  void PolyOp(const Field &src, Field &psi)
  {
    Field tmp(src.Grid());
    Field AtoN(src.Grid());
    AtoN = src;
    psi=AtoN*polynomial[0];
    for(int n=1;n<polynomial.size();n++){
      tmp = AtoN;
      this->Linop.Op(tmp,AtoN);               // iterate A^n
      psi = psi + polynomial[n]*AtoN;       // psi += poly_n A^n src
    }
  }

  // [PO TODO] debug this
  void PGCRsequence(const Field &src, Field &x)
  {
    Field Ap(src.Grid());
    Field r(src.Grid());
    // Field p(src.Grid());
    // p=src;
    std::vector<Field> p;
    p.push_back(src);
    r=src;
    x=Zero();
    x.Checkerboard()=src.Checkerboard();
    for(int k=0;k<ak.size();k++){
      x = x + ak[k]*p[k];
      this->Linop.Op(p[k], Ap);
      r = r - ak[k] * Ap;
      // p[k] = r;
      p.push_back(r);
      for (int i = 0; i < k; i++) {     // [PO TODO] check indices
        p[k+1] += bk[i, k+1] * p[i];
      }
      // p = r + bk[k] * p;
    }
  }

  void Solve(const Field &src, Field &psi)
  {
    psi=Zero();
    this->operator()(src, psi);
  }

  virtual void LogBegin(void)
  {
    std::cout << "PGCR::LogBegin() "<<std::endl;
    ak.resize(0);
    bk.resize(0);
    polynomial.resize(0);
    poly_Ap.push_back(0.0);     // start with (0.0); during first iteration should change to (0.0, 1.0)
    std::vector<ComplexD> p0_tmp;
    p0_tmp.push_back(1.0);
    poly_p.push_back(p0_tmp);
    poly_r.push_back(1.0);
  };

  // Updates vector psi and r and initializes vector p[k+1]
  virtual void LogIteration(int k, ComplexD a, std::vector<ComplexD> betas){
    std::cout << "PGCR::LogIteration(k = " << k << ")" << std::endl;
    ak.push_back(a);
    bk.push_back(betas);

    // update Ap by pushing p[k] to the right
    poly_Ap.push_back(0.0);   // need to pad the end with an element
    poly_Ap[0] = 0.0;         // technically this should be unnecessary, as the first component is never set
    for(int i = 0; i < k; i++){
      poly_Ap[i+1]=poly_p[k-1][i];        // A\vec{p} = (0, \vec{p}) bc A shifts components of p to the right
    }

    // update psi_{k+1} --> psi_k + a_k p_k
    polynomial.push_back(0.0);
    for(int i = 0; i < k; i++) {
      polynomial[i] += a * poly_p[k-1][i];
    }
    PF.data.push_back(polynomial);

    //  r_{k+1} --> r_k - a_k A p_k
    //  p_{k+1} --> r_k + \sum_{i=0}^k \beta_{ik} p_i, input betas = (\beta_{ik})_i
    poly_r.push_back(0.0);        // should be of size k+1 if we start with k = 1
    std::vector<ComplexD> p_next (k + 1, ComplexD(0.0));     // p_{k+1} = same size as r_{k+1}
    for(int i = 0; i < k + 1; i++){
      poly_r[i] = poly_r[i] - a * poly_Ap[i];     // update r_{k+1} --> r_k - \alpha_k A p_k
      p_next[i] = poly_r[i];                 // init new vector as r_{k+1}
    }

    // p_{k+1} --> p_{k+1} + \sum_i \beta_{ij} p_i
    int nbeta = betas.size();
    std::cout << "Betas: " << betas << std::endl;
    for (int j = 0; j < nbeta; j++) {
      for (int i = 0; i < j+1; i++) {
        p_next[i] += betas[j] * poly_p[j][i];
      }
    }
    poly_p.push_back(p_next);                 // add p_{k+1} to the list of p's
  }

  virtual void LogComplete(std::vector<ComplexD>& alphas, std::vector<std::vector<ComplexD>>& betas) {
    /** Logs all alphas and betas to complete the iterations. */
    std::cout << "PGCR::LogComplete() "<<std::endl;
    for (int i = 0; i < alphas.size(); i++) {
      PF.alphas.push_back(alphas[i]);
      PF.betas.push_back(betas[i]);
    }
  };

};

NAMESPACE_END(Grid);
#endif
