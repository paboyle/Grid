#pragma once

NAMESPACE_BEGIN(Grid);

    template<class Field>
    void HighBoundCheck(LinearOperatorBase<Field> &HermOp, 
			Field &Phi,
			RealD hi)
    {
      // Eigenvalue bound check at high end
      PowerMethod<Field> power_method;
      auto lambda_max = power_method(HermOp,Phi);
      std::cout << GridLogMessage << "Pseudofermion action lamda_max "<<lambda_max<<"( bound "<<hi<<")"<<std::endl;
      assert( (lambda_max < hi) && " High Bounds Check on operator failed" );
    }
      
    template<class Field> void InverseSqrtBoundsCheck(int MaxIter,double tol,
						       LinearOperatorBase<Field> &HermOp,
						       Field &GaussNoise,
						       MultiShiftFunction &PowerNegHalf) 
    {
      GridBase *FermionGrid = GaussNoise.Grid();

      Field X(FermionGrid);
      Field Y(FermionGrid);
      Field Z(FermionGrid);

      X=GaussNoise;
      RealD Nx = norm2(X);

      ConjugateGradientMultiShift<Field> msCG(MaxIter,PowerNegHalf);
      msCG(HermOp,X,Y);
      msCG(HermOp,Y,Z);

      RealD Nz = norm2(Z);

      HermOp.HermOp(Z,Y);
      RealD Ny = norm2(Y);

      X=X-Y;
      RealD Nd = norm2(X);
      std::cout << "************************* "<<std::endl;
      std::cout << " noise                         = "<<Nx<<std::endl;
      std::cout << " (MdagM^-1/2)^2  noise         = "<<Nz<<std::endl;
      std::cout << " MdagM (MdagM^-1/2)^2  noise   = "<<Ny<<std::endl;
      std::cout << " noise - MdagM (MdagM^-1/2)^2  noise   = "<<Nd<<std::endl;
      std::cout << "************************* "<<std::endl;
      assert( (std::sqrt(Nd/Nx)<tol) && " InverseSqrtBoundsCheck ");
    }

    /* For a HermOp = M^dag M, check the approximation of  HermOp^{-1/inv_pow}
       by computing   |X -    HermOp * [ Hermop^{-1/inv_pow} ]^{inv_pow} X|  < tol  
       for noise X (aka GaussNoise).
       ApproxNegPow should be the rational approximation for   X^{-1/inv_pow}
    */
    template<class Field> void InversePowerBoundsCheck(int inv_pow,
						       int MaxIter,double tol,
						       LinearOperatorBase<Field> &HermOp,
						       Field &GaussNoise,
						       MultiShiftFunction &ApproxNegPow) 
    {
      GridBase *FermionGrid = GaussNoise.Grid();

      Field X(FermionGrid);
      Field Y(FermionGrid);
      Field Z(FermionGrid);

      Field tmp1(FermionGrid), tmp2(FermionGrid);

      X=GaussNoise;
      RealD Nx = norm2(X);

      ConjugateGradientMultiShift<Field> msCG(MaxIter,ApproxNegPow);

      tmp1 = X;
      
      Field* in = &tmp1;
      Field* out = &tmp2;
      for(int i=0;i<inv_pow;i++){ //apply  [ Hermop^{-1/inv_pow}  ]^{inv_pow} X =   HermOp^{-1} X
	msCG(HermOp, *in, *out); //backwards conventions!
	if(i!=inv_pow-1) std::swap(in, out);
      }
      Z = *out;

      RealD Nz = norm2(Z);

      HermOp.HermOp(Z,Y);
      RealD Ny = norm2(Y);

      X=X-Y;
      RealD Nd = norm2(X);
      std::cout << "************************* "<<std::endl;
      std::cout << " noise                         = "<<Nx<<std::endl;
      std::cout << " (MdagM^-1/" << inv_pow << ")^" << inv_pow << " noise         = "<<Nz<<std::endl;
      std::cout << " MdagM (MdagM^-1/" << inv_pow << ")^" << inv_pow << " noise   = "<<Ny<<std::endl;
      std::cout << " noise - MdagM (MdagM^-1/" << inv_pow << ")^" << inv_pow << " noise   = "<<Nd<<std::endl;
      std::cout << "************************* "<<std::endl;
      assert( (std::sqrt(Nd/Nx)<tol) && " InversePowerBoundsCheck ");
    }


NAMESPACE_END(Grid);

