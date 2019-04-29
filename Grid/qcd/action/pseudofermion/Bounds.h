#pragma once

namespace Grid{
  namespace QCD{

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
      GridBase *FermionGrid = GaussNoise._grid;

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

  }
}
