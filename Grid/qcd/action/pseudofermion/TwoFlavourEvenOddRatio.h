    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/pseudofermion/TwoFlavourEvenOddRatio.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#ifndef QCD_PSEUDOFERMION_TWO_FLAVOUR_EVEN_ODD_RATIO_H
#define QCD_PSEUDOFERMION_TWO_FLAVOUR_EVEN_ODD_RATIO_H

NAMESPACE_BEGIN(Grid);

    ///////////////////////////////////////
    // Two flavour ratio
    ///////////////////////////////////////
    template<class Impl>
    class TwoFlavourEvenOddRatioPseudoFermionAction : public Action<typename Impl::GaugeField> {
    public:
      INHERIT_IMPL_TYPES(Impl);
      
    private:
      FermionOperator<Impl> & NumOp;// the basic operator
      FermionOperator<Impl> & DenOp;// the basic operator

      OperatorFunction<FermionField> &DerivativeSolver;
      OperatorFunction<FermionField> &ActionSolver;
      OperatorFunction<FermionField> &HeatbathSolver;

      FermionField PhiOdd;   // the pseudo fermion field for this trajectory
      FermionField PhiEven;  // the pseudo fermion field for this trajectory

      RealD RefreshAction;
      
    public:
      TwoFlavourEvenOddRatioPseudoFermionAction(FermionOperator<Impl>  &_NumOp, 
                                                FermionOperator<Impl>  &_DenOp, 
                                                OperatorFunction<FermionField> & DS,
                                                OperatorFunction<FermionField> & AS ) : 
      TwoFlavourEvenOddRatioPseudoFermionAction(_NumOp,_DenOp, DS,AS,AS) {};

      TwoFlavourEvenOddRatioPseudoFermionAction(FermionOperator<Impl>  &_NumOp, 
                                                FermionOperator<Impl>  &_DenOp, 
                                                OperatorFunction<FermionField> & DS,
                                                OperatorFunction<FermionField> & AS, OperatorFunction<FermionField> & HS) :
      NumOp(_NumOp), 
      DenOp(_DenOp), 
      DerivativeSolver(DS), 
      ActionSolver(AS),
      HeatbathSolver(HS),
      PhiEven(_NumOp.FermionRedBlackGrid()),
      PhiOdd(_NumOp.FermionRedBlackGrid()) 
        {
          conformable(_NumOp.FermionGrid(), _DenOp.FermionGrid());
          conformable(_NumOp.FermionRedBlackGrid(), _DenOp.FermionRedBlackGrid());
          conformable(_NumOp.GaugeGrid(), _DenOp.GaugeGrid());
          conformable(_NumOp.GaugeRedBlackGrid(), _DenOp.GaugeRedBlackGrid());
        };

      virtual std::string action_name(){
	std::stringstream sstream;
	sstream<<"TwoFlavourEvenOddRatioPseudoFermionAction det("<<DenOp.Mass()<<") / det("<<NumOp.Mass()<<")";
	return sstream.str();
      }

      virtual std::string LogParameters(){
	std::stringstream sstream;
	sstream<< GridLogMessage << "["<<action_name()<<"] -- No further parameters "<<std::endl;
	return sstream.str();
      } 

      
      const FermionField &getPhiOdd() const{ return PhiOdd; }

      virtual void refresh(const GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG& pRNG) {
        // P(eta_o) = e^{- eta_o^dag eta_o}
        //
        // e^{x^2/2 sig^2} => sig^2 = 0.5.
        // 
        RealD scale = std::sqrt(0.5);

        FermionField eta    (NumOp.FermionGrid());
        gaussian(pRNG,eta); eta = eta * scale;

	refresh(U,eta);
      }

      void refresh(const GaugeField &U, const FermionField &eta) {

        // P(phi) = e^{- phi^dag Vpc (MpcdagMpc)^-1 Vpcdag phi}
        //
        // NumOp == V
        // DenOp == M
        //
        FermionField etaOdd (NumOp.FermionRedBlackGrid());
        FermionField etaEven(NumOp.FermionRedBlackGrid());
        FermionField tmp    (NumOp.FermionRedBlackGrid());

        pickCheckerboard(Even,etaEven,eta);
        pickCheckerboard(Odd,etaOdd,eta);

        NumOp.ImportGauge(U);
        DenOp.ImportGauge(U);
	std::cout << " TwoFlavourRefresh:  Imported gauge "<<std::endl;

        SchurDifferentiableOperator<Impl> Mpc(DenOp);
        SchurDifferentiableOperator<Impl> Vpc(NumOp);

	std::cout << " TwoFlavourRefresh: Diff ops "<<std::endl;
        // Odd det factors
        Mpc.MpcDag(etaOdd,PhiOdd);
	std::cout << " TwoFlavourRefresh: MpcDag "<<std::endl;
        tmp=Zero();
	std::cout << " TwoFlavourRefresh: Zero() guess "<<std::endl;
        HeatbathSolver(Vpc,PhiOdd,tmp);
	std::cout << " TwoFlavourRefresh: Heatbath solver "<<std::endl;
        Vpc.Mpc(tmp,PhiOdd);            
	std::cout << " TwoFlavourRefresh: Mpc "<<std::endl;

        // Even det factors
        DenOp.MooeeDag(etaEven,tmp);
        NumOp.MooeeInvDag(tmp,PhiEven);
	std::cout << " TwoFlavourRefresh: Mee "<<std::endl;

	RefreshAction = norm2(etaEven)+norm2(etaOdd);
	std::cout << " refresh " <<action_name()<< " action "<<RefreshAction<<std::endl;
      };

      //////////////////////////////////////////////////////
      // S = phi^dag V (Mdag M)^-1 Vdag phi
      //////////////////////////////////////////////////////
      virtual RealD Sinitial(const GaugeField &U) {
	std::cout << GridLogMessage << "Returning stored two flavour refresh action "<<RefreshAction<<std::endl;
	return RefreshAction;
      }
      virtual RealD S(const GaugeField &U) {

        NumOp.ImportGauge(U);
        DenOp.ImportGauge(U);

        SchurDifferentiableOperator<Impl> Mpc(DenOp);
        SchurDifferentiableOperator<Impl> Vpc(NumOp);

        FermionField X(NumOp.FermionRedBlackGrid());
        FermionField Y(NumOp.FermionRedBlackGrid());

        Vpc.MpcDag(PhiOdd,Y);           // Y= Vdag phi
        X=Zero();
        ActionSolver(Mpc,Y,X);          // X= (MdagM)^-1 Vdag phi
        //Mpc.Mpc(X,Y);                   // Y=  Mdag^-1 Vdag phi
        // Multiply by Ydag
        RealD action = real(innerProduct(Y,X));

        //RealD action = norm2(Y);

        // The EE factorised block; normally can replace with zero if det is constant (gauge field indept)
        // Only really clover term that creates this. Leave the EE portion as a future to do to make most
        // rapid progresss on DWF for now.
        //
        NumOp.MooeeDag(PhiEven,X);
        DenOp.MooeeInvDag(X,Y);
        action = action + norm2(Y);

        return action;
      };

      //////////////////////////////////////////////////////
      // dS/du = phi^dag dV (Mdag M)^-1 V^dag  phi
      //       - phi^dag V (Mdag M)^-1 [ Mdag dM + dMdag M ]  (Mdag M)^-1 V^dag  phi
      //       + phi^dag V (Mdag M)^-1 dV^dag  phi
      //////////////////////////////////////////////////////
      virtual void deriv(const GaugeField &U,GaugeField & dSdU) {

        NumOp.ImportGauge(U);
        DenOp.ImportGauge(U);

        SchurDifferentiableOperator<Impl> Mpc(DenOp);
        SchurDifferentiableOperator<Impl> Vpc(NumOp);

        FermionField  X(NumOp.FermionRedBlackGrid());
        FermionField  Y(NumOp.FermionRedBlackGrid());

        // This assignment is necessary to be compliant with the HMC grids
	GaugeField force(dSdU.Grid());

        //Y=Vdag phi
        //X = (Mdag M)^-1 V^dag phi
        //Y = (Mdag)^-1 V^dag  phi
        Vpc.MpcDag(PhiOdd,Y);          // Y= Vdag phi
	std::cout << GridLogMessage <<" Y "<<norm2(Y)<<std::endl;
        X=Zero();
        DerivativeSolver(Mpc,Y,X);     // X= (MdagM)^-1 Vdag phi
	std::cout << GridLogMessage <<" X "<<norm2(X)<<std::endl;
        Mpc.Mpc(X,Y);                  // Y=  Mdag^-1 Vdag phi
	std::cout << GridLogMessage <<" Y "<<norm2(Y)<<std::endl;

        // phi^dag V (Mdag M)^-1 dV^dag  phi
        Vpc.MpcDagDeriv(force , X, PhiOdd );   dSdU = force;
	std::cout << GridLogMessage <<" deriv "<<norm2(force)<<std::endl;
  
        // phi^dag dV (Mdag M)^-1 V^dag  phi
        Vpc.MpcDeriv(force , PhiOdd, X );      dSdU = dSdU+force;
	std::cout << GridLogMessage <<" deriv "<<norm2(force)<<std::endl;

        //    -    phi^dag V (Mdag M)^-1 Mdag dM   (Mdag M)^-1 V^dag  phi
        //    -    phi^dag V (Mdag M)^-1 dMdag M   (Mdag M)^-1 V^dag  phi
        Mpc.MpcDeriv(force,Y,X);              dSdU = dSdU-force;
	std::cout << GridLogMessage <<" deriv "<<norm2(force)<<std::endl;
        Mpc.MpcDagDeriv(force,X,Y);           dSdU = dSdU-force;
	std::cout << GridLogMessage <<" deriv "<<norm2(force)<<std::endl;

        // FIXME No force contribution from EvenEven assumed here
        // Needs a fix for clover.
        assert(NumOp.ConstEE() == 1);
        assert(DenOp.ConstEE() == 1);

        dSdU = -dSdU;
        
      };
    };
NAMESPACE_END(Grid);
#endif
