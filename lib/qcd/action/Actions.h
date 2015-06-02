#ifndef GRID_QCD_ACTIONS_H
#define GRID_QCD_ACTIONS_H


// Some reorganisation likely required as both Chroma and IroIro
// are separating the concept of the operator from that of action.
//
// The FermAction contains methods to create 
//
// * Linear operators             (Hermitian and non-hermitian)  .. my LinearOperator
// * System solvers               (Hermitian and non-hermitian)  .. my OperatorFunction
// * MultiShift System solvers    (Hermitian and non-hermitian)  .. my OperatorFunction


////////////////////////////////////////////
// Abstract base interface
////////////////////////////////////////////
#include <qcd/action/fermion/FermionOperator.h>

////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////
#include <qcd/action/fermion/WilsonCompressor.h>     //used by all wilson type fermions
#include <qcd/action/fermion/WilsonKernels.h>        //used by all wilson type fermions

////////////////////////////////////////////
// 4D formulations
////////////////////////////////////////////
#include <qcd/action/fermion/WilsonFermion.h>
//#include <qcd/action/fermion/CloverFermion.h>

////////////////////////////////////////////
// 5D formulations
////////////////////////////////////////////
#include <qcd/action/fermion/WilsonFermion5D.h> // used by all 5d overlap types
#include <qcd/action/fermion/CayleyFermion5D.h>
#include <qcd/action/fermion/ContinuedFractionFermion5D.h>
//#include <qcd/action/fermion/PartialFraction.h>

#include <qcd/action/fermion/DomainWallFermion.h>
//#include <qcd/action/fermion/ScaledShamirCayleyTanh.h>


    // Chroma interface defining FermionAction
    /*
     template<typename T, typename P, typename Q>  class FermAct4D : public FermionAction<T,P,Q>
     virtual LinearOperator<T>* linOp(Handle< FermState<T,P,Q> > state) const = 0;
     virtual LinearOperator<T>* lMdagM(Handle< FermState<T,P,Q> > state) const = 0;
     virtual LinOpSystemSolver<T>* invLinOp(Handle< FermState<T,P,Q> > state,
     virtual MdagMSystemSolver<T>* invMdagM(Handle< FermState<T,P,Q> > state,
     virtual LinOpMultiSystemSolver<T>* mInvLinOp(Handle< FermState<T,P,Q> > state,
     virtual MdagMMultiSystemSolver<T>* mInvMdagM(Handle< FermState<T,P,Q> > state,
     virtual MdagMMultiSystemSolverAccumulate<T>* mInvMdagMAcc(Handle< FermState<T,P,Q> > state,
     virtual SystemSolver<T>* qprop(Handle< FermState<T,P,Q> > state,
     class DiffFermAct4D : public FermAct4D<T,P,Q>
     virtual DiffLinearOperator<T,Q,P>* linOp(Handle< FermState<T,P,Q> > state) const = 0;
     virtual DiffLinearOperator<T,Q,P>* lMdagM(Handle< FermState<T,P,Q> > state) const = 0;
    */


    // Chroma interface defining GaugeAction
    /*
      template<typename P, typename Q>   class GaugeAction
  virtual const CreateGaugeState<P,Q>& getCreateState() const = 0;
  virtual GaugeState<P,Q>* createState(const Q& q) const
  virtual const GaugeBC<P,Q>& getGaugeBC() const
  virtual const Set& getSet(void) const = 0;
  virtual void deriv(P& result, const Handle< GaugeState<P,Q> >& state) const 
  virtual Double S(const Handle< GaugeState<P,Q> >& state) const = 0;

  class LinearGaugeAction : public GaugeAction< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  typedef multi1d<LatticeColorMatrix>  P;
  typedef multi1d<LatticeColorMatrix>  Q;
  virtual void staple(LatticeColorMatrix& result,
		      const Handle< GaugeState<P,Q> >& state,
		      int mu, int cb) const = 0;
    */


#endif
