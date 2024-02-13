/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/scalar/CovariantLaplacian.h

Copyright (C) 2016

Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
			   /*  END LEGAL */
#pragma once 

NAMESPACE_BEGIN(Grid);

struct LaplacianParams : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(LaplacianParams, 
                                  RealD, lo, 
                                  RealD, hi, 
                                  int,   MaxIter, 
                                  RealD, tolerance, 
                                  int,   degree, 
                                  int,   precision);
  
  // constructor 
  LaplacianParams(RealD lo      = 0.0, 
                  RealD hi      = 1.0, 
                  int maxit     = 1000,
                  RealD tol     = 1.0e-8, 
                  int degree    = 10,
                  int precision = 64)
    : lo(lo),
      hi(hi),
      MaxIter(maxit),
      tolerance(tol),
      degree(degree),
      precision(precision){};
};

#define LEG_LOAD(Dir)						 \
  SE = st.GetEntry(ptype, Dir, ss);				 \
  if (SE->_is_local ) {						 \
    int perm= SE->_permute;					 \
    chi = coalescedReadPermute(in[SE->_offset],ptype,perm,lane); \
  } else {							 \
    chi = coalescedRead(buf[SE->_offset],lane);			 \
  }								 \
  acceleratorSynchronise();

const std::vector<int> directions4D   ({Xdir,Ydir,Zdir,Tdir,Xdir,Ydir,Zdir,Tdir});
const std::vector<int> displacements4D({1,1,1,1,-1,-1,-1,-1});

template<class Gimpl,class Field> class CovariantAdjointLaplacianStencil : public SparseMatrixBase<Field>
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);
//  RealD kappa;

  typedef typename Field::vector_object siteObject;

  template <typename vtype> using iImplDoubledGaugeField = iVector<iScalar<iMatrix<vtype, Nc> >, Nds>;
  typedef iImplDoubledGaugeField<Simd> SiteDoubledGaugeField;
  typedef Lattice<SiteDoubledGaugeField> DoubledGaugeField;
  typedef CartesianStencil<siteObject, siteObject, DefaultImplParams> StencilImpl;

  GridBase *grid;
  StencilImpl Stencil;
  SimpleCompressor<siteObject> Compressor;
  DoubledGaugeField Uds;

  CovariantAdjointLaplacianStencil( GridBase *_grid)
    : grid(_grid),
      Stencil    (grid,8,Even,directions4D,displacements4D),
      Uds(grid){}

  CovariantAdjointLaplacianStencil(GaugeField &Umu)
    :
      grid(Umu.Grid()),
      Stencil    (grid,8,Even,directions4D,displacements4D),
      Uds(grid)
  { GaugeImport(Umu); }

  void GaugeImport (const GaugeField &Umu)
  {
    assert(grid == Umu.Grid());
    for (int mu = 0; mu < Nd; mu++) {
      auto U = PeekIndex<LorentzIndex>(Umu, mu);
      PokeIndex<LorentzIndex>(Uds, U, mu );
      U = adj(Cshift(U, mu, -1));
      PokeIndex<LorentzIndex>(Uds, U, mu + 4);
    }
  };
  
  virtual GridBase *Grid(void) { return grid; };
  virtual void  Morig(const Field &_in, Field &_out)
  {
    ///////////////////////////////////////////////
    // Halo exchange for this geometry of stencil
    ///////////////////////////////////////////////
    Stencil.HaloExchange(_in, Compressor);

    ///////////////////////////////////
    // Arithmetic expressions
    ///////////////////////////////////
//    auto st = Stencil.View(AcceleratorRead);
    autoView( st     , Stencil    , AcceleratorRead);
    auto buf = st.CommBuf();

    autoView( in     , _in    , AcceleratorRead);
    autoView( out    , _out   , AcceleratorWrite);
    autoView( U     , Uds    , AcceleratorRead);

    typedef typename Field::vector_object        vobj;
    typedef decltype(coalescedRead(in[0]))    calcObj;
    typedef decltype(coalescedRead(U[0](0))) calcLink;

    const int      Nsimd = vobj::Nsimd();
    const uint64_t NN = grid->oSites();

    accelerator_for( ss, NN, Nsimd, {

	StencilEntry *SE;
	
	const int lane=acceleratorSIMTlane(Nsimd);

	calcObj chi;
	calcObj res;
	calcObj Uchi;
	calcObj Utmp;
	calcObj Utmp2;
	calcLink UU;
	calcLink Udag;
	int ptype;

	res                 = coalescedRead(in[ss])*(-8.0);

#define LEG_LOAD_MULT(leg,polarisation)			\
	UU = coalescedRead(U[ss](polarisation));	\
	Udag = adj(UU);					\
	LEG_LOAD(leg);					\
	mult(&Utmp(), &UU, &chi());			\
	Utmp2 = adj(Utmp);				\
	mult(&Utmp(), &UU, &Utmp2());			\
	Uchi = adj(Utmp);				\
	res = res + Uchi;
	
	LEG_LOAD_MULT(0,Xp);
	LEG_LOAD_MULT(1,Yp);
	LEG_LOAD_MULT(2,Zp);
	LEG_LOAD_MULT(3,Tp);
	LEG_LOAD_MULT(4,Xm);
	LEG_LOAD_MULT(5,Ym);
	LEG_LOAD_MULT(6,Zm);
	LEG_LOAD_MULT(7,Tm);

	coalescedWrite(out[ss], res,lane);
    });

  };
  virtual void  Mnew (const Field &_in, Field &_out)
  {
    ///////////////////////////////////////////////
    // Halo exchange for this geometry of stencil
    ///////////////////////////////////////////////
//    Stencil.HaloExchange(_in, Compressor);
      std::vector<std::vector<CommsRequest_t> > requests;
      Stencil.Prepare();
  {
    GRID_TRACE("Laplace Gather");
    Stencil.HaloGather(_in,Compressor);
  }

  tracePush("Laplace Communication");
  Stencil.CommunicateBegin(requests);
  {
    GRID_TRACE("MergeSHM");
    Stencil.CommsMergeSHM(Compressor);
  }
    

    ///////////////////////////////////
    // Arithmetic expressions
    ///////////////////////////////////
//    auto st = Stencil.View(AcceleratorRead);
    autoView( st     , Stencil    , AcceleratorRead);
    auto buf = st.CommBuf();

    autoView( in     , _in    , AcceleratorRead);
    autoView( out    , _out   , AcceleratorWrite);
    autoView( U     , Uds    , AcceleratorRead);

    typedef typename Field::vector_object        vobj;
    typedef decltype(coalescedRead(in[0]))    calcObj;
    typedef decltype(coalescedRead(U[0](0))) calcLink;

    const int      Nsimd = vobj::Nsimd();
    const uint64_t NN = grid->oSites();

    accelerator_for( ss, NN, Nsimd, {

	StencilEntry *SE;
	
	const int lane=acceleratorSIMTlane(Nsimd);

	calcObj chi;
	calcObj res;
	calcObj Uchi;
	calcObj Utmp;
	calcObj Utmp2;
	calcLink UU;
	calcLink Udag;
	int ptype;

	res                 = coalescedRead(in[ss])*(-8.0);


        SE = st.GetEntry(ptype, 0, ss);				 
        if (SE->_is_local ) {
	LEG_LOAD_MULT(0,Xp);
	}
        SE = st.GetEntry(ptype, 1, ss);				 
        if (SE->_is_local ) {
	LEG_LOAD_MULT(1,Yp);
	}
        SE = st.GetEntry(ptype, 2, ss);				 
        if (SE->_is_local ) {
	LEG_LOAD_MULT(2,Zp);
	}
        SE = st.GetEntry(ptype, 3, ss);				 
        if (SE->_is_local ) {
	LEG_LOAD_MULT(3,Tp);
	}
        SE = st.GetEntry(ptype, 4, ss);				 
        if (SE->_is_local ) {
	LEG_LOAD_MULT(4,Xm);
	}
        SE = st.GetEntry(ptype, 5, ss);				 
        if (SE->_is_local ) {
	LEG_LOAD_MULT(5,Ym);
	}
        SE = st.GetEntry(ptype, 6, ss);				 
        if (SE->_is_local ) {
	LEG_LOAD_MULT(6,Zm);
	}
        SE = st.GetEntry(ptype, 7, ss);				 
        if (SE->_is_local ) {
	LEG_LOAD_MULT(7,Tm);
	}

	coalescedWrite(out[ss], res,lane);
    });

    Stencil.CommunicateComplete(requests);
  tracePop("Communication");

  {
    GRID_TRACE("Merge");
    Stencil.CommsMerge(Compressor);
  }


    accelerator_for( ss, NN, Nsimd, {

	StencilEntry *SE;
	
	const int lane=acceleratorSIMTlane(Nsimd);

	calcObj chi;
	calcObj res;
	calcObj Uchi;
	calcObj Utmp;
	calcObj Utmp2;
	calcLink UU;
	calcLink Udag;
	int ptype;

//	res                 = coalescedRead(in[ss])*(-8.0);
	res                 = coalescedRead(out[ss]);

        SE = st.GetEntry(ptype, 0, ss);				 
        if ((SE->_is_local )==0){
	LEG_LOAD_MULT(0,Xp);
	}
        SE = st.GetEntry(ptype, 1, ss);				 
        if ((SE->_is_local )==0){
	LEG_LOAD_MULT(1,Yp);
	}
        SE = st.GetEntry(ptype, 2, ss);				 
        if ((SE->_is_local )==0){
	LEG_LOAD_MULT(2,Zp);
	}
        SE = st.GetEntry(ptype, 3, ss);
        if ((SE->_is_local )==0){
	LEG_LOAD_MULT(3,Tp);
	}
        SE = st.GetEntry(ptype, 4, ss);
        if ((SE->_is_local )==0){
	LEG_LOAD_MULT(4,Xm);
	}
        SE = st.GetEntry(ptype, 5, ss);
        if ((SE->_is_local )==0){
	LEG_LOAD_MULT(5,Ym);
	}
        SE = st.GetEntry(ptype, 6, ss);
        if ((SE->_is_local )==0){
	LEG_LOAD_MULT(6,Zm);
	}
        SE = st.GetEntry(ptype, 7, ss);
        if ((SE->_is_local )==0){
	LEG_LOAD_MULT(7,Tm);
	}

	coalescedWrite(out[ss], res,lane);
    });
  };
  virtual void  M(const Field &in, Field &out) {Mnew(in,out);};
  virtual void  Mdag (const Field &in, Field &out) { M(in,out);}; // Laplacian is hermitian
  virtual  void Mdiag    (const Field &in, Field &out)                  {assert(0);}; // Unimplemented need only for multigrid
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp){assert(0);}; // Unimplemented need only for multigrid
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out)     {assert(0);}; // Unimplemented need only for multigrid
};

#undef LEG_LOAD_MULT
#undef LEG_LOAD

////////////////////////////////////////////////////////////
// Laplacian operator L on adjoint fields
//
// phi: adjoint field
// L: D_mu^dag D_mu
//
// L phi(x) = Sum_mu [ U_mu(x)phi(x+mu)U_mu(x)^dag + 
//                     U_mu(x-mu)^dag phi(x-mu)U_mu(x-mu)
//                     -2phi(x)]
//
// Operator designed to be encapsulated by
// an HermitianLinearOperator<.. , ..>
////////////////////////////////////////////////////////////

template <class Impl>
class LaplacianAdjointField: public Metric<typename Impl::Field> {
  OperatorFunction<typename Impl::Field> &Solver;
  LaplacianParams param;
  MultiShiftFunction PowerHalf;    
  MultiShiftFunction PowerInvHalf;    
//template<class Gimpl,class Field> class CovariantAdjointLaplacianStencil : public SparseMatrixBase<Field>
  CovariantAdjointLaplacianStencil<Impl,typename Impl::LinkField> LapStencil;

public:
  INHERIT_GIMPL_TYPES(Impl);

  LaplacianAdjointField(GridBase* grid, OperatorFunction<GaugeField>& S, LaplacianParams& p, const RealD k = 1.0, bool if_remez=true)
    : U(Nd, grid), Solver(S), param(p), kappa(k)
	,LapStencil(grid){
    AlgRemez remez(param.lo,param.hi,param.precision);
    std::cout<<GridLogMessage << "Generating degree "<<param.degree<<" for x^(1/2)"<<std::endl;
    if(if_remez){
    remez.generateApprox(param.degree,1,2);
    PowerHalf.Init(remez,param.tolerance,false);
    PowerInvHalf.Init(remez,param.tolerance,true);
    }
    this->triv=0;
        

  };
  LaplacianAdjointField(){this->triv=0; printf("triv=%d\n",this->Trivial());}
  void Mdir(const GaugeField&, GaugeField&, int, int){ assert(0);}
  void MdirAll(const GaugeField&, std::vector<GaugeField> &){ assert(0);}
  void Mdiag(const GaugeField&, GaugeField&){ assert(0);}

  void ImportGauge(const GaugeField& _U) {
    RealD total=0.;
    for (int mu = 0; mu < Nd; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(_U, mu);
      total += norm2(U[mu]);
    }
    LapStencil.GaugeImport (_U);

    std::cout << GridLogDebug <<"ImportGauge:norm2(U _U) = "<<total<<std::endl;
  }

  void M(const GaugeField& in, GaugeField& out) {
    // in is an antihermitian matrix
    // test
    //GaugeField herm = in + adj(in);
    //std::cout << "AHermiticity: " << norm2(herm) << std::endl;
//    std::cout << GridLogDebug <<"M:Kappa = "<<kappa<<std::endl;

    GaugeLinkField sum(in.Grid());
#if 0
    GaugeLinkField tmp(in.Grid());
    GaugeLinkField tmp2(in.Grid());

    for (int nu = 0; nu < Nd; nu++) {
      sum = Zero();
      GaugeLinkField in_nu = PeekIndex<LorentzIndex>(in, nu);
      GaugeLinkField out_nu(out.Grid());
      for (int mu = 0; mu < Nd; mu++) {
        tmp = U[mu] * Cshift(in_nu, mu, +1) * adj(U[mu]);
        tmp2 = adj(U[mu]) * in_nu * U[mu];
        sum += tmp + Cshift(tmp2, mu, -1) - 2.0 * in_nu;
      }
      out_nu = (1.0 - kappa) * in_nu - kappa / (double(4 * Nd)) * sum;
      PokeIndex<LorentzIndex>(out, out_nu, nu);
    }
#else
    for (int nu = 0; nu < Nd; nu++) {
      GaugeLinkField in_nu = PeekIndex<LorentzIndex>(in, nu);
      GaugeLinkField out_nu(out.Grid());
      LapStencil.M(in_nu,sum);
      out_nu = (1.0 - kappa) * in_nu - kappa / (double(4 * Nd)) * sum;
      PokeIndex<LorentzIndex>(out, out_nu, nu);
    }
#endif
//    std::cout << GridLogDebug <<"M:norm2(out) = "<<norm2(out)<<std::endl;
  }

#if 0
  void Quad(const GaugeField& in, GaugeField& out,RealD a0,RealD a1,RealD a2) {

    GaugeLinkField tmp(in.Grid());
    GaugeLinkField tmp2(in.Grid());
#if 0
    std::vector<GaugeLinkField> sum(in.Grid(),Nd);
    std::vector<GaugeLinkField> sum2(in.Grid(),Nd);
    std::vector<GaugeLinkField> in_nu(in.Grid(),Nd);
    std::vector<GaugeLinkField> out_nu(in.Grid(),Nd);

    for (int nu = 0; nu < Nd; nu++) {
      sum[nu] = Zero();
      in_nu[nu] = PeekIndex<LorentzIndex>(in, nu);
      out_nu[nu] = a0*in_nu[nu];
      for (int mu = 0; mu < Nd; mu++) {
        tmp = U[mu] * Cshift(in_nu[nu], mu, +1) * adj(U[mu]);
        tmp2 = adj(U[mu]) * in_nu[nu] * U[mu];
        sum[nu] += tmp + Cshift(tmp2, mu, -1) - 2.0 * in_nu;
      }
      out_nu[nu] +=  a1*  1. / (double(4 * Nd)) * sum[nu];
      sum2[nu] = Zero();
      for (int mu = 0; mu < Nd; mu++) {
        tmp = U[mu] * Cshift(sum[nu], mu, +1) * adj(U[mu]);
        tmp2 = adj(U[mu]) * in_nu * U[mu];
        sum2[nu] += tmp + Cshift(tmp2, mu, -1) - 2.0 * in_nu;
      }
      out_nu[nu] +=  a2* ( 1. / (double(4 * Nd)))^2 * sum[nu];
      PokeIndex<LorentzIndex>(out, out_nu[nu], nu);
    }
#else
    for (int nu = 0; nu < Nd; nu++) {
      GaugeLinkField in_nu = PeekIndex<LorentzIndex>(in, nu);
      GaugeLinkField out_nu(out.Grid());
      GaugeLinkField sum(out.Grid());
      GaugeLinkField sum2(out.Grid());
      out_nu=a0*in_nu;
      LapStencil.M(in_nu,sum);
      out_nu +=  a1*  1. / (double(4 * Nd)) * sum;
      LapStencil.M(sum,sum2);
      out_nu +=  a2* ( 1. / (double(4 * Nd)))^2 * sum2;
//      out_nu += (1.0 - kappa) * in_nu - kappa / (double(4 * Nd)) * sum;
      PokeIndex<LorentzIndex>(out, out_nu, nu);
    }
#endif
  }
#endif

  void MDeriv(const GaugeField& in, GaugeField& der) {
    // in is anti-hermitian
//    std::cout << GridLogDebug <<"MDeriv:Kappa = "<<kappa<<std::endl;
    RealD factor = -kappa / (double(4 * Nd));
    
    for (int mu = 0; mu < Nd; mu++){
      GaugeLinkField der_mu(der.Grid());
      der_mu = Zero();
      for (int nu = 0; nu < Nd; nu++){
        GaugeLinkField in_nu = PeekIndex<LorentzIndex>(in, nu);
        der_mu += U[mu] * Cshift(in_nu, mu, 1) * adj(U[mu]) * in_nu;
      }
      // the minus sign comes by using the in_nu instead of the
      // adjoint in the last multiplication
      PokeIndex<LorentzIndex>(der,  -2.0 * factor * der_mu, mu);
    } 
    std::cout << GridLogDebug <<"MDeriv: Kappa= "<< kappa << " norm2(der) = "<<norm2(der)<<std::endl;
  }

  // separating this temporarily
  void MDeriv(const GaugeField& left, const GaugeField& right,
              GaugeField& der) {
    // in is anti-hermitian
    RealD factor = -kappa / (double(4 * Nd));

    for (int mu = 0; mu < Nd; mu++) {
      GaugeLinkField der_mu(der.Grid());
      der_mu = Zero();
      for (int nu = 0; nu < Nd; nu++) {
        GaugeLinkField left_nu = PeekIndex<LorentzIndex>(left, nu);
        GaugeLinkField right_nu = PeekIndex<LorentzIndex>(right, nu);
        der_mu += U[mu] * Cshift(left_nu, mu, 1) * adj(U[mu]) * right_nu;
        der_mu += U[mu] * Cshift(right_nu, mu, 1) * adj(U[mu]) * left_nu;
      }
      PokeIndex<LorentzIndex>(der, -factor * der_mu, mu);
    }
    std::cout << GridLogDebug <<"MDeriv: Kappa= "<< kappa << " norm2(der) = "<<norm2(der)<<std::endl;
  }

  void Minv(const GaugeField& in, GaugeField& inverted){
    HermitianLinearOperator<LaplacianAdjointField<Impl>,GaugeField> HermOp(*this);
    Solver(HermOp, in, inverted);
    std::cout << GridLogDebug <<"Minv:norm2(inverted) = "<<norm2(inverted)<<std::endl;
  }


  void MinvDeriv(const GaugeField& in, GaugeField& der) {
    GaugeField X(in.Grid());
    Minv(in,X);
    MDeriv(X,der);
    der *=-1.0;
    std::cout << GridLogDebug <<"MinvDeriv:norm2(der) = "<<norm2(der)<<std::endl;
  }

  void MSquareRoot(GaugeField& P){
    GaugeField Gp(P.Grid());
    HermitianLinearOperator<LaplacianAdjointField<Impl>,GaugeField> HermOp(*this);
    ConjugateGradientMultiShift<GaugeField> msCG(param.MaxIter,PowerHalf);
    msCG(HermOp,P,Gp);
    P = Gp; 
    std::cout << GridLogDebug <<"MSquareRoot:norm2(P) = "<<norm2(P)<<std::endl;
  }

  void MInvSquareRoot(GaugeField& P){
    GaugeField Gp(P.Grid());
    HermitianLinearOperator<LaplacianAdjointField<Impl>,GaugeField> HermOp(*this);
    ConjugateGradientMultiShift<GaugeField> msCG(param.MaxIter,PowerInvHalf);
    msCG(HermOp,P,Gp);
    P = Gp; 
    std::cout << GridLogDebug <<"MInvSquareRoot:norm2(P) = "<<norm2(P)<<std::endl;
  }



private:
  RealD kappa;
  std::vector<GaugeLinkField> U;
};

NAMESPACE_END(Grid);
