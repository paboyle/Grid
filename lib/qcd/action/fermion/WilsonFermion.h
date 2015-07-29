#ifndef  GRID_QCD_WILSON_FERMION_H
#define  GRID_QCD_WILSON_FERMION_H

namespace Grid {

  namespace QCD {

    class WilsonFermion : public FermionOperator<LatticeFermion,LatticeGaugeField>
    {
    public:

      ///////////////////////////////////////////////////////////////
      // Implement the abstract base
      ///////////////////////////////////////////////////////////////
      GridBase *GaugeGrid(void)              { return _grid ;}
      GridBase *GaugeRedBlackGrid(void)      { return _cbgrid ;}
      GridBase *FermionGrid(void)            { return _grid;}
      GridBase *FermionRedBlackGrid(void)    { return _cbgrid;}

      // override multiply
      virtual RealD  M    (const LatticeFermion &in, LatticeFermion &out);
      virtual RealD  Mdag (const LatticeFermion &in, LatticeFermion &out);

      // half checkerboard operaions
      void   Meooe       (const LatticeFermion &in, LatticeFermion &out);
      void   MeooeDag    (const LatticeFermion &in, LatticeFermion &out);

      virtual void   Mooee       (const LatticeFermion &in, LatticeFermion &out); // remain virtual so we 
      virtual void   MooeeDag    (const LatticeFermion &in, LatticeFermion &out); // can derive Clover
      virtual void   MooeeInv    (const LatticeFermion &in, LatticeFermion &out); // from Wilson base
      virtual void   MooeeInvDag (const LatticeFermion &in, LatticeFermion &out);

      ////////////////////////
      // 
      // Force term: d/dtau S = 0
      //
      // It is simplest to consider the two flavour force term 
      // 
      //       S[U,phi] = phidag (MdagM)^-1 phi
      //
      // But simplify even this to
      //
      //       S[U,phi] = phidag MdagM phi
      //
      // (other options exist depending on nature of action fragment.)
      // 
      // Require momentum be traceless anti-hermitian to move within group manifold [ P = i P^a T^a ]
      //
      // Define the HMC hamiltonian
      //
      //       H = 1/2 Tr P^2 + S(U,phi)
      // 
      // .
      // U =  P U    (lorentz & color indices multiplied)
      //
      // Hence
      //
      // .c    c  c       c
      // U =  U  P   = - U  P        (c == dagger)
      //
      // So, taking some liberty with implicit indices
      //                  .      .          .c      c
      // dH/dt = 0 = Tr P P +Tr[ U  dS/dU + U  dS/dU  ]
      //
      //                .                          c      c
      //           = Tr P P + i Tr[  P U dS/dU  - U   P dS/dU  ]
      //
      //                   .                        c  c
      //           = Tr P (P + i ( U dS/dU - P dS/dU  U ]
      //
      //              .                        c  c
      //           => P  = -i [ U dS/dU - dS/dU  U ]      generates HMC EoM
      //
      // Simple case work this out using S = phi^dag MdagM phi for wilson:
      //                               c       c
      // dSdt     = dU_xdt  dSdUx  + dUxdt dSdUx
      //          
      //         = Tr i P U_x [ (\phi^\dag)_x (1+g) (M \phi)_x+\mu   +(\phi^\dag M^\dag)_x (1-g) \phi_{x+\mu} ]
      //                 c
      //            - i U_x P [ (\phi^\dag)_x+mu (1-g) (M \phi)_x    +(\phi^\dag M^\dag)_(x+\mu) (1+g) \phi_{x} ]
      //
      //         = i [(\phi^\dag)_x      ]_j P_jk  [U_x(1+g) (M \phi)_x+\mu]_k        (1)
      //         + i [(\phi^\dagM^\dag)_x]_j P_jk  [U_x(1-g) (\phi)_x+\mu]_k          (2)
      //         - i [(\phi^\dag)_x+mu (1-g) U^dag_x]_j P_jk  [(M \phi)_xk            (3)
      //         - i [(\phi^\dagM^\dag)_x+mu (1+g) U^dag_x]_j P_jk  [ \phi]_xk        (4)
      //
      // Observe that (1)* = (4)
      //              (2)* = (3)
      //
      // Write as    .
      //             P_{kj}  = - i (  [U_x(1+g) (M \phi)_x+\mu] (x) [(\phi^\dag)_x] + [U_x(1-g) (\phi)_x+\mu] (x) [(\phi^\dagM^\dag)_x] - h.c )
      //
      // where (x) denotes outer product in colour and spins are traced.
      //
      // Need only evaluate (1) and (2)  [Chroma] or (2) and (4) [IroIro] and take the 
      // traceless anti hermitian part (of term in brackets w/o the "i")
      // 
      // Generalisation to S=phi^dag (MdagM)^{-1} phi is simple:
      //
      // For more complicated DWF etc... apply product rule in differentiation
      //
      ////////////////////////
      void DhopDeriv  (LatticeGaugeField &mat,const LatticeFermion &U,const LatticeFermion &V,int dag);
      void DhopDerivEO(LatticeGaugeField &mat,const LatticeFermion &U,const LatticeFermion &V,int dag);
      void DhopDerivOE(LatticeGaugeField &mat,const LatticeFermion &U,const LatticeFermion &V,int dag);

      // Extra support internal
      void DerivInternal(CartesianStencil & st,
			 LatticeDoubledGaugeField & U,
			 LatticeGaugeField &mat,
			 const LatticeFermion &A,
			 const LatticeFermion &B,
			 int dag);


      // non-hermitian hopping term; half cb or both
      void Dhop  (const LatticeFermion &in, LatticeFermion &out,int dag);
      void DhopOE(const LatticeFermion &in, LatticeFermion &out,int dag);
      void DhopEO(const LatticeFermion &in, LatticeFermion &out,int dag);

      // Multigrid assistance
      void   Mdir (const LatticeFermion &in, LatticeFermion &out,int dir,int disp);
      void DhopDir(const LatticeFermion &in, LatticeFermion &out,int dir,int disp);
      void DhopDirDisp(const LatticeFermion &in, LatticeFermion &out,int dirdisp,int gamma,int dag);

      ///////////////////////////////////////////////////////////////
      // Extra methods added by derived
      ///////////////////////////////////////////////////////////////
      void DhopInternal(CartesianStencil & st,
			LatticeDoubledGaugeField &U,
			const LatticeFermion &in, 
			LatticeFermion &out,
			int dag);

      // Constructor
      WilsonFermion(LatticeGaugeField &_Umu,GridCartesian &Fgrid,GridRedBlackCartesian &Hgrid,RealD _mass);

      // DoubleStore
      virtual void ImportGauge(const LatticeGaugeField &_Umu);
      void DoubleStore(LatticeDoubledGaugeField &Uds,const LatticeGaugeField &Umu);

      ///////////////////////////////////////////////////////////////
      // Data members require to support the functionality
      ///////////////////////////////////////////////////////////////
      static int HandOptDslash; // these are a temporary hack
      static int MortonOrder;

      //    protected:
    public:

      RealD                        mass;

      GridBase                     *    _grid; 
      GridBase                     *  _cbgrid;

      static const int npoint=8;
      static const std::vector<int> directions   ;
      static const std::vector<int> displacements;

      //Defines the stencils for even and odd
      CartesianStencil Stencil; 
      CartesianStencil StencilEven; 
      CartesianStencil StencilOdd; 

      // Copy of the gauge field , with even and odd subsets
      LatticeDoubledGaugeField Umu;
      LatticeDoubledGaugeField UmuEven;
      LatticeDoubledGaugeField UmuOdd;

      // Comms buffer
      std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  comm_buf;

      
    };

  }
}
#endif
