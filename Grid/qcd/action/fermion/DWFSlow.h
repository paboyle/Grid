/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/DWFSlow.h

Copyright (C) 2022

Author: Peter Boyle <pboyle@bnl.gov>

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

template <class Impl>
class DWFSlowFermion : public FermionOperator<Impl>
{
public:
  INHERIT_IMPL_TYPES(Impl);

  ///////////////////////////////////////////////////////////////
  // Implement the abstract base
  ///////////////////////////////////////////////////////////////
  GridBase *GaugeGrid(void) { return _grid4; }
  GridBase *GaugeRedBlackGrid(void) { return _cbgrid4; }
  GridBase *FermionGrid(void) { return _grid; }
  GridBase *FermionRedBlackGrid(void) { return _cbgrid; }

  FermionField _tmp;
  FermionField &tmp(void) { return _tmp; }

  //////////////////////////////////////////////////////////////////
  // override multiply; cut number routines if pass dagger argument
  // and also make interface more uniformly consistent
  //////////////////////////////////////////////////////////////////
  virtual void  M(const FermionField &in, FermionField &out)
  {
    FermionField tmp(_grid);
    out = (5.0 - M5) * in;
    Dhop(in,tmp,DaggerNo);
    out = out + tmp;
  }
  virtual void  Mdag(const FermionField &in, FermionField &out)
  {
    FermionField tmp(_grid);
    out = (5.0 - M5) * in;
    Dhop(in,tmp,DaggerYes);
    out = out + tmp;
  };

  /////////////////////////////////////////////////////////
  // half checkerboard operations 5D redblack so just site identiy
  /////////////////////////////////////////////////////////
  void Meooe(const FermionField &in, FermionField &out)
  {
    if ( in.Checkerboard() == Odd ) {
      this->DhopEO(in,out,DaggerNo);
    } else {
      this->DhopOE(in,out,DaggerNo);
    }
  }
  void MeooeDag(const FermionField &in, FermionField &out)
  {
    if ( in.Checkerboard() == Odd ) {
      this->DhopEO(in,out,DaggerYes);
    } else {
      this->DhopOE(in,out,DaggerYes);
    }
  };

  // allow override for twisted mass and clover
  virtual void Mooee(const FermionField &in, FermionField &out)
  {
    out = (5.0 - M5) * in;
  }
  virtual void MooeeDag(const FermionField &in, FermionField &out)
  {
    out = (5.0 - M5) * in;
  }
  virtual void MooeeInv(const FermionField &in, FermionField &out)
  {
    out = (1.0/(5.0 - M5)) * in;
  };
  virtual void MooeeInvDag(const FermionField &in, FermionField &out)
  {
    out = (1.0/(5.0 - M5)) * in;
  };

  virtual void  MomentumSpacePropagator(FermionField &out,const FermionField &in,RealD _mass,std::vector<double> twist) {} ;

  ////////////////////////
  // Derivative interface
  ////////////////////////
  // Interface calls an internal routine
  void DhopDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag)  { assert(0);};
  void DhopDerivOE(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){ assert(0);};
  void DhopDerivEO(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){ assert(0);};

  ///////////////////////////////////////////////////////////////
  // non-hermitian hopping term; half cb or both
  ///////////////////////////////////////////////////////////////
  void Dhop(const FermionField &in, FermionField &out, int dag)
  {
    FermionField tmp(in.Grid());
    Dhop5(in,out,MassField,MassField,dag );
    for(int mu=0;mu<4;mu++){
      DhopDirU(in,Umu[mu],Umu[mu],tmp,mu,dag );    out = out + tmp;
    }
  };
  void DhopOE(const FermionField &in, FermionField &out, int dag)
  {
    FermionField tmp(in.Grid());
    assert(in.Checkerboard()==Even);
    Dhop5(in,out,MassFieldOdd,MassFieldEven,dag);
    for(int mu=0;mu<4;mu++){
      DhopDirU(in,UmuOdd[mu],UmuEven[mu],tmp,mu,dag );    out = out + tmp;
    }
  };
  void DhopEO(const FermionField &in, FermionField &out, int dag)
  {
    FermionField tmp(in.Grid());
    assert(in.Checkerboard()==Odd);
    Dhop5(in,out, MassFieldEven,MassFieldOdd ,dag );  
    for(int mu=0;mu<4;mu++){
      DhopDirU(in,UmuEven[mu],UmuOdd[mu],tmp,mu,dag );    out = out + tmp;
    }
  };

  ///////////////////////////////////////////////////////////////
  // Multigrid assistance; force term uses too
  ///////////////////////////////////////////////////////////////
  void Mdir(const FermionField &in, FermionField &out, int dir, int disp){ assert(0);};
  void MdirAll(const FermionField &in, std::vector<FermionField> &out)   { assert(0);};
  void DhopDir(const FermionField &in, FermionField &out, int dir, int disp) { assert(0);};
  void DhopDirAll(const FermionField &in, std::vector<FermionField> &out)    { assert(0);};
  void DhopDirCalc(const FermionField &in, FermionField &out, int dirdisp,int gamma, int dag) { assert(0);};

  void DhopDirU(const FermionField &in, const GaugeLinkField &U5e, const GaugeLinkField &U5o, FermionField &out, int mu, int dag)
  {
    RealD     sgn= 1.0;
    if (dag ) sgn=-1.0;

    Gamma::Algebra Gmu [] = {
			 Gamma::Algebra::GammaX,
			 Gamma::Algebra::GammaY,
			 Gamma::Algebra::GammaZ,
			 Gamma::Algebra::GammaT
    };

    //    mass is  1,1,1,1,-m has to multiply the round the world term
    FermionField tmp (in.Grid());
    tmp = U5e * Cshift(in,mu+1,1);
    out = tmp - Gamma(Gmu[mu])*tmp*sgn;
    
    tmp = Cshift(adj(U5o)*in,mu+1,-1);
    out = out + tmp + Gamma(Gmu[mu])*tmp*sgn;

    out = -0.5*out;
  };

  void Dhop5(const FermionField &in, FermionField &out, ComplexField &massE, ComplexField &massO, int dag)
  {
    // Mass term.... must multiple the round world with mass = 1,1,1,1, -m
    RealD     sgn= 1.0;
    if (dag ) sgn=-1.0;

    Gamma G5(Gamma::Algebra::Gamma5);

    FermionField tmp (in.Grid());
    tmp = massE*Cshift(in,0,1);
    out = tmp - G5*tmp*sgn;
    
    tmp = Cshift(massO*in,0,-1);
    out = out + tmp + G5*tmp*sgn;
    out = -0.5*out;
  };

  // Constructor
  DWFSlowFermion(GaugeField &_Umu, GridCartesian &Fgrid,
		 GridRedBlackCartesian &Hgrid, RealD _mass, RealD _M5)
    :
    _grid(&Fgrid),
    _cbgrid(&Hgrid),
    _grid4(_Umu.Grid()),
    Umu(Nd,&Fgrid),
    UmuEven(Nd,&Hgrid),
    UmuOdd(Nd,&Hgrid),
    MassField(&Fgrid),
    MassFieldEven(&Hgrid),
    MassFieldOdd(&Hgrid),
    M5(_M5),
    mass(_mass),
    _tmp(&Hgrid)
    {
      Ls=Fgrid._fdimensions[0];
      ImportGauge(_Umu);

      typedef typename FermionField::scalar_type scalar;

      Lattice<iScalar<vInteger> > coor(&Fgrid);
      LatticeCoordinate(coor, 0); // Scoor
      ComplexField one(&Fgrid);
      MassField =scalar(-mass);
      one       =scalar(1.0);
      MassField =where(coor==Integer(Ls-1),MassField,one);
      for(int mu=0;mu<Nd;mu++){
	pickCheckerboard(Even,UmuEven[mu],Umu[mu]);
	pickCheckerboard(Odd ,UmuOdd[mu],Umu[mu]);
      }
      pickCheckerboard(Even,MassFieldEven,MassField);
      pickCheckerboard(Odd ,MassFieldOdd,MassField);
      
    }
  
  // DoubleStore impl dependent
  void ImportGauge(const GaugeField &_Umu4)
  {
    GaugeLinkField U4(_grid4);
    for(int mu=0;mu<Nd;mu++){
      U4 = PeekIndex<LorentzIndex>(_Umu4, mu);
      for(int s=0;s<this->Ls;s++){
	InsertSlice(U4,Umu[mu],s,0);
      }
    }
  }

  ///////////////////////////////////////////////////////////////
  // Data members require to support the functionality
  ///////////////////////////////////////////////////////////////

public:
  virtual RealD Mass(void) { return mass; }
  virtual int   isTrivialEE(void) { return 1; };
  RealD mass;
  RealD M5;
  int Ls;

  GridBase *_grid4;
  GridBase *_grid;
  GridBase *_cbgrid4;
  GridBase *_cbgrid;

  // Copy of the gauge field , with even and odd subsets
  std::vector<GaugeLinkField> Umu;
  std::vector<GaugeLinkField> UmuEven;
  std::vector<GaugeLinkField> UmuOdd;
  ComplexField MassField;
  ComplexField MassFieldEven;
  ComplexField MassFieldOdd;

  ///////////////////////////////////////////////////////////////
  // Conserved current utilities
  ///////////////////////////////////////////////////////////////
  void ContractConservedCurrent(PropagatorField &q_in_1,
                                PropagatorField &q_in_2,
                                PropagatorField &q_out,
                                PropagatorField &phys_src,
                                Current curr_type,
                                unsigned int mu){}
  void SeqConservedCurrent(PropagatorField &q_in,
                           PropagatorField &q_out,
                           PropagatorField &phys_src,
                           Current curr_type,
                           unsigned int mu,
                           unsigned int tmin,
			   unsigned int tmax,
			   ComplexField &lattice_cmplx){}
};

typedef DWFSlowFermion<WilsonImplF> DWFSlowFermionF;
typedef DWFSlowFermion<WilsonImplD> DWFSlowFermionD;

NAMESPACE_END(Grid);
