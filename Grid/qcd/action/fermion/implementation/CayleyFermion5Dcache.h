/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/CayleyFermion5D.cc

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Gianluca Filaci <g.filaci@ed.ac.uk>

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

#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/CayleyFermion5D.h>


NAMESPACE_BEGIN(Grid);

// Pminus fowards
// Pplus  backwards..
template<class Impl>  
void
CayleyFermion5D<Impl>::M5D(const FermionField &psi_i,
			       const FermionField &phi_i, 
			       FermionField &chi_i,
			       Vector<Coeff_t> &lower,
			       Vector<Coeff_t> &diag,
			       Vector<Coeff_t> &upper)
{
  
  chi_i.Checkerboard()=psi_i.Checkerboard();
  GridBase *grid=psi_i.Grid();
  autoView(psi , psi_i,AcceleratorRead);
  autoView(phi , phi_i,AcceleratorRead);
  autoView(chi , chi_i,AcceleratorWrite);
  assert(phi.Checkerboard() == psi.Checkerboard());

  auto pdiag = &diag[0];
  auto pupper = &upper[0];
  auto plower = &lower[0];

  int Ls =this->Ls;

  // 10 = 3 complex mult + 2 complex add
  // Flops = 10.0*(Nc*Ns) *Ls*vol (/2 for red black counting)
  M5Dcalls++;
  M5Dtime-=usecond();

  uint64_t nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{
    uint64_t ss= sss*Ls;
    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp1, tmp2;
    for(int s=0;s<Ls;s++){
      uint64_t idx_u = ss+((s+1)%Ls);
      uint64_t idx_l = ss+((s+Ls-1)%Ls);
      spProj5m(tmp1,psi(idx_u));
      spProj5p(tmp2,psi(idx_l));
      coalescedWrite(chi[ss+s],pdiag[s]*phi(ss+s)+pupper[s]*tmp1+plower[s]*tmp2);
    }
  });
  M5Dtime+=usecond();
}

template<class Impl>  
void
CayleyFermion5D<Impl>::M5Ddag(const FermionField &psi_i,
			      const FermionField &phi_i, 
			      FermionField &chi_i,
			      Vector<Coeff_t> &lower,
			      Vector<Coeff_t> &diag,
			      Vector<Coeff_t> &upper)
{
  chi_i.Checkerboard()=psi_i.Checkerboard();
  GridBase *grid=psi_i.Grid();
  autoView(psi , psi_i,AcceleratorRead);
  autoView(phi , phi_i,AcceleratorRead);
  autoView(chi , chi_i,AcceleratorWrite);
  assert(phi.Checkerboard() == psi.Checkerboard());

  auto pdiag = &diag[0];
  auto pupper = &upper[0];
  auto plower = &lower[0];

  int Ls=this->Ls;

  // Flops = 6.0*(Nc*Ns) *Ls*vol
  M5Dcalls++;
  M5Dtime-=usecond();

  uint64_t nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{
    uint64_t ss=sss*Ls;
    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp1,tmp2;
    for(int s=0;s<Ls;s++){
      uint64_t idx_u = ss+((s+1)%Ls);
      uint64_t idx_l = ss+((s+Ls-1)%Ls);
      spProj5p(tmp1,psi(idx_u));
      spProj5m(tmp2,psi(idx_l));
      coalescedWrite(chi[ss+s],pdiag[s]*phi(ss+s)+pupper[s]*tmp1+plower[s]*tmp2);
    }
  });
  M5Dtime+=usecond();
}

template<class Impl>
void
CayleyFermion5D<Impl>::MooeeInv    (const FermionField &psi_i, FermionField &chi_i)
{
  chi_i.Checkerboard()=psi_i.Checkerboard();
  GridBase *grid=psi_i.Grid();

  autoView(psi , psi_i,AcceleratorRead);
  autoView(chi , chi_i,AcceleratorWrite);

  int Ls=this->Ls;

  auto plee  = & lee [0];
  auto pdee  = & dee [0];
  auto puee  = & uee [0];
  auto pleem = & leem[0];
  auto pueem = & ueem[0];

  MooeeInvCalls++;
  MooeeInvTime-=usecond();
  uint64_t nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{
    uint64_t ss=sss*Ls;
    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp, acc, res;

    // X = Nc*Ns
    // flops = 2X + (Ls-2)(4X + 4X) + 6X + 1 + 2X + (Ls-1)(10X + 1) = -16X + Ls(1+18X) = -192 + 217*Ls flops
    // Apply (L^{\prime})^{-1} L_m^{-1}
    res = psi(ss);
    spProj5m(tmp,res);
    acc = pleem[0]*tmp;
    spProj5p(tmp,res);
    coalescedWrite(chi[ss],res);
    
    for(int s=1;s<Ls-1;s++){
      res = psi(ss+s);
      res -= plee[s-1]*tmp;
      spProj5m(tmp,res);
      acc += pleem[s]*tmp;
      spProj5p(tmp,res);
      coalescedWrite(chi[ss+s],res);
    }
    res = psi(ss+Ls-1) - plee[Ls-2]*tmp - acc;
    
    // Apply U_m^{-1} D^{-1} U^{-1}
    res = (1.0/pdee[Ls-1])*res;
    coalescedWrite(chi[ss+Ls-1],res);
    spProj5p(acc,res);
    spProj5m(tmp,res);
    for (int s=Ls-2;s>=0;s--){
      res = (1.0/pdee[s])*chi(ss+s) - puee[s]*tmp - pueem[s]*acc;
      spProj5m(tmp,res);
      coalescedWrite(chi[ss+s],res);
    }
  });

  MooeeInvTime+=usecond();
  
}

template<class Impl>
void
CayleyFermion5D<Impl>::MooeeInvDag (const FermionField &psi_i, FermionField &chi_i)
{
  chi_i.Checkerboard()=psi_i.Checkerboard();
  GridBase *grid=psi_i.Grid();
  int Ls=this->Ls;

  autoView(psi , psi_i,AcceleratorRead);
  autoView(chi , chi_i,AcceleratorWrite);

  auto plee  = & lee [0];
  auto pdee  = & dee [0];
  auto puee  = & uee [0];
  auto pleem = & leem[0];
  auto pueem = & ueem[0];

  assert(psi.Checkerboard() == psi.Checkerboard());

  MooeeInvCalls++;
  MooeeInvTime-=usecond();


  uint64_t nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{
    uint64_t ss=sss*Ls;
    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp, acc, res;

    // X = Nc*Ns
    // flops = 2X + (Ls-2)(4X + 4X) + 6X + 1 + 2X + (Ls-1)(10X + 1) = -16X + Ls(1+18X) = -192 + 217*Ls flops
    // Apply (U^{\prime})^{-dagger} U_m^{-\dagger}
    res = psi(ss);
    spProj5p(tmp,res);
    acc = conjugate(pueem[0])*tmp;
    spProj5m(tmp,res);
    coalescedWrite(chi[ss],res);
    
    for(int s=1;s<Ls-1;s++){
      res = psi(ss+s);
      res -= conjugate(puee[s-1])*tmp;
      spProj5p(tmp,res);
      acc += conjugate(pueem[s])*tmp;
      spProj5m(tmp,res);
      coalescedWrite(chi[ss+s],res);
    }
    res = psi(ss+Ls-1) - conjugate(puee[Ls-2])*tmp - acc;
    
    // Apply L_m^{-\dagger} D^{-dagger} L^{-dagger}
    res = conjugate(1.0/pdee[Ls-1])*res;
    coalescedWrite(chi[ss+Ls-1],res);
    spProj5m(acc,res);
    spProj5p(tmp,res);
    for (int s=Ls-2;s>=0;s--){
      res = conjugate(1.0/pdee[s])*chi(ss+s) - conjugate(plee[s])*tmp - conjugate(pleem[s])*acc;
      spProj5p(tmp,res);
      coalescedWrite(chi[ss+s],res);
    }
  });
  MooeeInvTime+=usecond();

}

NAMESPACE_END(Grid);
