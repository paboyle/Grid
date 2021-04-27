/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: ./lib/qcd/smearing/StoutSmearing.h
 
 Copyright (C) 2019
 
 Author: unknown
 Author: Felix Erben <ferben@ed.ac.uk>
 Author: Michael Marshall <Michael.Marshall@ed.ac.uk>

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
/*
  @file StoutSmearing.h
  @brief Declares Stout smearing class
*/
#pragma once

NAMESPACE_BEGIN(Grid);

/*!  @brief Stout smearing of link variable. */
template <class Gimpl>
class Smear_Stout : public Smear<Gimpl> {
 private:
  int OrthogDim = -1;
  const std::vector<double> SmearRho;
  // Smear<Gimpl>* ownership semantics:
  //    Smear<Gimpl>* passed in to constructor are owned by caller, so we don't delete them here
  //    Smear<Gimpl>* created within constructor need to be deleted as part of the destructor
  const std::unique_ptr<Smear<Gimpl>> OwnedBase; // deleted at destruction
  const Smear<Gimpl>* SmearBase; // Not owned by this object, so not deleted at destruction

  // only anticipated to be used from default constructor
  inline static std::vector<double> rho3D(double rho, int orthogdim){
    std::vector<double> rho3d(Nd*Nd);
    for (int mu=0; mu<Nd; mu++)
      for (int nu=0; nu<Nd; nu++)
        rho3d[mu + Nd * nu] = (mu == nu || mu == orthogdim || nu == orthogdim) ? 0.0 : rho;
    return rho3d;
  };
  
public:
  INHERIT_GIMPL_TYPES(Gimpl)

  /*! Stout smearing with base explicitly specified */
  Smear_Stout(Smear<Gimpl>* base) : SmearBase{base} {
    assert(Nc == 3 && "Stout smearing currently implemented only for Nc==3");
  }

  /*! Construct stout smearing object from explicitly specified rho matrix */
  Smear_Stout(const std::vector<double>& rho_)
    : OwnedBase{new Smear_APE<Gimpl>(rho_)}, SmearBase{OwnedBase.get()} {
    std::cout << GridLogDebug << "Stout smearing constructor : Smear_Stout(const std::vector<double>& " << rho_ << " )" << std::endl
    assert(Nc == 3 && "Stout smearing currently implemented only for Nc==3");
    }

  /*! Default constructor. rho is constant in all directions, optionally except for orthogonal dimension */
  Smear_Stout(double rho = 1.0, int orthogdim = -1)
  : OrthogDim{orthogdim}, SmearRho{ rho3D(rho,orthogdim) }, OwnedBase{ new Smear_APE<Gimpl>(SmearRho) }, SmearBase{OwnedBase.get()} {
    assert(Nc == 3 && "Stout smearing currently implemented only for Nc==3");
  }

  ~Smear_Stout() {}  // delete SmearBase...

  void smear(GaugeField& u_smr, const GaugeField& U) const {
    GaugeField C(U.Grid());
    GaugeLinkField tmp(U.Grid()), iq_mu(U.Grid()), Umu(U.Grid());

    std::cout << GridLogDebug << "Stout smearing started\n";

    // C contains the staples multiplied by some rho
    u_smr = U ; // set the smeared field to the current gauge field
    SmearBase->smear(C, U);

    for (int mu = 0; mu < Nd; mu++) {
      if( mu == OrthogDim ) continue ;
      // u_smr = exp(iQ_mu)*U_mu apart from Orthogdim
      Umu = peekLorentz(U, mu);
      tmp = peekLorentz(C, mu);
      iq_mu = Ta( tmp * adj(Umu));  
      exponentiate_iQ(tmp, iq_mu);
      pokeLorentz(u_smr, tmp * Umu, mu);
    }
    std::cout << GridLogDebug << "Stout smearing completed\n";
  };

  void derivative(GaugeField& SigmaTerm, const GaugeField& iLambda,
                  const GaugeField& Gauge) const {
    SmearBase->derivative(SigmaTerm, iLambda, Gauge);
  };

  void BaseSmear(GaugeField& C, const GaugeField& U) const {
    SmearBase->smear(C, U);
  };


  // Repetion of code here (use the Tensor_exp.h function)
  void exponentiate_iQ(GaugeLinkField& e_iQ, const GaugeLinkField& iQ) const {
    // Put this outside
    // only valid for SU(3) matrices

    // only one Lorentz direction at a time

    // notice that it actually computes
    // exp ( input matrix )
    // the i sign is coming from outside
    // input matrix is anti-hermitian NOT hermitian

    GridBase* grid = iQ.Grid();
    GaugeLinkField unity(grid);
    unity = 1.0;

    GaugeLinkField iQ2(grid), iQ3(grid);
    LatticeComplex u(grid), w(grid);
    LatticeComplex f0(grid), f1(grid), f2(grid);

    iQ2 = iQ * iQ;
    iQ3 = iQ * iQ2;

    //We should check sgn(c0) here already and then apply eq (34) from 0311018
    set_uw(u, w, iQ2, iQ3);
    set_fj(f0, f1, f2, u, w);

    e_iQ = f0 * unity + timesMinusI(f1) * iQ - f2 * iQ2;
  };

  void set_uw(LatticeComplex& u, LatticeComplex& w, GaugeLinkField& iQ2,
              GaugeLinkField& iQ3) const {
    Complex one_over_three = 1.0 / 3.0;
    Complex one_over_two = 1.0 / 2.0;

    GridBase* grid = u.Grid();
    LatticeComplex c0(grid), c1(grid), tmp(grid), c0max(grid), theta(grid);

    // sign in c0 from the conventions on the Ta
    c0 = -imag(trace(iQ3)) * one_over_three;  
    c1 = -real(trace(iQ2)) * one_over_two;

    // Cayley Hamilton checks to machine precision, tested
    tmp = c1 * one_over_three;
    c0max = 2.0 * pow(tmp, 1.5);

    theta = acos(c0 / c0max) *
      one_over_three;  // divide by three here, now leave as it is
    u = sqrt(tmp) * cos(theta);
    w = sqrt(c1) * sin(theta);
  }

  void set_fj(LatticeComplex& f0, LatticeComplex& f1, LatticeComplex& f2,
              const LatticeComplex& u, const LatticeComplex& w) const {
    GridBase* grid = u.Grid();
    LatticeComplex xi0(grid), u2(grid), w2(grid), cosw(grid);
    LatticeComplex fden(grid);
    LatticeComplex h0(grid), h1(grid), h2(grid);
    LatticeComplex e2iu(grid), emiu(grid), ixi0(grid), qt(grid);
    LatticeComplex unity(grid);
    unity = 1.0;

    xi0 = func_xi0(w);
    u2 = u * u;
    w2 = w * w;
    cosw = cos(w);

    ixi0 = timesI(xi0);
    emiu = cos(u) - timesI(sin(u));
    e2iu = cos(2.0 * u) + timesI(sin(2.0 * u));

    h0 = e2iu * (u2 - w2) +
      emiu * ((8.0 * u2 * cosw) + (2.0 * u * (3.0 * u2 + w2) * ixi0));
    h1 = e2iu * (2.0 * u) - emiu * ((2.0 * u * cosw) - (3.0 * u2 - w2) * ixi0);
    h2 = e2iu - emiu * (cosw + (3.0 * u) * ixi0);

    fden = unity / (9.0 * u2 - w2);  // reals
    f0 = h0 * fden;
    f1 = h1 * fden;
    f2 = h2 * fden;
  }

  LatticeComplex func_xi0(const LatticeComplex& w) const {
    // Definition from arxiv 0311018
    //if (abs(w) < 0.05) {w2 = w*w; return 1.0 - w2/6.0 * (1.0-w2/20.0 * (1.0-w2/42.0));}
    return sin(w) / w;
  }

  LatticeComplex func_xi1(const LatticeComplex& w) const {
    // Define a function to do the check
    // if( w < 1e-4 ) std::cout << GridLogWarning << "[Smear_stout] w too small:
    // "<< w <<"\n";
    return cos(w) / (w * w) - sin(w) / (w * w * w);
  }
};

NAMESPACE_END(Grid);
