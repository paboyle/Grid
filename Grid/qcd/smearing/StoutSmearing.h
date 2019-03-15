/*
  @file stoutSmear.hpp
  @brief Declares Stout smearing class
*/
#ifndef STOUT_SMEAR_
#define STOUT_SMEAR_

namespace Grid {
namespace QCD {

/*!  @brief Stout smearing of link variable. */
template <class Gimpl>
class Smear_Stout : public Smear<Gimpl> {
 private:
  const Smear<Gimpl>* SmearBase;

 public:
  INHERIT_GIMPL_TYPES(Gimpl)

  Smear_Stout(Smear<Gimpl>* base) : SmearBase(base) {
    assert(Nc == 3);//                  "Stout smearing currently implemented only for Nc==3");
  }

  /*! Default constructor */
  Smear_Stout(double rho = 1.0) : SmearBase(new Smear_APE<Gimpl>(rho)) {
    assert(Nc == 3);//                  "Stout smearing currently implemented only for Nc==3");
  }

  ~Smear_Stout() {}  // delete SmearBase...

  void smear(GaugeField& u_smr, const GaugeField& U) const {
    GaugeField C(U._grid);
    GaugeLinkField tmp(U._grid), iq_mu(U._grid), Umu(U._grid);

    std::cout << GridLogDebug << "Stout smearing started\n";

    // Smear the configurations
    SmearBase->smear(C, U);

    for (int mu = 0; mu < Nd; mu++) {
      tmp = peekLorentz(C, mu);
      Umu = peekLorentz(U, mu);
      iq_mu = Ta(
          tmp *
          adj(Umu));  // iq_mu = Ta(Omega_mu) to match the signs with the paper
      exponentiate_iQ(tmp, iq_mu);
      pokeLorentz(u_smr, tmp * Umu, mu);  // u_smr = exp(iQ_mu)*U_mu
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

    GridBase* grid = iQ._grid;
    GaugeLinkField unity(grid);
    unity = 1.0;

    GaugeLinkField iQ2(grid), iQ3(grid);
    LatticeComplex u(grid), w(grid);
    LatticeComplex f0(grid), f1(grid), f2(grid);

    iQ2 = iQ * iQ;
    iQ3 = iQ * iQ2;

    set_uw(u, w, iQ2, iQ3);
    set_fj(f0, f1, f2, u, w);

    e_iQ = f0 * unity + timesMinusI(f1) * iQ - f2 * iQ2;
  };

  void set_uw(LatticeComplex& u, LatticeComplex& w, GaugeLinkField& iQ2,
              GaugeLinkField& iQ3) const {
    Complex one_over_three = 1.0 / 3.0;
    Complex one_over_two = 1.0 / 2.0;

    GridBase* grid = u._grid;
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
    GridBase* grid = u._grid;
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
    // Define a function to do the check
    // if( w < 1e-4 ) std::cout << GridLogWarning<< "[Smear_stout] w too small:
    // "<< w <<"\n";
    return sin(w) / w;
  }

  LatticeComplex func_xi1(const LatticeComplex& w) const {
    // Define a function to do the check
    // if( w < 1e-4 ) std::cout << GridLogWarning << "[Smear_stout] w too small:
    // "<< w <<"\n";
    return cos(w) / (w * w) - sin(w) / (w * w * w);
  }
};
}
}

#endif
