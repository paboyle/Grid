#include <Grid.h>

namespace Grid {
namespace QCD {

    /*
     * BF sequence
     *
      void bfmbase<Float>::MooeeInv(Fermion_t psi, 
			       Fermion_t chi, 
			      int dag, int cb)

    double m    = this->mass;
    double tm   = this->twistedmass;
    double mtil = 4.0+this->mass;

    double sq = mtil*mtil + tm*tm;

    double a = mtil/sq;
    double b = -tm /sq;
    if(dag) b=-b;
    axpibg5x(chi,psi,a,b);

      void bfmbase<Float>::Mooee(Fermion_t psi, 
			   Fermion_t chi, 
			   int dag,int cb)
    double a = 4.0+this->mass;
    double b = this->twistedmass;
    if(dag) b=-b;
    axpibg5x(chi,psi,a,b);
    */

  template<class Impl>
  void WilsonTMFermion<Impl>::Mooee(const FermionField &in, FermionField &out) {
    RealD a = 4.0+this->mass;
    RealD b = this->mu;
    out.checkerboard = in.checkerboard;
    axpibg5x(out,in,a,b);
  }
  template<class Impl>
  void WilsonTMFermion<Impl>::MooeeDag(const FermionField &in, FermionField &out) {
    RealD a = 4.0+this->mass;
    RealD b = -this->mu;
    out.checkerboard = in.checkerboard;
    axpibg5x(out,in,a,b);
  }
  template<class Impl>
  void WilsonTMFermion<Impl>::MooeeInv(const FermionField &in, FermionField &out) {
    RealD m    = this->mass;
    RealD tm   = this->mu;
    RealD mtil = 4.0+this->mass;
    RealD sq   = mtil*mtil+tm*tm;
    RealD a    = mtil/sq;
    RealD b    = -tm /sq;
    axpibg5x(out,in,a,b);
  }
  template<class Impl>
  void WilsonTMFermion<Impl>::MooeeInvDag(const FermionField &in, FermionField &out) {
    RealD m    = this->mass;
    RealD tm   = this->mu;
    RealD mtil = 4.0+this->mass;
    RealD sq   = mtil*mtil+tm*tm;
    RealD a    = mtil/sq;
    RealD b    = tm /sq;
    axpibg5x(out,in,a,b);
  }

  FermOpTemplateInstantiate(WilsonTMFermion);

}
}
