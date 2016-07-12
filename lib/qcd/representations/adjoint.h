/*
 *	Policy classes for the HMC
 *	Author: Guido Cossu
*/	

#ifndef ADJOINT_H
#define ADJOINT_H


namespace Grid {
namespace QCD {

/*
* This is an helper class for the HMC
* Should contain only the data for the adjoint representation
* and the facility to convert from the fundamental -> adjoint
*/

template <int ncolour>
class AdjointRep {
 public:
  typename SU_Adjoint<ncolour>::LatticeAdjMatrix U;
  const int Dimension = ncolour * ncolour - 1;

  explicit AdjointRep(GridBase* grid):U(grid) {}
  void update_representation(const LatticeGaugeField& Uin) {
    // Uin is in the fundamental representation
    // get the U in AdjointRep
    // (U_adj)_B = tr[e^a U e^b U^dag]
    // e^a = t^a/sqrt(T_F)
    // where t^a is the generator in the fundamental
    // T_F is 1/2 for the fundamental representation
    conformable(U, Uin);
    U = zero;
    LatticeGaugeField tmp(Uin._grid);

    Vector<typename SU<ncolour>::Matrix > ta(ncolour * ncolour - 1);

    // FIXME probably not very efficient to get all the generators everytime
    for (int a = 0; a < Dimension; a++) SU<ncolour>::generator(a, ta[a]);

    for (int a = 0; a < Dimension; a++) {
    	tmp = 2.0 * adj(Uin) * ta[a] * Uin;
      for (int b = 0; b < (ncolour * ncolour - 1); b++) {
        auto Tr = TensorRemove(trace(tmp * ta[b]));
        pokeColour(U, Tr, a,b);
      }
    }  	 

  }
};

typedef	 AdjointRep<Nc> AdjointRepresentation;

}
}




#endif