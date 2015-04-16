#ifndef GRID_QCD_H
#define GRID_QCD_H
namespace Grid{
namespace QCD {

    static const int Nc=3;
    static const int Ns=4;
    static const int CbRed  =0;
    static const int CbBlack=1;
    
    // QCD iMatrix types
    template<typename vtype> using iSinglet          = iScalar<iScalar<vtype> > ;
    template<typename vtype> using iSpinMatrix       = iMatrix<iScalar<vtype>, Ns>;
    template<typename vtype> using iSpinColourMatrix = iMatrix<iMatrix<vtype, Nc>, Ns>;
    template<typename vtype> using iColourMatrix     = iScalar<iMatrix<vtype, Nc>> ;

    template<typename vtype> using iSpinVector       = iVector<iScalar<vtype>, Ns>;
    template<typename vtype> using iColourVector     = iScalar<iVector<vtype, Nc> >;
    template<typename vtype> using iSpinColourVector = iVector<iVector<vtype, Nc>, Ns>;

    typedef iSinglet<Complex >          TComplex;    // This is painful. Tensor singlet complex type.
    typedef iSinglet<vComplex >         vTComplex;   // what if we don't know the tensor structure
    typedef iSinglet<Real >             TReal;       // Shouldn't need these.
    typedef iSinglet<vInteger >         vTInteger;

    typedef iSpinMatrix<Complex >       SpinMatrix;
    typedef iColourMatrix<Complex >     ColourMatrix;
    typedef iSpinColourMatrix<Complex > SpinColourMatrix;

    typedef iSpinVector<Complex >       SpinVector;
    typedef iColourVector<Complex >     ColourVector;
    typedef iSpinColourVector<Complex > SpinColourVector;

    
    typedef iSpinMatrix<vComplex >       vSpinMatrix;
    typedef iColourMatrix<vComplex >     vColourMatrix;
    typedef iSpinColourMatrix<vComplex > vSpinColourMatrix;
    
    typedef iSpinVector<vComplex >       vSpinVector;
    typedef iColourVector<vComplex >     vColourVector;
    typedef iSpinColourVector<vComplex > vSpinColourVector;
    
    typedef Lattice<vTComplex>            LatticeComplex;
    typedef Lattice<vInteger>            LatticeInteger; // Predicates for "where"
    
    typedef Lattice<vColourMatrix>     LatticeColourMatrix;
    typedef Lattice<vSpinMatrix>       LatticeSpinMatrix;
    typedef Lattice<vSpinColourMatrix> LatticePropagator;
    typedef LatticePropagator LatticeSpinColourMatrix;

    typedef Lattice<vSpinColourVector> LatticeFermion;
    typedef Lattice<vSpinColourVector> LatticeSpinColourVector;
    typedef Lattice<vSpinVector>       LatticeSpinVector;
    typedef Lattice<vColourVector>     LatticeColourVector;


    // FIXME for debug; deprecate this
   inline void LatticeCoordinate(LatticeInteger &l,int mu){
      GridBase *grid = l._grid;
      int Nsimd = grid->iSites();
      std::vector<int> gcoor;
      std::vector<Integer> mergebuf(Nsimd);
      std::vector<Integer *> mergeptr(Nsimd);
      for(int o=0;o<grid->oSites();o++){
	for(int i=0;i<grid->iSites();i++){
	  //	  RankIndexToGlobalCoor(grid->ThisRank(),o,i,gcoor);
	  grid->RankIndexToGlobalCoor(0,o,i,gcoor);
	  mergebuf[i]=gcoor[mu];
	  mergeptr[i]=&mergebuf[i];
	}
	merge(l._odata[o],mergeptr);
      }
    };

#include <Grid_predicated.h>

#if 0

#endif

}   //namespace QCD
} // Grid
#endif
