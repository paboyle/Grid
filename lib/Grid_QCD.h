#ifndef GRID_QCD_H
#define GRID_QCD_H
namespace Grid{

namespace QCD {

    static const int Nc=3;
    static const int Ns=4;
    static const int Nd=4;

    static const int CbRed  =0;
    static const int CbBlack=1;
    
    //////////////////////////////////////////////////////////////////////////////
    // QCD iMatrix types
    // Index conventions:                            Lorentz x Spin x Colour
    //////////////////////////////////////////////////////////////////////////////
    static const int ColourIndex = 1;
    static const int SpinIndex   = 2;
    static const int LorentzIndex= 3;
    

    // ChrisK very keen to add extra space for Gparity doubling.
    //
    // Also add domain wall index, in a way where Wilson operator 
    // naturally distributes across the 5th dimensions.
    //
    // That probably makes for GridRedBlack4dCartesian grid.

    template<typename vtype> using iSinglet          = iScalar<iScalar<iScalar<vtype> > >;
    template<typename vtype> using iSpinMatrix       = iScalar<iMatrix<iScalar<vtype>, Ns> >;
    template<typename vtype> using iSpinColourMatrix = iScalar<iMatrix<iMatrix<vtype, Nc>, Ns> >;
    template<typename vtype> using iColourMatrix     = iScalar<iScalar<iMatrix<vtype, Nc> > > ;
    template<typename vtype> using iLorentzColourMatrix = iVector<iScalar<iMatrix<vtype, Nc> >, Nd > ;


    template<typename vtype> using iSpinVector       = iScalar<iVector<iScalar<vtype>, Ns> >;
    template<typename vtype> using iColourVector     = iScalar<iScalar<iVector<vtype, Nc> > >;
    template<typename vtype> using iSpinColourVector = iScalar<iVector<iVector<vtype, Nc>, Ns> >;

    typedef iSpinMatrix<Complex >          SpinMatrix;
    typedef iColourMatrix<Complex >        ColourMatrix;
    typedef iSpinColourMatrix<Complex >    SpinColourMatrix;
    typedef iLorentzColourMatrix<Complex > LorentzColourMatrix;

    typedef iSpinVector<Complex >       SpinVector;
    typedef iColourVector<Complex >     ColourVector;
    typedef iSpinColourVector<Complex > SpinColourVector;

    
    typedef iSpinMatrix<vComplex >          vSpinMatrix;
    typedef iColourMatrix<vComplex >        vColourMatrix;
    typedef iSpinColourMatrix<vComplex >    vSpinColourMatrix;
    typedef iLorentzColourMatrix<vComplex > vLorentzColourMatrix;
    
    typedef iSpinVector<vComplex >       vSpinVector;
    typedef iColourVector<vComplex >     vColourVector;
    typedef iSpinColourVector<vComplex > vSpinColourVector;
    
    typedef iSinglet<Complex >          TComplex;    // This is painful. Tensor singlet complex type.
    typedef iSinglet<vComplex >         vTComplex;   // what if we don't know the tensor structure
    typedef iSinglet<Real >             TReal;       // Shouldn't need these; can I make it work without?
    typedef iSinglet<vReal >            vTReal;      
    typedef iSinglet<vInteger >         vTInteger;
    typedef iSinglet<Integer >          TInteger;

    typedef Lattice<vTReal>              LatticeReal;
    typedef Lattice<vTComplex>           LatticeComplex;
    typedef Lattice<vInteger>            LatticeInteger; // Predicates for "where"
    
    typedef Lattice<vColourMatrix>     LatticeColourMatrix;
    typedef Lattice<vSpinMatrix>       LatticeSpinMatrix;
    typedef Lattice<vSpinColourMatrix> LatticeSpinColourMatrix;

    typedef Lattice<vSpinColourVector> LatticeSpinColourVector;
    typedef Lattice<vSpinVector>       LatticeSpinVector;
    typedef Lattice<vColourVector>     LatticeColourVector;

    ///////////////////////////////////////////
    // Physical names for things
    ///////////////////////////////////////////
    typedef Lattice<vSpinColourVector> LatticeFermion;
    typedef Lattice<vSpinColourMatrix> LatticePropagator;
    typedef Lattice<vLorentzColourMatrix> LatticeGaugeField;


    //////////////////////////////////////////////////////////////////////////////
    // Peek and Poke named after physics attributes
    //////////////////////////////////////////////////////////////////////////////
    //spin
    template<class vobj> auto peekSpin(const vobj &rhs,int i) -> decltype(peekIndex<SpinIndex>(rhs,0))
    {
      return peekIndex<SpinIndex>(rhs,i);
    }
    template<class vobj> auto peekSpin(const vobj &rhs,int i,int j) -> decltype(peekIndex<SpinIndex>(rhs,0,0))
    {
      return peekIndex<SpinIndex>(rhs,i,j);
    }
    template<class vobj> auto peekSpin(const Lattice<vobj> &rhs,int i) -> decltype(peekIndex<SpinIndex>(rhs,0))
    {
      return peekIndex<SpinIndex>(rhs,i);
    }
    template<class vobj> auto peekSpin(const Lattice<vobj> &rhs,int i,int j) -> decltype(peekIndex<SpinIndex>(rhs,0,0))
    {
      return peekIndex<SpinIndex>(rhs,i,j);
    }
    //colour
    template<class vobj> auto peekColour(const vobj &rhs,int i) -> decltype(peekIndex<ColourIndex>(rhs,0))
    {
      return peekIndex<ColourIndex>(rhs,i);
    }
    template<class vobj> auto peekColour(const vobj &rhs,int i,int j) -> decltype(peekIndex<ColourIndex>(rhs,0,0))
    {
      return peekIndex<ColourIndex>(rhs,i,j);
    }
    template<class vobj> auto peekColour(const Lattice<vobj> &rhs,int i) -> decltype(peekIndex<ColourIndex>(rhs,0))
    {
      return peekIndex<ColourIndex>(rhs,i);
    }
    template<class vobj> auto peekColour(const Lattice<vobj> &rhs,int i,int j) -> decltype(peekIndex<ColourIndex>(rhs,0,0))
    {
      return peekIndex<ColourIndex>(rhs,i,j);
    }
    //lorentz
    template<class vobj> auto peekLorentz(const vobj &rhs,int i) -> decltype(peekIndex<LorentzIndex>(rhs,0))
    {
      return peekIndex<LorentzIndex>(rhs,i);
    }
    template<class vobj> auto peekLorentz(const vobj &rhs,int i,int j) -> decltype(peekIndex<LorentzIndex>(rhs,0,0))
    {
      return peekIndex<LorentzIndex>(rhs,i,j);
    }
    template<class vobj> auto peekLorentz(const Lattice<vobj> &rhs,int i) -> decltype(peekIndex<LorentzIndex>(rhs,0))
    {
      return peekIndex<LorentzIndex>(rhs,i);
    }
    template<class vobj> auto peekLorentz(const Lattice<vobj> &rhs,int i,int j) -> decltype(peekIndex<LorentzIndex>(rhs,0,0))
    {
      return peekIndex<LorentzIndex>(rhs,i,j);
    }

}   //namespace QCD
} // Grid
#endif
