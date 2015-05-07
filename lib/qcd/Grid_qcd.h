#ifndef GRID_QCD_H
#define GRID_QCD_H
namespace Grid{

namespace QCD {

    static const int Nc=3;
    static const int Ns=4;
    static const int Nd=4;
    static const int Nhs=2; // half spinor
    static const int Nds=8; // double stored gauge field

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

    // s,sp,c,spc,lc
    template<typename vtype> using iSinglet                   = iScalar<iScalar<iScalar<vtype> > >;
    template<typename vtype> using iSpinMatrix                = iScalar<iMatrix<iScalar<vtype>, Ns> >;
    template<typename vtype> using iColourMatrix              = iScalar<iScalar<iMatrix<vtype, Nc> > > ;
    template<typename vtype> using iSpinColourMatrix          = iScalar<iMatrix<iMatrix<vtype, Nc>, Ns> >;
    template<typename vtype> using iLorentzColourMatrix       = iVector<iScalar<iMatrix<vtype, Nc> >, Nd > ;
    template<typename vtype> using iDoubleStoredColourMatrix  = iVector<iScalar<iMatrix<vtype, Nc> >, Nds > ;
    template<typename vtype> using iSpinVector                = iScalar<iVector<iScalar<vtype>, Ns> >;
    template<typename vtype> using iColourVector              = iScalar<iScalar<iVector<vtype, Nc> > >;
    template<typename vtype> using iSpinColourVector          = iScalar<iVector<iVector<vtype, Nc>, Ns> >;
    template<typename vtype> using iHalfSpinVector            = iScalar<iVector<iScalar<vtype>, Nhs> >;
    template<typename vtype> using iHalfSpinColourVector      = iScalar<iVector<iVector<vtype, Nc>, Nhs> >;

    // Spin matrix
    typedef iSpinMatrix<Complex  >          SpinMatrix;
    typedef iSpinMatrix<ComplexF >          SpinMatrixF;
    typedef iSpinMatrix<ComplexD >          SpinMatrixD;

    typedef iSpinMatrix<vComplex >          vSpinMatrix;
    typedef iSpinMatrix<vComplexF>          vSpinMatrixF;
    typedef iSpinMatrix<vComplexD>          vSpinMatrixD;

    // Colour Matrix
    typedef iColourMatrix<Complex  >        ColourMatrix;
    typedef iColourMatrix<ComplexF >        ColourMatrixF;
    typedef iColourMatrix<ComplexD >        ColourMatrixD;

    typedef iColourMatrix<vComplex >        vColourMatrix;
    typedef iColourMatrix<vComplexF>        vColourMatrixF;
    typedef iColourMatrix<vComplexD>        vColourMatrixD;

    // SpinColour matrix
    typedef iSpinColourMatrix<Complex  >    SpinColourMatrix;
    typedef iSpinColourMatrix<ComplexF >    SpinColourMatrixF;
    typedef iSpinColourMatrix<ComplexD >    SpinColourMatrixD;

    typedef iSpinColourMatrix<vComplex >    vSpinColourMatrix;
    typedef iSpinColourMatrix<vComplexF>    vSpinColourMatrixF;
    typedef iSpinColourMatrix<vComplexD>    vSpinColourMatrixD;

    // LorentzColour
    typedef iLorentzColourMatrix<Complex  > LorentzColourMatrix;
    typedef iLorentzColourMatrix<ComplexF > LorentzColourMatrixF;
    typedef iLorentzColourMatrix<ComplexD > LorentzColourMatrixD;

    typedef iLorentzColourMatrix<vComplex > vLorentzColourMatrix;
    typedef iLorentzColourMatrix<vComplexF> vLorentzColourMatrixF;
    typedef iLorentzColourMatrix<vComplexD> vLorentzColourMatrixD;

    // DoubleStored gauge field
    typedef iDoubleStoredColourMatrix<Complex  > DoubleStoredColourMatrix;
    typedef iDoubleStoredColourMatrix<ComplexF > DoubleStoredColourMatrixF;
    typedef iDoubleStoredColourMatrix<ComplexD > DoubleStoredColourMatrixD;

    typedef iDoubleStoredColourMatrix<vComplex > vDoubleStoredColourMatrix;
    typedef iDoubleStoredColourMatrix<vComplexF> vDoubleStoredColourMatrixF;
    typedef iDoubleStoredColourMatrix<vComplexD> vDoubleStoredColourMatrixD;

    // Spin vector
    typedef iSpinVector<Complex >           SpinVector;
    typedef iSpinVector<ComplexF>           SpinVectorF;
    typedef iSpinVector<ComplexD>           SpinVectorD;

    typedef iSpinVector<vComplex >           vSpinVector;
    typedef iSpinVector<vComplexF>           vSpinVectorF;
    typedef iSpinVector<vComplexD>           vSpinVectorD;

    // Colour vector
    typedef iColourVector<Complex >         ColourVector;
    typedef iColourVector<ComplexF>         ColourVectorF;
    typedef iColourVector<ComplexD>         ColourVectorD;

    typedef iColourVector<vComplex >         vColourVector;
    typedef iColourVector<vComplexF>         vColourVectorF;
    typedef iColourVector<vComplexD>         vColourVectorD;

    // SpinColourVector
    typedef iSpinColourVector<Complex >     SpinColourVector;
    typedef iSpinColourVector<ComplexF>     SpinColourVectorF;
    typedef iSpinColourVector<ComplexD>     SpinColourVectorD;

    typedef iSpinColourVector<vComplex >     vSpinColourVector;
    typedef iSpinColourVector<vComplexF>     vSpinColourVectorF;
    typedef iSpinColourVector<vComplexD>     vSpinColourVectorD;

    // HalfSpin vector
    typedef iHalfSpinVector<Complex >       HalfSpinVector;
    typedef iHalfSpinVector<ComplexF>       HalfSpinVectorF;
    typedef iHalfSpinVector<ComplexD>       HalfSpinVectorD;

    typedef iHalfSpinVector<vComplex >       vHalfSpinVector;
    typedef iHalfSpinVector<vComplexF>       vHalfSpinVectorF;
    typedef iHalfSpinVector<vComplexD>       vHalfSpinVectorD;

    // HalfSpinColour vector
    typedef iHalfSpinColourVector<Complex > HalfSpinColourVector;
    typedef iHalfSpinColourVector<ComplexF> HalfSpinColourVectorF;
    typedef iHalfSpinColourVector<ComplexD> HalfSpinColourVectorD;
    
    typedef iHalfSpinColourVector<vComplex > vHalfSpinColourVector;
    typedef iHalfSpinColourVector<vComplexF> vHalfSpinColourVectorF;
    typedef iHalfSpinColourVector<vComplexD> vHalfSpinColourVectorD;
    
    // singlets
    typedef iSinglet<Complex >         TComplex;    // FIXME This is painful. Tensor singlet complex type.
    typedef iSinglet<ComplexF>         TComplexF;    // FIXME This is painful. Tensor singlet complex type.
    typedef iSinglet<ComplexD>         TComplexD;    // FIXME This is painful. Tensor singlet complex type.

    typedef iSinglet<vComplex >        vTComplex ;   // what if we don't know the tensor structure
    typedef iSinglet<vComplexF>        vTComplexF;   // what if we don't know the tensor structure
    typedef iSinglet<vComplexD>        vTComplexD;   // what if we don't know the tensor structure

    typedef iSinglet<Real >            TReal;       // Shouldn't need these; can I make it work without?
    typedef iSinglet<RealF>            TRealF;       // Shouldn't need these; can I make it work without?
    typedef iSinglet<RealD>            TRealD;       // Shouldn't need these; can I make it work without?

    typedef iSinglet<vReal >           vTReal;      
    typedef iSinglet<vRealF>           vTRealF;      
    typedef iSinglet<vRealD>           vTRealD;      

    typedef iSinglet<vInteger>         vTInteger;
    typedef iSinglet<Integer >         TInteger;


    // Lattices of these
    typedef Lattice<vColourMatrix>          LatticeColourMatrix;
    typedef Lattice<vColourMatrixF>         LatticeColourMatrixF;
    typedef Lattice<vColourMatrixD>         LatticeColourMatrixD;

    typedef Lattice<vSpinMatrix>            LatticeSpinMatrix;
    typedef Lattice<vSpinMatrixF>           LatticeSpinMatrixF;
    typedef Lattice<vSpinMatrixD>           LatticeSpinMatrixD;

    typedef Lattice<vSpinColourMatrix>      LatticeSpinColourMatrix;
    typedef Lattice<vSpinColourMatrixF>     LatticeSpinColourMatrixF;
    typedef Lattice<vSpinColourMatrixD>     LatticeSpinColourMatrixD;


    typedef Lattice<vLorentzColourMatrix>  LatticeLorentzColourMatrix;
    typedef Lattice<vLorentzColourMatrixF> LatticeLorentzColourMatrixF;
    typedef Lattice<vLorentzColourMatrixD> LatticeLorentzColourMatrixD;

    // DoubleStored gauge field
    typedef Lattice<vDoubleStoredColourMatrix>  LatticeDoubleStoredColourMatrix;
    typedef Lattice<vDoubleStoredColourMatrixF> LatticeDoubleStoredColourMatrixF;
    typedef Lattice<vDoubleStoredColourMatrixD> LatticeDoubleStoredColourMatrixD;

    typedef Lattice<vSpinVector>            LatticeSpinVector;
    typedef Lattice<vSpinVectorF>           LatticeSpinVectorF;
    typedef Lattice<vSpinVectorD>           LatticeSpinVectorD;

    typedef Lattice<vColourVector>          LatticeColourVector;
    typedef Lattice<vColourVectorF>         LatticeColourVectorF;
    typedef Lattice<vColourVectorD>         LatticeColourVectorD;

    typedef Lattice<vSpinColourVector>      LatticeSpinColourVector;
    typedef Lattice<vSpinColourVectorF>     LatticeSpinColourVectorF;
    typedef Lattice<vSpinColourVectorD>     LatticeSpinColourVectorD;

    typedef Lattice<vHalfSpinVector>        LatticeHalfSpinVector;
    typedef Lattice<vHalfSpinVectorF>       LatticeHalfSpinVectorF;
    typedef Lattice<vHalfSpinVectorD>       LatticeHalfSpinVectorD;

    typedef Lattice<vHalfSpinColourVector>  LatticeHalfSpinColourVector;
    typedef Lattice<vHalfSpinColourVectorF> LatticeHalfSpinColourVectorF;
    typedef Lattice<vHalfSpinColourVectorD> LatticeHalfSpinColourVectorD;

    typedef Lattice<vTReal>            LatticeReal;
    typedef Lattice<vTRealF>           LatticeRealF;
    typedef Lattice<vTRealD>           LatticeRealD;

    typedef Lattice<vTComplex>         LatticeComplex;
    typedef Lattice<vTComplexF>        LatticeComplexF;
    typedef Lattice<vTComplexD>        LatticeComplexD;

    typedef Lattice<vTInteger>         LatticeInteger; // Predicates for "where"


    ///////////////////////////////////////////
    // Physical names for things
    ///////////////////////////////////////////
    typedef LatticeHalfSpinColourVector  LatticeHalfFermion;
    typedef LatticeHalfSpinColourVectorF LatticeHalfFermionF;
    typedef LatticeHalfSpinColourVectorF LatticeHalfFermionD;

    typedef LatticeSpinColourVector      LatticeFermion;
    typedef LatticeSpinColourVectorF     LatticeFermionF;
    typedef LatticeSpinColourVectorD     LatticeFermionD;

    typedef LatticeSpinColourMatrix                LatticePropagator;
    typedef LatticeSpinColourMatrixF               LatticePropagatorF;
    typedef LatticeSpinColourMatrixD               LatticePropagatorD;

    typedef LatticeLorentzColourMatrix             LatticeGaugeField;
    typedef LatticeLorentzColourMatrixF            LatticeGaugeFieldF;
    typedef LatticeLorentzColourMatrixD            LatticeGaugeFieldD;

    typedef LatticeDoubleStoredColourMatrix        LatticeDoubledGaugeField;
    typedef LatticeDoubleStoredColourMatrixF       LatticeDoubledGaugeFieldF;
    typedef LatticeDoubleStoredColourMatrixD       LatticeDoubledGaugeFieldD;

    // Uhgg... typing this hurt  ;)
    // (my keyboard got burning hot when I typed this, must be the anti-Fermion)
    typedef Lattice<vColourVector>          LatticeStaggeredFermion;    
    typedef Lattice<vColourVectorF>         LatticeStaggeredFermionF;    
    typedef Lattice<vColourVectorD>         LatticeStaggeredFermionD;    

    typedef Lattice<vColourMatrix>          LatticeStaggeredPropagator; 
    typedef Lattice<vColourMatrixF>         LatticeStaggeredPropagatorF; 
    typedef Lattice<vColourMatrixD>         LatticeStaggeredPropagatorD; 

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

    // FIXME transpose Colour, transpose Spin, traceColour traceSpin

}   //namespace QCD
} // Grid

#include <qcd/Grid_qcd_dirac.h>
#include <qcd/Grid_qcd_2spinor.h>
//#include <qcd/Grid_qcd_pauli.h>
#include <qcd/Grid_qcd_wilson_dop.h>

#endif
