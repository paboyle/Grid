    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/QCD.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: neo <cossu@post.kek.jp>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#ifndef GRID_QCD_H
#define GRID_QCD_H
namespace Grid{

namespace QCD {


    static const int Xp = 0;
    static const int Yp = 1;
    static const int Zp = 2;
    static const int Tp = 3;
    static const int Xm = 4;
    static const int Ym = 5;
    static const int Zm = 6;
    static const int Tm = 7;

    static const int Nc=3;
    static const int Ns=4;
    static const int Nd=4;
    static const int Nhs=2; // half spinor
    static const int Nds=8; // double stored gauge field
    static const int Ngp=2; // gparity index range

    //////////////////////////////////////////////////////////////////////////////
    // QCD iMatrix types
    // Index conventions:                            Lorentz x Spin x Colour
    //////////////////////////////////////////////////////////////////////////////
    static const int ColourIndex = 2;
    static const int SpinIndex   = 1;
    static const int LorentzIndex= 0;

    // Useful traits is this a spin index
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;

    const int SpinorIndex = 2;
    template<typename T> struct isSpinor {
      static const bool value = (SpinorIndex==T::TensorLevel);
    };
    template <typename T> using IfSpinor    = Invoke<std::enable_if< isSpinor<T>::value,int> > ;
    template <typename T> using IfNotSpinor = Invoke<std::enable_if<!isSpinor<T>::value,int> > ;

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

    template<typename vtype> using iGparitySpinColourVector       = iVector<iVector<iVector<vtype, Nc>, Ns>, Ngp >;
    template<typename vtype> using iGparityHalfSpinColourVector   = iVector<iVector<iVector<vtype, Nc>, Nhs>, Ngp >;

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
    typedef iSinglet<Complex >         TComplex;     // FIXME This is painful. Tensor singlet complex type.
    typedef iSinglet<ComplexF>         TComplexF;    // FIXME This is painful. Tensor singlet complex type.
    typedef iSinglet<ComplexD>         TComplexD;    // FIXME This is painful. Tensor singlet complex type.

    typedef iSinglet<vComplex >        vTComplex ;   // what if we don't know the tensor structure
    typedef iSinglet<vComplexF>        vTComplexF;   // what if we don't know the tensor structure
    typedef iSinglet<vComplexD>        vTComplexD;   // what if we don't know the tensor structure

    typedef iSinglet<Real >            TReal;        // Shouldn't need these; can I make it work without?
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

    template<class GF> using LorentzScalar = Lattice<iScalar<typename GF::vector_object::element> >;

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
    template<class vobj> auto peekSpin(const vobj &rhs,int i) -> decltype(PeekIndex<SpinIndex>(rhs,0))
    {
      return PeekIndex<SpinIndex>(rhs,i);
    }
    template<class vobj> auto peekSpin(const vobj &rhs,int i,int j) -> decltype(PeekIndex<SpinIndex>(rhs,0,0))
    {
      return PeekIndex<SpinIndex>(rhs,i,j);
    }
    template<class vobj> auto peekSpin(const Lattice<vobj> &rhs,int i) -> decltype(PeekIndex<SpinIndex>(rhs,0))
    {
      return PeekIndex<SpinIndex>(rhs,i);
    }
    template<class vobj> auto peekSpin(const Lattice<vobj> &rhs,int i,int j) -> decltype(PeekIndex<SpinIndex>(rhs,0,0))
    {
      return PeekIndex<SpinIndex>(rhs,i,j);
    }
    //colour
    template<class vobj> auto peekColour(const vobj &rhs,int i) -> decltype(PeekIndex<ColourIndex>(rhs,0))
    {
      return PeekIndex<ColourIndex>(rhs,i);
    }
    template<class vobj> auto peekColour(const vobj &rhs,int i,int j) -> decltype(PeekIndex<ColourIndex>(rhs,0,0))
    {
      return PeekIndex<ColourIndex>(rhs,i,j);
    }
    template<class vobj> auto peekColour(const Lattice<vobj> &rhs,int i) -> decltype(PeekIndex<ColourIndex>(rhs,0))
    {
      return PeekIndex<ColourIndex>(rhs,i);
    }
    template<class vobj> auto peekColour(const Lattice<vobj> &rhs,int i,int j) -> decltype(PeekIndex<ColourIndex>(rhs,0,0))
    {
      return PeekIndex<ColourIndex>(rhs,i,j);
    }
    //lorentz
    template<class vobj> auto peekLorentz(const vobj &rhs,int i) -> decltype(PeekIndex<LorentzIndex>(rhs,0))
    {
      return PeekIndex<LorentzIndex>(rhs,i);
    }
    template<class vobj> auto peekLorentz(const Lattice<vobj> &rhs,int i) -> decltype(PeekIndex<LorentzIndex>(rhs,0))
    {
      return PeekIndex<LorentzIndex>(rhs,i);
    }

    //////////////////////////////////////////////
    // Poke lattice
    //////////////////////////////////////////////
    template<class vobj> 
      void pokeColour(Lattice<vobj> &lhs,
		      const Lattice<decltype(peekIndex<ColourIndex>(lhs._odata[0],0))> & rhs,
		      int i)
    {
      PokeIndex<ColourIndex>(lhs,rhs,i);
    }
    template<class vobj> 
      void pokeColour(Lattice<vobj> &lhs,
		      const Lattice<decltype(peekIndex<ColourIndex>(lhs._odata[0],0,0))> & rhs,
		      int i,int j)
    {
      PokeIndex<ColourIndex>(lhs,rhs,i,j);
    }
    template<class vobj> 
      void pokeSpin(Lattice<vobj> &lhs,
		      const Lattice<decltype(peekIndex<SpinIndex>(lhs._odata[0],0))> & rhs,
		      int i)
    {
      PokeIndex<SpinIndex>(lhs,rhs,i);
    }
    template<class vobj> 
      void pokeSpin(Lattice<vobj> &lhs,
		      const Lattice<decltype(peekIndex<SpinIndex>(lhs._odata[0],0,0))> & rhs,
		      int i,int j)
    {
      PokeIndex<SpinIndex>(lhs,rhs,i,j);
    }
    template<class vobj> 
      void pokeLorentz(Lattice<vobj> &lhs,
		      const Lattice<decltype(peekIndex<LorentzIndex>(lhs._odata[0],0))> & rhs,
		      int i)
    {
      PokeIndex<LorentzIndex>(lhs,rhs,i);
    }

    //////////////////////////////////////////////
    // Poke scalars
    //////////////////////////////////////////////
    template<class vobj> void pokeSpin(vobj &lhs,const decltype(peekIndex<SpinIndex>(lhs,0)) & rhs,int i)
    {
      pokeIndex<SpinIndex>(lhs,rhs,i);
    }
    template<class vobj> void pokeSpin(vobj &lhs,const decltype(peekIndex<SpinIndex>(lhs,0,0)) & rhs,int i,int j)
    {
      pokeIndex<SpinIndex>(lhs,rhs,i,j);
    }

    template<class vobj> void pokeColour(vobj &lhs,const decltype(peekIndex<ColourIndex>(lhs,0)) & rhs,int i)
    {
      pokeIndex<ColourIndex>(lhs,rhs,i);
    }
    template<class vobj> void pokeColour(vobj &lhs,const decltype(peekIndex<ColourIndex>(lhs,0,0)) & rhs,int i,int j)
    {
      pokeIndex<ColourIndex>(lhs,rhs,i,j);
    }

    template<class vobj> void pokeLorentz(vobj &lhs,const decltype(peekIndex<LorentzIndex>(lhs,0)) & rhs,int i)
    {
      pokeIndex<LorentzIndex>(lhs,rhs,i);
    }

    //////////////////////////////////////////////
    // Fermion <-> propagator assignements
    //////////////////////////////////////////////
    template <class Prop, class Ferm>
    void FermToProp(Prop &p, const Ferm &f, const int s, const int c)
    {
        for(int j = 0; j < Ns; ++j)
        {
            auto pjs = peekSpin(p, j, s);
            auto fj  = peekSpin(f, j);
            
            for(int i = 0; i < Nc; ++i)
            {
                pokeColour(pjs, peekColour(fj, i), i, c);
            }
            pokeSpin(p, pjs, j, s);
        }
    }
    
    template <class Prop, class Ferm>
    void PropToFerm(Ferm &f, const Prop &p, const int s, const int c)
    {
        for(int j = 0; j < Ns; ++j)
        {
            auto pjs = peekSpin(p, j, s);
            auto fj  = peekSpin(f, j);
            
            for(int i = 0; i < Nc; ++i)
            {
                pokeColour(fj, peekColour(pjs, i, c), i);
            }
            pokeSpin(f, fj, j);
        }
    }
    
    //////////////////////////////////////////////
    // transpose array and scalar
    //////////////////////////////////////////////
    template<int Index,class vobj> inline Lattice<vobj> transposeSpin(const Lattice<vobj> &lhs){
      return transposeIndex<SpinIndex>(lhs);
    }
    template<int Index,class vobj> inline Lattice<vobj> transposeColour(const Lattice<vobj> &lhs){
      return transposeIndex<ColourIndex>(lhs);
    }
    template<int Index,class vobj> inline vobj transposeSpin(const vobj &lhs){
      return transposeIndex<SpinIndex>(lhs);
    }
    template<int Index,class vobj> inline vobj transposeColour(const vobj &lhs){
      return transposeIndex<ColourIndex>(lhs);
    }

    //////////////////////////////////////////
    // Trace lattice and non-lattice
    //////////////////////////////////////////
    template<int Index,class vobj>
    inline auto traceSpin(const Lattice<vobj> &lhs) -> Lattice<decltype(traceIndex<SpinIndex>(lhs._odata[0]))>
    {
      return traceIndex<SpinIndex>(lhs);
    }
    template<int Index,class vobj>
    inline auto traceColour(const Lattice<vobj> &lhs) -> Lattice<decltype(traceIndex<ColourIndex>(lhs._odata[0]))>
    {
      return traceIndex<ColourIndex>(lhs);
    }
    template<int Index,class vobj>
    inline auto traceSpin(const vobj &lhs) -> Lattice<decltype(traceIndex<SpinIndex>(lhs))>
    {
      return traceIndex<SpinIndex>(lhs);
    }
    template<int Index,class vobj>
    inline auto traceColour(const vobj &lhs) -> Lattice<decltype(traceIndex<ColourIndex>(lhs))>
    {
      return traceIndex<ColourIndex>(lhs);
    }

}   //namespace QCD
} // Grid

#include <qcd/utils/SpaceTimeGrid.h>
#include <qcd/spin/Dirac.h>
#include <qcd/spin/TwoSpinor.h>
#include <qcd/utils/LinalgUtils.h>
#include <qcd/utils/CovariantCshift.h>
#include <qcd/utils/SUn.h>
#include <qcd/action/Actions.h>
#include <qcd/hmc/integrators/Integrator.h>
#include <qcd/hmc/integrators/Integrator_algorithm.h>
#include <qcd/hmc/HMC.h>


#endif
