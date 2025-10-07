#ifndef GRID_QCD_PAULI_H
#define GRID_QCD_PAULI_H

#include <array>

NAMESPACE_BEGIN(Grid);
//
/*
 * Pauli basis
 * sx        sy       sz       ident
 * (0 1)  , (0 -i) , ( 1 0 )
 * (1 0)    (i  0)   ( 0 -1)
 *
 * These are hermitian.
 *
 * Also supply wilson "projectors" (1+/-sx), (1+/-sy), (1+/-sz)
 *
 * spPauliProjXm
 * spPauliProjYm etc...
 */
class Pauli {
  public:
    GRID_SERIALIZABLE_ENUM(Algebra, undef,
                           SigmaX           , 0,
                           MinusSigmaX      , 1,
                           SigmaY           , 2,
                           MinusSigmaY      , 3,
                           SigmaZ           , 4,
                           MinusSigmaZ      , 5,
                           Identity         , 6,
			   MinusIdentity    , 7);
  
    static constexpr unsigned int nPauli = 8;
    static const std::array<const char *, nPauli>                name;
    static const std::array<std::array<Algebra, nPauli>, nPauli> mul;
    static const std::array<Algebra, nPauli>                     adj;
    static const std::array<const Pauli, 4>                      gmu;
    static const std::array<const Pauli, 16>                     gall;
    Algebra                                                      g;
  public:
  accelerator Pauli(Algebra initg): g(initg) {}  
};

#define CopyImplementation(iTemplate,multPauli,multFlavour)	\
  template<class vtype>							\
  accelerator_inline void multPauli(iTemplate<vtype, Nhs> &ret, const iTemplate<vtype, Nhs> &rhs) {	\
    multFlavour(ret,rhs);						\
}

CopyImplementation(iVector,multPauliSigmaX,multFlavourSigmaX);
CopyImplementation(iMatrix,lmultPauliSigmaX,lmultFlavourSigmaX);
CopyImplementation(iMatrix,rmultPauliSigmaX,rmultFlavourSigmaX);

CopyImplementation(iVector,multPauliMinusSigmaX ,multFlavourMinusSigmaX);
CopyImplementation(iMatrix,lmultPauliMinusSigmaX,lmultFlavourMinusSigmaX);
CopyImplementation(iMatrix,rmultPauliMinusSigmaX,rmultFlavourMinusSigmaX);

CopyImplementation(iVector,multPauliSigmaY,multFlavourSigmaY);
CopyImplementation(iMatrix,lmultPauliSigmaY,lmultFlavourSigmaY);
CopyImplementation(iMatrix,rmultPauliSigmaY,rmultFlavourSigmaY);

CopyImplementation(iVector,multPauliMinusSigmaY ,multFlavourMinusSigmaY);
CopyImplementation(iMatrix,lmultPauliMinusSigmaY,lmultFlavourMinusSigmaY);
CopyImplementation(iMatrix,rmultPauliMinusSigmaY,rmultFlavourMinusSigmaY);

CopyImplementation(iVector,multPauliSigmaZ,multFlavourSigmaZ);
CopyImplementation(iMatrix,lmultPauliSigmaZ,lmultFlavourSigmaZ);
CopyImplementation(iMatrix,rmultPauliSigmaZ,rmultFlavourSigmaZ);

CopyImplementation(iVector,multPauliMinusSigmaZ ,multFlavourMinusSigmaZ);
CopyImplementation(iMatrix,lmultPauliMinusSigmaZ,lmultFlavourMinusSigmaZ);
CopyImplementation(iMatrix,rmultPauliMinusSigmaZ,rmultFlavourMinusSigmaZ);

CopyImplementation(iVector,multPauliIdentity,multFlavourIdentity);
CopyImplementation(iMatrix,lmultPauliIdentity,lmultFlavourIdentity);
CopyImplementation(iMatrix,rmultPauliIdentity,rmultFlavourIdentity);

CopyImplementation(iVector,multPauliMinusIdentity ,multFlavourMinusIdentity);
CopyImplementation(iMatrix,lmultPauliMinusIdentity,lmultFlavourMinusIdentity);
CopyImplementation(iMatrix,rmultPauliMinusIdentity,rmultFlavourMinusIdentity);

/*
 * sx        sy       sz       ident
 * (0 1)  , (0 -i) , ( 1 0 )
 * (1 0)    (i  0)   ( 0 -1)
 */
template<class vtype,IfSpinor<iVector<vtype,Nhs> > = 0> accelerator_inline void pauliProjXp (iVector<vtype,Nhs> &hspin,const iVector<vtype,Nhs> &fspin)
{
  hspin(0)=fspin(0)+fspin(1);
  hspin(1)=fspin(1)+fspin(0);
}
template<class vtype,IfSpinor<iVector<vtype,Nhs> > = 0> accelerator_inline void pauliProjXm (iVector<vtype,Nhs> &hspin,const iVector<vtype,Nhs> &fspin)
{
  hspin(0)=fspin(0)-fspin(1);
  hspin(1)=fspin(1)-fspin(0);
}

template<class vtype,IfSpinor<iVector<vtype,Nhs> > = 0> accelerator_inline void pauliProjYp (iVector<vtype,Nhs> &hspin,const iVector<vtype,Nhs> &fspin)
{
  hspin(0)=fspin(0)-timesI(fspin(1));
  hspin(1)=fspin(1)+timesI(fspin(0));
}
template<class vtype,IfSpinor<iVector<vtype,Nhs> > = 0> accelerator_inline void pauliProjYm (iVector<vtype,Nhs> &hspin,const iVector<vtype,Nhs> &fspin)
{
  hspin(0)=fspin(0)+timesI(fspin(1));
  hspin(1)=fspin(1)-timesI(fspin(0));
}
template<class vtype,IfSpinor<iVector<vtype,Nhs> > = 0> accelerator_inline void pauliProjZp (iVector<vtype,Nhs> &hspin,const iVector<vtype,Nhs> &fspin)
{
  hspin(0)=fspin(0)+fspin(0);
  hspin(1)=Zero();
}
template<class vtype,IfSpinor<iVector<vtype,Nhs> > = 0> accelerator_inline void pauliProjZm (iVector<vtype,Nhs> &hspin,const iVector<vtype,Nhs> &fspin)
{
  hspin(0)=Zero();
  hspin(1)=fspin(1)+fspin(1);
}
template<class vtype,IfSpinor<iVector<vtype,Nhs> > = 0> accelerator_inline void pauliAssign(iVector<vtype,Nhs> &fspin,const iVector<vtype,Nhs> &hspin)
{
  fspin = hspin;
}
template<class vtype,IfSpinor<iVector<vtype,Nhs> > = 0> accelerator_inline void pauliAdd   (iVector<vtype,Nhs> &fspin,const iVector<vtype,Nhs> &hspin)
{
  fspin = fspin + hspin;
}

template<class vtype> 
accelerator_inline auto operator*(const Pauli &G, const iVector<vtype, Nhs> &arg)
->typename std::enable_if<matchGridTensorIndex<iVector<vtype, Nhs>, PauliIndex>::value, iVector<vtype, Nhs>>::type
{
  iVector<vtype, Nhs> ret;

  switch (G.g) 
  {
  case Pauli::Algebra::SigmaX:
    multPauliSigmaX(ret, arg); break;
  case Pauli::Algebra::MinusSigmaX:
    multPauliMinusSigmaX(ret, arg); break;
  case Pauli::Algebra::SigmaY:
    multPauliSigmaY(ret, arg); break;
  case Pauli::Algebra::MinusSigmaY:
    multPauliMinusSigmaY(ret, arg); break;
  case Pauli::Algebra::SigmaZ:
    multPauliSigmaZ(ret, arg); break;
  case Pauli::Algebra::MinusSigmaZ:
    multPauliMinusSigmaZ(ret, arg); break;
  case Pauli::Algebra::Identity:
    multPauliIdentity(ret, arg); break;
  case Pauli::Algebra::MinusIdentity:
    multPauliMinusIdentity(ret, arg); break;
  default: assert(0);
  }
 
  return ret;
}

template<class vtype> 
accelerator_inline auto operator*(const Pauli &G, const iMatrix<vtype, Nhs> &arg)
->typename std::enable_if<matchGridTensorIndex<iMatrix<vtype, Nhs>, PauliIndex>::value, iMatrix<vtype, Nhs>>::type
{
  iMatrix<vtype, Nhs> ret;

  switch (G.g) 
  {
  case Pauli::Algebra::SigmaX:
    lmultPauliSigmaX(ret, arg); break;
  case Pauli::Algebra::MinusSigmaX:
    lmultPauliMinusSigmaX(ret, arg); break;
  case Pauli::Algebra::SigmaY:
    lmultPauliSigmaY(ret, arg); break;
  case Pauli::Algebra::MinusSigmaY:
    lmultPauliMinusSigmaY(ret, arg); break;
  case Pauli::Algebra::SigmaZ:
    lmultPauliSigmaZ(ret, arg); break;
  case Pauli::Algebra::MinusSigmaZ:
    lmultPauliMinusSigmaZ(ret, arg); break;
  case Pauli::Algebra::Identity:
    lmultPauliIdentity(ret, arg); break;
  case Pauli::Algebra::MinusIdentity:
    lmultPauliMinusIdentity(ret, arg); break;
  default: assert(0);
  }
  
  return ret;
}

template<class vtype> 
accelerator_inline auto operator*(const iMatrix<vtype, Nhs> &arg, const Pauli &G)
->typename std::enable_if<matchGridTensorIndex<iMatrix<vtype, Nhs>, PauliIndex>::value, iMatrix<vtype, Nhs>>::type
{
  iMatrix<vtype, Nhs> ret;

  switch (G.g) 
  {
  case Pauli::Algebra::SigmaX:
    rmultPauliSigmaX(ret, arg); break;
  case Pauli::Algebra::MinusSigmaX:
    rmultPauliMinusSigmaX(ret, arg); break;
  case Pauli::Algebra::SigmaY:
    rmultPauliSigmaY(ret, arg); break;
  case Pauli::Algebra::MinusSigmaY:
    rmultPauliMinusSigmaY(ret, arg); break;
  case Pauli::Algebra::SigmaZ:
    rmultPauliSigmaZ(ret, arg); break;
  case Pauli::Algebra::MinusSigmaZ:
    rmultPauliMinusSigmaZ(ret, arg); break;
  case Pauli::Algebra::Identity:
    rmultPauliIdentity(ret, arg); break;
  case Pauli::Algebra::MinusIdentity:
    rmultPauliMinusIdentity(ret, arg); break;
  default: assert(0);
  }

  return ret;
}


NAMESPACE_END(Grid);

#endif // GRID_QCD_GAMMA_H
