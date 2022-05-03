#ifndef GRID_QCD_GPARITY_FLAVOUR_H
#define GRID_QCD_GPARITY_FLAVOUR_H

//Support for flavour-matrix operations acting on the G-parity flavour index

#include <array>

NAMESPACE_BEGIN(Grid);

class GparityFlavour {
  public:
    GRID_SERIALIZABLE_ENUM(Algebra, undef,
                           SigmaX, 0,
			   MinusSigmaX, 1,
                           SigmaY, 2,
			   MinusSigmaY, 3,
                           SigmaZ, 4,
			   MinusSigmaZ, 5,
			   Identity, 6,
			   MinusIdentity, 7,
			   ProjPlus, 8,
			   MinusProjPlus, 9,
			   ProjMinus, 10,
			   MinusProjMinus, 11
			   );
    static constexpr unsigned int nSigma = 12;
    static const std::array<const char *, nSigma>                name;
    static const std::array<const GparityFlavour, 3>             sigma_mu;
    static const std::array<const GparityFlavour, 6>            sigma_all;
    Algebra                                                      g;
  public:
  accelerator GparityFlavour(Algebra initg): g(initg) {}  
};



// 0 1  x   vector
// 1 0
template<class vtype>
accelerator_inline void multFlavourSigmaX(iVector<vtype, Ngp> &ret, const iVector<vtype, Ngp> &rhs)
{
  ret(0) = rhs(1);
  ret(1) = rhs(0);
};
template<class vtype>
accelerator_inline void lmultFlavourSigmaX(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = rhs(1,0);
  ret(0,1) = rhs(1,1);
  ret(1,0) = rhs(0,0);
  ret(1,1) = rhs(0,1);
};
template<class vtype>
accelerator_inline void rmultFlavourSigmaX(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = rhs(0,1);
  ret(0,1) = rhs(0,0);
  ret(1,0) = rhs(1,1);
  ret(1,1) = rhs(1,0);
};


template<class vtype>
accelerator_inline void multFlavourMinusSigmaX(iVector<vtype, Ngp> &ret, const iVector<vtype, Ngp> &rhs)
{
  ret(0) = -rhs(1);
  ret(1) = -rhs(0);
};
template<class vtype>
accelerator_inline void lmultFlavourMinusSigmaX(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = -rhs(1,0);
  ret(0,1) = -rhs(1,1);
  ret(1,0) = -rhs(0,0);
  ret(1,1) = -rhs(0,1);
};
template<class vtype>
accelerator_inline void rmultFlavourMinusSigmaX(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = -rhs(0,1);
  ret(0,1) = -rhs(0,0);
  ret(1,0) = -rhs(1,1);
  ret(1,1) = -rhs(1,0);
};





// 0 -i  x   vector
// i 0
template<class vtype>
accelerator_inline void multFlavourSigmaY(iVector<vtype, Ngp> &ret, const iVector<vtype, Ngp> &rhs)
{
  ret(0) = timesMinusI(rhs(1));
  ret(1) = timesI(rhs(0));
};
template<class vtype>
accelerator_inline void lmultFlavourSigmaY(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = timesMinusI(rhs(1,0));
  ret(0,1) = timesMinusI(rhs(1,1));
  ret(1,0) = timesI(rhs(0,0));
  ret(1,1) = timesI(rhs(0,1));
};
template<class vtype>
accelerator_inline void rmultFlavourSigmaY(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = timesI(rhs(0,1));
  ret(0,1) = timesMinusI(rhs(0,0));
  ret(1,0) = timesI(rhs(1,1));
  ret(1,1) = timesMinusI(rhs(1,0));
};

template<class vtype>
accelerator_inline void multFlavourMinusSigmaY(iVector<vtype, Ngp> &ret, const iVector<vtype, Ngp> &rhs)
{
  ret(0) = timesI(rhs(1));
  ret(1) = timesMinusI(rhs(0));
};
template<class vtype>
accelerator_inline void lmultFlavourMinusSigmaY(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = timesI(rhs(1,0));
  ret(0,1) = timesI(rhs(1,1));
  ret(1,0) = timesMinusI(rhs(0,0));
  ret(1,1) = timesMinusI(rhs(0,1));
};
template<class vtype>
accelerator_inline void rmultFlavourMinusSigmaY(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = timesMinusI(rhs(0,1));
  ret(0,1) = timesI(rhs(0,0));
  ret(1,0) = timesMinusI(rhs(1,1));
  ret(1,1) = timesI(rhs(1,0));
};





// 1 0  x   vector
// 0 -1
template<class vtype>
accelerator_inline void multFlavourSigmaZ(iVector<vtype, Ngp> &ret, const iVector<vtype, Ngp> &rhs)
{
  ret(0) = rhs(0);
  ret(1) = -rhs(1);
};
template<class vtype>
accelerator_inline void lmultFlavourSigmaZ(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = rhs(0,0);
  ret(0,1) = rhs(0,1);
  ret(1,0) = -rhs(1,0);
  ret(1,1) = -rhs(1,1);
};
template<class vtype>
accelerator_inline void rmultFlavourSigmaZ(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = rhs(0,0);
  ret(0,1) = -rhs(0,1);
  ret(1,0) = rhs(1,0);
  ret(1,1) = -rhs(1,1);
};


template<class vtype>
accelerator_inline void multFlavourMinusSigmaZ(iVector<vtype, Ngp> &ret, const iVector<vtype, Ngp> &rhs)
{
  ret(0) = -rhs(0);
  ret(1) = rhs(1);
};
template<class vtype>
accelerator_inline void lmultFlavourMinusSigmaZ(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = -rhs(0,0);
  ret(0,1) = -rhs(0,1);
  ret(1,0) = rhs(1,0);
  ret(1,1) = rhs(1,1);
};
template<class vtype>
accelerator_inline void rmultFlavourMinusSigmaZ(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = -rhs(0,0);
  ret(0,1) = rhs(0,1);
  ret(1,0) = -rhs(1,0);
  ret(1,1) = rhs(1,1);
};






template<class vtype>
accelerator_inline void multFlavourIdentity(iVector<vtype, Ngp> &ret, const iVector<vtype, Ngp> &rhs)
{
  ret(0) = rhs(0);
  ret(1) = rhs(1);
};
template<class vtype>
accelerator_inline void lmultFlavourIdentity(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = rhs(0,0);
  ret(0,1) = rhs(0,1);
  ret(1,0) = rhs(1,0);
  ret(1,1) = rhs(1,1);
};
template<class vtype>
accelerator_inline void rmultFlavourIdentity(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = rhs(0,0);
  ret(0,1) = rhs(0,1);
  ret(1,0) = rhs(1,0);
  ret(1,1) = rhs(1,1);
};

template<class vtype>
accelerator_inline void multFlavourMinusIdentity(iVector<vtype, Ngp> &ret, const iVector<vtype, Ngp> &rhs)
{
  ret(0) = -rhs(0);
  ret(1) = -rhs(1);
};
template<class vtype>
accelerator_inline void lmultFlavourMinusIdentity(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = -rhs(0,0);
  ret(0,1) = -rhs(0,1);
  ret(1,0) = -rhs(1,0);
  ret(1,1) = -rhs(1,1);
};
template<class vtype>
accelerator_inline void rmultFlavourMinusIdentity(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = -rhs(0,0);
  ret(0,1) = -rhs(0,1);
  ret(1,0) = -rhs(1,0);
  ret(1,1) = -rhs(1,1);
};





//G-parity flavour projection 1/2(1+\sigma_2)
//1 -i
//i  1
template<class vtype>
accelerator_inline void multFlavourProjPlus(iVector<vtype, Ngp> &ret, const iVector<vtype, Ngp> &rhs)
{
  ret(0) = 0.5*rhs(0) + 0.5*timesMinusI(rhs(1));
  ret(1) = 0.5*timesI(rhs(0)) + 0.5*rhs(1);
};
template<class vtype>
accelerator_inline void lmultFlavourProjPlus(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = 0.5*rhs(0,0) + 0.5*timesMinusI(rhs(1,0));
  ret(0,1) = 0.5*rhs(0,1) + 0.5*timesMinusI(rhs(1,1));
  ret(1,0) = 0.5*timesI(rhs(0,0)) + 0.5*rhs(1,0);
  ret(1,1) = 0.5*timesI(rhs(0,1)) + 0.5*rhs(1,1);
};
template<class vtype>
accelerator_inline void rmultFlavourProjPlus(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = 0.5*rhs(0,0) + 0.5*timesI(rhs(0,1));
  ret(0,1) = 0.5*timesMinusI(rhs(0,0)) + 0.5*rhs(0,1);
  ret(1,0) = 0.5*rhs(1,0) + 0.5*timesI(rhs(1,1));
  ret(1,1) = 0.5*timesMinusI(rhs(1,0)) + 0.5*rhs(1,1);
};


template<class vtype>
accelerator_inline void multFlavourMinusProjPlus(iVector<vtype, Ngp> &ret, const iVector<vtype, Ngp> &rhs)
{
  ret(0) = -0.5*rhs(0) + 0.5*timesI(rhs(1));
  ret(1) = 0.5*timesMinusI(rhs(0)) - 0.5*rhs(1);
};
template<class vtype>
accelerator_inline void lmultFlavourMinusProjPlus(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = -0.5*rhs(0,0) + 0.5*timesI(rhs(1,0));
  ret(0,1) = -0.5*rhs(0,1) + 0.5*timesI(rhs(1,1));
  ret(1,0) = 0.5*timesMinusI(rhs(0,0)) - 0.5*rhs(1,0);
  ret(1,1) = 0.5*timesMinusI(rhs(0,1)) - 0.5*rhs(1,1);
};
template<class vtype>
accelerator_inline void rmultFlavourMinusProjPlus(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = -0.5*rhs(0,0) + 0.5*timesMinusI(rhs(0,1));
  ret(0,1) = 0.5*timesI(rhs(0,0)) - 0.5*rhs(0,1);
  ret(1,0) = -0.5*rhs(1,0) + 0.5*timesMinusI(rhs(1,1));
  ret(1,1) = 0.5*timesI(rhs(1,0)) - 0.5*rhs(1,1);
};





//G-parity flavour projection 1/2(1-\sigma_2)
//1 i
//-i  1
template<class vtype>
accelerator_inline void multFlavourProjMinus(iVector<vtype, Ngp> &ret, const iVector<vtype, Ngp> &rhs)
{
  ret(0) = 0.5*rhs(0) + 0.5*timesI(rhs(1));
  ret(1) = 0.5*timesMinusI(rhs(0)) + 0.5*rhs(1);
};
template<class vtype>
accelerator_inline void lmultFlavourProjMinus(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = 0.5*rhs(0,0) + 0.5*timesI(rhs(1,0));
  ret(0,1) = 0.5*rhs(0,1) + 0.5*timesI(rhs(1,1));
  ret(1,0) = 0.5*timesMinusI(rhs(0,0)) + 0.5*rhs(1,0);
  ret(1,1) = 0.5*timesMinusI(rhs(0,1)) + 0.5*rhs(1,1);
};
template<class vtype>
accelerator_inline void rmultFlavourProjMinus(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = 0.5*rhs(0,0) + 0.5*timesMinusI(rhs(0,1));
  ret(0,1) = 0.5*timesI(rhs(0,0)) + 0.5*rhs(0,1);
  ret(1,0) = 0.5*rhs(1,0) + 0.5*timesMinusI(rhs(1,1));
  ret(1,1) = 0.5*timesI(rhs(1,0)) + 0.5*rhs(1,1);
};


template<class vtype>
accelerator_inline void multFlavourMinusProjMinus(iVector<vtype, Ngp> &ret, const iVector<vtype, Ngp> &rhs)
{
  ret(0) = -0.5*rhs(0) + 0.5*timesMinusI(rhs(1));
  ret(1) = 0.5*timesI(rhs(0)) - 0.5*rhs(1);
};
template<class vtype>
accelerator_inline void lmultFlavourMinusProjMinus(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = -0.5*rhs(0,0) + 0.5*timesMinusI(rhs(1,0));
  ret(0,1) = -0.5*rhs(0,1) + 0.5*timesMinusI(rhs(1,1));
  ret(1,0) = 0.5*timesI(rhs(0,0)) - 0.5*rhs(1,0);
  ret(1,1) = 0.5*timesI(rhs(0,1)) - 0.5*rhs(1,1);
};
template<class vtype>
accelerator_inline void rmultFlavourMinusProjMinus(iMatrix<vtype, Ngp> &ret, const iMatrix<vtype, Ngp> &rhs)
{
  ret(0,0) = -0.5*rhs(0,0) + 0.5*timesI(rhs(0,1));
  ret(0,1) = 0.5*timesMinusI(rhs(0,0)) - 0.5*rhs(0,1);
  ret(1,0) = -0.5*rhs(1,0) + 0.5*timesI(rhs(1,1));
  ret(1,1) = 0.5*timesMinusI(rhs(1,0)) - 0.5*rhs(1,1);
};










template<class vtype> 
accelerator_inline auto operator*(const GparityFlavour &G, const iVector<vtype, Ngp> &arg)
->typename std::enable_if<matchGridTensorIndex<iVector<vtype, Ngp>, GparityFlavourTensorIndex>::value, iVector<vtype, Ngp>>::type
{
  iVector<vtype, Ngp> ret;

  switch (G.g) 
  {
  case GparityFlavour::Algebra::SigmaX:
    multFlavourSigmaX(ret, arg); break;
  case GparityFlavour::Algebra::MinusSigmaX:
    multFlavourMinusSigmaX(ret, arg); break;
  case GparityFlavour::Algebra::SigmaY:
    multFlavourSigmaY(ret, arg); break;
  case GparityFlavour::Algebra::MinusSigmaY:
    multFlavourMinusSigmaY(ret, arg); break;
  case GparityFlavour::Algebra::SigmaZ:
    multFlavourSigmaZ(ret, arg); break;
  case GparityFlavour::Algebra::MinusSigmaZ:
    multFlavourMinusSigmaZ(ret, arg); break;
  case GparityFlavour::Algebra::Identity:
    multFlavourIdentity(ret, arg); break;
  case GparityFlavour::Algebra::MinusIdentity:
    multFlavourMinusIdentity(ret, arg); break;
  case GparityFlavour::Algebra::ProjPlus:
    multFlavourProjPlus(ret, arg); break;
  case GparityFlavour::Algebra::MinusProjPlus:
    multFlavourMinusProjPlus(ret, arg); break;
  case GparityFlavour::Algebra::ProjMinus:
    multFlavourProjMinus(ret, arg); break;
  case GparityFlavour::Algebra::MinusProjMinus:
    multFlavourMinusProjMinus(ret, arg); break;
  default: assert(0);
  }
 
  return ret;
}

template<class vtype> 
accelerator_inline auto operator*(const GparityFlavour &G, const iMatrix<vtype, Ngp> &arg)
->typename std::enable_if<matchGridTensorIndex<iMatrix<vtype, Ngp>, GparityFlavourTensorIndex>::value, iMatrix<vtype, Ngp>>::type
{
  iMatrix<vtype, Ngp> ret;

  switch (G.g) 
  {
  case GparityFlavour::Algebra::SigmaX:
    lmultFlavourSigmaX(ret, arg); break;
  case GparityFlavour::Algebra::MinusSigmaX:
    lmultFlavourMinusSigmaX(ret, arg); break;
  case GparityFlavour::Algebra::SigmaY:
    lmultFlavourSigmaY(ret, arg); break;
  case GparityFlavour::Algebra::MinusSigmaY:
    lmultFlavourMinusSigmaY(ret, arg); break;
  case GparityFlavour::Algebra::SigmaZ:
    lmultFlavourSigmaZ(ret, arg); break;
  case GparityFlavour::Algebra::MinusSigmaZ:
    lmultFlavourMinusSigmaZ(ret, arg); break;
  case GparityFlavour::Algebra::Identity:
    lmultFlavourIdentity(ret, arg); break;
  case GparityFlavour::Algebra::MinusIdentity:
    lmultFlavourMinusIdentity(ret, arg); break;
  case GparityFlavour::Algebra::ProjPlus:
    lmultFlavourProjPlus(ret, arg); break;
  case GparityFlavour::Algebra::MinusProjPlus:
    lmultFlavourMinusProjPlus(ret, arg); break;
  case GparityFlavour::Algebra::ProjMinus:
    lmultFlavourProjMinus(ret, arg); break;
  case GparityFlavour::Algebra::MinusProjMinus:
    lmultFlavourMinusProjMinus(ret, arg); break;  
  default: assert(0);
  }
  
  return ret;
}

template<class vtype> 
accelerator_inline auto operator*(const iMatrix<vtype, Ngp> &arg, const GparityFlavour &G)
->typename std::enable_if<matchGridTensorIndex<iMatrix<vtype, Ngp>, GparityFlavourTensorIndex>::value, iMatrix<vtype, Ngp>>::type
{
  iMatrix<vtype, Ngp> ret;

  switch (G.g) 
  {
  case GparityFlavour::Algebra::SigmaX:
    rmultFlavourSigmaX(ret, arg); break;
  case GparityFlavour::Algebra::MinusSigmaX:
    rmultFlavourMinusSigmaX(ret, arg); break;
  case GparityFlavour::Algebra::SigmaY:
    rmultFlavourSigmaY(ret, arg); break;
  case GparityFlavour::Algebra::MinusSigmaY:
    rmultFlavourMinusSigmaY(ret, arg); break;
  case GparityFlavour::Algebra::SigmaZ:
    rmultFlavourSigmaZ(ret, arg); break;
  case GparityFlavour::Algebra::MinusSigmaZ:
    rmultFlavourMinusSigmaZ(ret, arg); break;
  case GparityFlavour::Algebra::Identity:
    rmultFlavourIdentity(ret, arg); break;
  case GparityFlavour::Algebra::MinusIdentity:
    rmultFlavourMinusIdentity(ret, arg); break;
  case GparityFlavour::Algebra::ProjPlus:
    rmultFlavourProjPlus(ret, arg); break;
  case GparityFlavour::Algebra::MinusProjPlus:
    rmultFlavourMinusProjPlus(ret, arg); break;
  case GparityFlavour::Algebra::ProjMinus:
    rmultFlavourProjMinus(ret, arg); break;
  case GparityFlavour::Algebra::MinusProjMinus:
    rmultFlavourMinusProjMinus(ret, arg); break;
  default: assert(0);
  }

  return ret;
}

NAMESPACE_END(Grid);

#endif // include guard
