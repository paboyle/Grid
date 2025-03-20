// This file is #included into the body of the class template definition of
// GaugeGroup. So, image there to be
//
// template <int ncolour, class group_name>
// class GaugeGroup {
//
// around it.
//
// Please note that the unconventional file extension makes sure that it
// doesn't get found by the scripts/filelist during bootstrapping.

private:
template <ONLY_IF_Sp>
static int su2subgroups(GroupName::Sp) { return (ncolour/2 * (ncolour/2 - 1)) / 2; }

// Sp(2N) has N(2N+1) = 2N^2+N generators
//
// normalise the generators such that
// Trace ( Ta Tb) = 1/2 delta_ab
//
// N generators in the cartan, 2N^2 off
// off diagonal:
//     there are 6 types named a,b,c,d and w,z
//     abcd are N(N-1)/2 each while wz are N each

template <class cplx, ONLY_IF_Sp>
static void generator(int lieIndex, iGroupMatrix<cplx> &ta, GroupName::Sp) {
  // map lie index into type of generators: diagonal, abcd type, wz type

  const int nsp = ncolour/2;
  int diagIndex;
  int aIndex, bIndex, cIndex, dIndex;
  int wIndex, zIndex;  // a,b,c,d are N(N-1)/2 and w,z are N
  const int mod = nsp * (nsp - 1) * 0.5;
  const int offdiag =
      2 * nsp * nsp;  // number of generators not in the cartan subalgebra
  const int wmod = 4 * mod;
  const int zmod = wmod + nsp;
  if (lieIndex >= offdiag) {
    diagIndex = lieIndex - offdiag;  // 0, ... ,N-1
    // std::cout << GridLogMessage << "diag type " << std::endl;
    generatorDiagtype(diagIndex, ta);
    return;
  }
  if ((lieIndex >= wmod) && (lieIndex < zmod)) {
    // std::cout << GridLogMessage << "w type " << std::endl;
    wIndex = lieIndex - wmod;  // 0, ... ,N-1
    generatorWtype(wIndex, ta);
    return;
  }
  if ((lieIndex >= zmod) && (lieIndex < offdiag)) {
    // std::cout << GridLogMessage << "z type " << std::endl;
    // std::cout << GridLogMessage << "lie index " << lieIndex << std::endl;
    // std::cout << GridLogMessage << "z mod " << zmod << std::endl;
    zIndex = lieIndex - zmod;  // 0, ... ,N-1
    generatorZtype(zIndex, ta);
    return;
  }
  if (lieIndex < mod) {  // atype 0, ... , N(N-1)/2=mod
    // std::cout << GridLogMessage << "a type " << std::endl;
    aIndex = lieIndex;
    // std::cout << GridLogMessage << "a indx " << aIndex << std::endl;
    generatorAtype(aIndex, ta);
    return;
  }
  if ((lieIndex >= mod) && lieIndex < 2 * mod) {  // btype mod, ... , 2mod-1
    // std::cout << GridLogMessage << "b type " << std::endl;
    bIndex = lieIndex - mod;
    generatorBtype(bIndex, ta);
    return;
  }
  if ((lieIndex >= 2 * mod) &&
      lieIndex < 3 * mod) {  // ctype 2mod, ... , 3mod-1
    // std::cout << GridLogMessage << "c type " << std::endl;
    cIndex = lieIndex - 2 * mod;
    generatorCtype(cIndex, ta);
    return;
  }
  if ((lieIndex >= 3 * mod) &&
      lieIndex < wmod) {  // ctype 3mod, ... , 4mod-1 = wmod-1
    // std::cout << GridLogMessage << "d type " << std::endl;
    dIndex = lieIndex - 3 * mod;
    generatorDtype(dIndex, ta);
    return;
  }

}  // end of generator

template <class cplx, ONLY_IF_Sp>
static void generatorDiagtype(int diagIndex, iGroupMatrix<cplx> &ta) {
  // ta(i,i) = - ta(i+N,i+N) = 1/2 for each i index of the cartan subalgebra

  const int nsp=ncolour/2;
  ta = Zero();
  RealD nrm = 1.0 / 2;

  ta()()(diagIndex, diagIndex) = nrm;
  ta()()(diagIndex + nsp, diagIndex + nsp) = -nrm;
}

template <class cplx, ONLY_IF_Sp>
static void generatorAtype(int aIndex, iGroupMatrix<cplx> &ta) {
  // ta(i,j) = ta(j,i) = -ta(i+N,j+N) = -ta(j+N,i+N) = 1 / 2 sqrt(2)
  // with i<j and i=0,...,N-2
  // follows that j=i+1, ... , N
  int i1, i2;
  const int nsp=ncolour/2;
  ta = Zero();
  RealD nrm = 1 / (2 * std::sqrt(2));

  su2SubGroupIndex(i1, i2, aIndex);
  ta()()(i1, i2) = 1;
  ta()()(i2, i1) = 1;
  ta()()(i1 + nsp, i2 + nsp) = -1;
  ta()()(i2 + nsp, i1 + nsp) = -1;

  ta = ta * nrm;
}

template <class cplx, ONLY_IF_Sp>
static void generatorBtype(int bIndex, iGroupMatrix<cplx> &ta) {
  // ta(i,j) = -ta(j,i) = ta(i+N,j+N) = -ta(j+N,i+N) = i / 1/ 2 sqrt(2)
  // with i<j and i=0,...,N-2
  // follows that j=i+1, ... , N-1

  const int nsp=ncolour/2;
  int i1, i2;
  ta = Zero();
  cplx i(0.0, 1.0);
  RealD nrm = 1 / (2 * std::sqrt(2));
  su2SubGroupIndex(i1, i2, bIndex);

  ta()()(i1, i2) = i;
  ta()()(i2, i1) = -i;
  ta()()(i1 + nsp, i2 + nsp) = i;
  ta()()(i2 + nsp, i1 + nsp) = -i;

  ta = ta * nrm;
}

template <class cplx, ONLY_IF_Sp>
static void generatorCtype(int cIndex, iGroupMatrix<cplx> &ta) {
  // ta(i,j+N) = ta(j,i+N) = ta(i+N,j) = ta(j+N,i) = 1 / 2 sqrt(2)

  const int nsp=ncolour/2;
  int i1, i2;
  ta = Zero();
  RealD nrm = 1 / (2 * std::sqrt(2));
  su2SubGroupIndex(i1, i2, cIndex);

  ta()()(i1, i2 + nsp) = 1;
  ta()()(i2, i1 + nsp) = 1;
  ta()()(i1 + nsp, i2) = 1;
  ta()()(i2 + nsp, i1) = 1;

  ta = ta * nrm;
}

template <class cplx, ONLY_IF_Sp>
static void generatorDtype(int dIndex, iGroupMatrix<cplx> &ta) {
  // ta(i,j+N) = ta(j,i+N) = -ta(i+N,j) = -ta(j+N,i) = i /  2 sqrt(2)

  const int nsp=ncolour/2;
  int i1, i2;
  ta = Zero();
  cplx i(0.0, 1.0);
  RealD nrm = 1 / (2 * std::sqrt(2));
  su2SubGroupIndex(i1, i2, dIndex);

  ta()()(i1, i2 + nsp) = i;
  ta()()(i2, i1 + nsp) = i;
  ta()()(i1 + nsp, i2) = -i;
  ta()()(i2 + nsp, i1) = -i;

  ta = ta * nrm;
}

template <class cplx, ONLY_IF_Sp>
static void generatorWtype(int wIndex, iGroupMatrix<cplx> &ta) {
  // ta(i,i+N) =  ta(i+N,i) = 1/2

  const int nsp=ncolour/2;
  ta = Zero();
  RealD nrm = 1.0 / 2;  // check

  ta()()(wIndex, wIndex + nsp) = 1;
  ta()()(wIndex + nsp, wIndex) = 1;

  ta = ta * nrm;
}

template <class cplx, ONLY_IF_Sp>
static void generatorZtype(int zIndex, iGroupMatrix<cplx> &ta) {
  // ta(i,i+N) = - ta(i+N,i) = i/2

  const int nsp=ncolour/2;
  ta = Zero();
  RealD nrm = 1.0 / 2;  // check
  cplx i(0.0, 1.0);
  ta()()(zIndex, zIndex + nsp) = i;
  ta()()(zIndex + nsp, zIndex) = -i;

  ta = ta * nrm;
}

////////////////////////////////////////////////////////////////////////
// Map a su2 subgroup number to the pair of rows that are non zero
////////////////////////////////////////////////////////////////////////
template <ONLY_IF_Sp>
static accelerator_inline void su2SubGroupIndex(int &i1, int &i2, int su2_index, GroupName::Sp) {
  const int nsp=ncolour/2;
  assert((su2_index >= 0) && (su2_index < (nsp * (nsp - 1)) / 2));

  int spare = su2_index;
  for (i1 = 0; spare >= (nsp - 1 - i1); i1++) {
    spare = spare - (nsp - 1 - i1);  // remove the Nc-1-i1 terms
  }
  i2 = i1 + 1 + spare;
}

static void testGenerators(GroupName::Sp) {
  Matrix ta;
  Matrix tb;
  std::cout << GridLogMessage
            << "Fundamental - Checking trace ta tb is 0.5 delta_ab "
            << std::endl;
  for (int a = 0; a < AlgebraDimension; a++) {
    for (int b = 0; b < AlgebraDimension; b++) {
      generator(a, ta);
      generator(b, tb);
      Complex tr = TensorRemove(trace(ta * tb));
      std::cout << GridLogMessage << "(" << a << "," << b << ") =  " << tr
                << std::endl;
      if (a == b) assert(abs(tr - Complex(0.5)) < 1.0e-6);
      if (a != b) assert(abs(tr) < 1.0e-6);
    }
  }
  std::cout << GridLogMessage << std::endl;
  std::cout << GridLogMessage << "Fundamental - Checking if hermitian"
            << std::endl;
  for (int a = 0; a < AlgebraDimension; a++) {
    generator(a, ta);
    std::cout << GridLogMessage << a << std::endl;
    assert(norm2(ta - adj(ta)) < 1.0e-6);
  }
  std::cout << GridLogMessage << std::endl;
  std::cout << GridLogMessage << "Fundamental - Checking if traceless"
            << std::endl;
  for (int a = 0; a < AlgebraDimension; a++) {
    generator(a, ta);
    Complex tr = TensorRemove(trace(ta));
    std::cout << GridLogMessage << a << std::endl;
    assert(abs(tr) < 1.0e-6);
  }
}

template <int N>
static Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > >
ProjectOnGeneralGroup(const Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > &Umu, GroupName::Sp) {
  return ProjectOnSpGroup(Umu);
}

template <class vtype>
accelerator_inline static iScalar<vtype> ProjectOnGeneralGroup(const iScalar<vtype> &r, GroupName::Sp) {
  return ProjectOnSpGroup(r);
}

template <class vtype, int N>
accelerator_inline static iVector<vtype,N> ProjectOnGeneralGroup(const iVector<vtype,N> &r, GroupName::Sp) {
  return ProjectOnSpGroup(r);
}

template <class vtype,int N, typename std::enable_if< GridTypeMapper<vtype>::TensorLevel == 0 >::type * =nullptr>
accelerator_inline static iMatrix<vtype,N> ProjectOnGeneralGroup(const iMatrix<vtype,N> &arg, GroupName::Sp) {
  return ProjectOnSpGroup(arg);
}

template <typename LatticeMatrixType>   
static void taProj(const LatticeMatrixType &in, LatticeMatrixType &out, GroupName::Sp) {
  out = SpTa(in);
}

public:

template <ONLY_IF_Sp>
static void Omega(LatticeColourMatrixD &in) {
  const int nsp=ncolour/2;
  LatticeColourMatrixD OmegaLatt(in.Grid());
  LatticeColourMatrixD identity(in.Grid());
  ColourMatrix Omega;

  OmegaLatt = Zero();
  Omega = Zero();
  identity = 1.;

  for (int i = 0; i < nsp; i++) {
    Omega()()(i, nsp + i) = 1.;
    Omega()()(nsp + i, i) = -1;
  }
  OmegaLatt = OmegaLatt + (identity * Omega);
  in = OmegaLatt;
}

template <ONLY_IF_Sp, class vtype, int N>
static void Omega(iScalar<iScalar<iMatrix<vtype, N> > > &in) {
  const int nsp=ncolour/2;
    
  iScalar<iScalar<iMatrix<vtype, N> > > Omega;
  Omega = Zero();

  for (int i = 0; i < nsp; i++) {
    Omega()()(i, nsp + i) = 1.;
    Omega()()(nsp + i, i) = -1;
  }
    
  in = Omega;
}
