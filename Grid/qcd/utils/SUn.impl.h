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

template <ONLY_IF_SU>
static int su2subgroups(GroupName::SU) { return (ncolour * (ncolour - 1)) / 2; }
////////////////////////////////////////////////////////////////////////
// There are N^2-1 generators for SU(N).
//
// We take a traceless hermitian generator basis as follows
//
// * Normalisation: trace ta tb = 1/2 delta_ab = T_F delta_ab
//   T_F = 1/2  for SU(N) groups
//
// * Off diagonal
//    - pairs of rows i1,i2 behaving like pauli matrices signma_x, sigma_y
//
//    - there are (Nc-1-i1) slots for i2 on each row [ x  0  x ]
//      direct count off each row
//
//    - Sum of all pairs is Nc(Nc-1)/2: proof arithmetic series
//
//      (Nc-1) + (Nc-2)+...  1      ==> Nc*(Nc-1)/2
//      1+ 2+          +   + Nc-1
//
//    - There are 2 x Nc (Nc-1)/ 2 of these = Nc^2 - Nc
//
//    - We enumerate the row-col pairs.
//    - for each row col pair there is a (sigma_x) and a (sigma_y) like
//    generator
//
//
//   t^a_ij = { in 0.. Nc(Nc-1)/2 -1} =>  1/2(delta_{i,i1} delta_{j,i2} +
//   delta_{i,i1} delta_{j,i2})
//   t^a_ij = { in Nc(Nc-1)/2 ... Nc(Nc-1) - 1} =>  i/2( delta_{i,i1}
//   delta_{j,i2} - i delta_{i,i1} delta_{j,i2})
//
// * Diagonal; must be traceless and normalised
//   - Sequence is
//   N  (1,-1,0,0...)
//   N  (1, 1,-2,0...)
//   N  (1, 1, 1,-3,0...)
//   N  (1, 1, 1, 1,-4,0...)
//
//   where 1/2 = N^2 (1+.. m^2)etc.... for the m-th diagonal generator
//   NB this gives the famous SU3 result for su2 index 8
//
//   N= sqrt(1/2 . 1/6 ) = 1/2 . 1/sqrt(3)
//
//   ( 1      )
//   (    1   ) / sqrt(3) /2  = 1/2 lambda_8
//   (      -2)
//
////////////////////////////////////////////////////////////////////////
template <class cplx, ONLY_IF_SU>
static void generator(int lieIndex, iGroupMatrix<cplx> &ta, GroupName::SU) {
  // map lie index to which type of generator
  int diagIndex;
  int su2Index;
  int sigxy;
  int NNm1 = ncolour * (ncolour - 1);
  if (lieIndex >= NNm1) {
    diagIndex = lieIndex - NNm1;
    generatorDiagonal(diagIndex, ta);
    return;
  }
  sigxy = lieIndex & 0x1;  // even or odd
  su2Index = lieIndex >> 1;
  if (sigxy)
    generatorSigmaY(su2Index, ta);
  else
    generatorSigmaX(su2Index, ta);
}

template <class cplx, ONLY_IF_SU>
static void generatorSigmaY(int su2Index, iGroupMatrix<cplx> &ta) {
  ta = Zero();
  int i1, i2;
  su2SubGroupIndex(i1, i2, su2Index);
  ta()()(i1, i2) = 1.0;
  ta()()(i2, i1) = 1.0;
  ta = ta * 0.5;
}

template <class cplx, ONLY_IF_SU>
static void generatorSigmaX(int su2Index, iGroupMatrix<cplx> &ta) {
  ta = Zero();
  cplx i(0.0, 1.0);
  int i1, i2;
  su2SubGroupIndex(i1, i2, su2Index);
  ta()()(i1, i2) = i;
  ta()()(i2, i1) = -i;
  ta = ta * 0.5;
}

template <class cplx, ONLY_IF_SU>
static void generatorDiagonal(int diagIndex, iGroupMatrix<cplx> &ta) {
  // diag ({1, 1, ..., 1}(k-times), -k, 0, 0, ...)
  ta = Zero();
  int k = diagIndex + 1;                  // diagIndex starts from 0
  for (int i = 0; i <= diagIndex; i++) {  // k iterations
    ta()()(i, i) = 1.0;
  }
  ta()()(k, k) = -k;  // indexing starts from 0
  RealD nrm = 1.0 / std::sqrt(2.0 * k * (k + 1));
  ta = ta * nrm;
}

////////////////////////////////////////////////////////////////////////
// Map a su2 subgroup number to the pair of rows that are non zero
////////////////////////////////////////////////////////////////////////
static accelerator_inline void su2SubGroupIndex(int &i1, int &i2, int su2_index, GroupName::SU) {
  assert((su2_index >= 0) && (su2_index < (ncolour * (ncolour - 1)) / 2));

  int spare = su2_index;
  for (i1 = 0; spare >= (ncolour - 1 - i1); i1++) {
    spare = spare - (ncolour - 1 - i1);  // remove the Nc-1-i1 terms
  }
  i2 = i1 + 1 + spare;
}

public:
//////////////////////////////////////////////////////////////////////////////////////////
// Pull out a subgroup and project on to real coeffs x pauli basis
//////////////////////////////////////////////////////////////////////////////////////////
template <class vcplx, ONLY_IF_SU>
static void su2Extract(Lattice<iSinglet<vcplx> > &Determinant,
                       Lattice<iSU2Matrix<vcplx> > &subgroup,
                       const Lattice<iGroupMatrix<vcplx> > &source,
                       int su2_index) {
  GridBase *grid(source.Grid());
  conformable(subgroup, source);
  conformable(subgroup, Determinant);
  int i0, i1;
  su2SubGroupIndex(i0, i1, su2_index);

  autoView(subgroup_v, subgroup, AcceleratorWrite);
  autoView(source_v, source, AcceleratorRead);
  autoView(Determinant_v, Determinant, AcceleratorWrite);
  accelerator_for(ss, grid->oSites(), 1, {
    subgroup_v[ss]()()(0, 0) = source_v[ss]()()(i0, i0);
    subgroup_v[ss]()()(0, 1) = source_v[ss]()()(i0, i1);
    subgroup_v[ss]()()(1, 0) = source_v[ss]()()(i1, i0);
    subgroup_v[ss]()()(1, 1) = source_v[ss]()()(i1, i1);

    iSU2Matrix<vcplx> Sigma = subgroup_v[ss];

    Sigma = Sigma - adj(Sigma) + trace(adj(Sigma));

    subgroup_v[ss] = Sigma;

    // this should be purely real
    Determinant_v[ss] =
        Sigma()()(0, 0) * Sigma()()(1, 1) - Sigma()()(0, 1) * Sigma()()(1, 0);
  });
}

//////////////////////////////////////////////////////////////////////////////////////////
// Set matrix to one and insert a pauli subgroup
//////////////////////////////////////////////////////////////////////////////////////////
template <class vcplx, ONLY_IF_SU>
static void su2Insert(const Lattice<iSU2Matrix<vcplx> > &subgroup,
                      Lattice<iGroupMatrix<vcplx> > &dest, int su2_index) {
  GridBase *grid(dest.Grid());
  conformable(subgroup, dest);
  int i0, i1;
  su2SubGroupIndex(i0, i1, su2_index);

  dest = 1.0;  // start out with identity
  autoView(dest_v, dest, AcceleratorWrite);
  autoView(subgroup_v, subgroup, AcceleratorRead);
  accelerator_for(ss, grid->oSites(), 1, {
    dest_v[ss]()()(i0, i0) = subgroup_v[ss]()()(0, 0);
    dest_v[ss]()()(i0, i1) = subgroup_v[ss]()()(0, 1);
    dest_v[ss]()()(i1, i0) = subgroup_v[ss]()()(1, 0);
    dest_v[ss]()()(i1, i1) = subgroup_v[ss]()()(1, 1);
  });
}

///////////////////////////////////////////////
// Generate e^{ Re Tr Staple Link} dlink
//
// *** Note Staple should be appropriate linear compbination between all
// staples.
// *** If already by beta pass coefficient 1.0.
// *** This routine applies the additional 1/Nc factor that comes after trace
// in action.
//
///////////////////////////////////////////////
template <ONLY_IF_SU>
static void SubGroupHeatBath(
    GridSerialRNG &sRNG, GridParallelRNG &pRNG,
    RealD beta,  // coeff multiplying staple in action (with no 1/Nc)
    LatticeMatrix &link,
    const LatticeMatrix &barestaple,  // multiplied by action coeffs so th
    int su2_subgroup, int nheatbath, LatticeInteger &wheremask) {
  GridBase *grid = link.Grid();

  const RealD twopi = 2.0 * M_PI;

  LatticeMatrix staple(grid);

  staple = barestaple * (beta / ncolour);

  LatticeMatrix V(grid);
  V = link * staple;

  // Subgroup manipulation in the lie algebra space
  LatticeSU2Matrix u(
      grid);  // Kennedy pendleton "u" real projected normalised Sigma
  LatticeSU2Matrix uinv(grid);
  LatticeSU2Matrix ua(grid);  // a in pauli form
  LatticeSU2Matrix b(grid);   // rotated matrix after hb

  // Some handy constant fields
  LatticeComplex ones(grid);
  ones = 1.0;
  LatticeComplex zeros(grid);
  zeros = Zero();
  LatticeReal rones(grid);
  rones = 1.0;
  LatticeReal rzeros(grid);
  rzeros = Zero();
  LatticeComplex udet(grid);  // determinant of real(staple)
  LatticeInteger mask_true(grid);
  mask_true = 1;
  LatticeInteger mask_false(grid);
  mask_false = 0;

  /*
    PLB 156 P393 (1985) (Kennedy and Pendleton)

    Note: absorb "beta" into the def of sigma compared to KP paper; staple
    passed to this routine has "beta" already multiplied in

    Action linear in links h and of form:

    beta S = beta  Sum_p (1 - 1/Nc Re Tr Plaq )

    Writing Sigma = 1/Nc (beta Sigma') where sum over staples is "Sigma' "

    beta S = const - beta/Nc Re Tr h Sigma'
    = const - Re Tr h Sigma

    Decompose h and Sigma into (1, sigma_j) ; h_i real, h^2=1, Sigma_i complex
    arbitrary.

    Tr h Sigma = h_i Sigma_j Tr (sigma_i sigma_j)  = h_i Sigma_j 2 delta_ij
    Re Tr h Sigma = 2 h_j Re Sigma_j

    Normalised re Sigma_j = xi u_j

    With u_j a unit vector and U can be in SU(2);

    Re Tr h Sigma = 2 h_j Re Sigma_j = 2 xi (h.u)

    4xi^2 = Det [ Sig - Sig^dag  + 1 Tr Sigdag]
    u   = 1/2xi [ Sig - Sig^dag  + 1 Tr Sigdag]

    xi = sqrt(Det)/2;

    Write a= u h in SU(2); a has pauli decomp a_j;

    Note: Product b' xi is unvariant because scaling Sigma leaves
    normalised vector "u" fixed; Can rescale Sigma so b' = 1.
  */

  ////////////////////////////////////////////////////////
  // Real part of Pauli decomposition
  // Note a subgroup can project to zero in cold start
  ////////////////////////////////////////////////////////
  su2Extract(udet, u, V, su2_subgroup);

  //////////////////////////////////////////////////////
  // Normalising this vector if possible; else identity
  //////////////////////////////////////////////////////
  LatticeComplex xi(grid);

  LatticeSU2Matrix lident(grid);

  SU2Matrix ident = Complex(1.0);
  SU2Matrix pauli1;
  GaugeGroup<2, GroupName::SU>::generator(0, pauli1);
  SU2Matrix pauli2;
  GaugeGroup<2, GroupName::SU>::generator(1, pauli2);
  SU2Matrix pauli3;
  GaugeGroup<2, GroupName::SU>::generator(2, pauli3);
  pauli1 = timesI(pauli1) * 2.0;
  pauli2 = timesI(pauli2) * 2.0;
  pauli3 = timesI(pauli3) * 2.0;

  LatticeComplex cone(grid);
  LatticeReal adet(grid);
  adet = abs(toReal(udet));
  lident = Complex(1.0);
  cone = Complex(1.0);
  Real machine_epsilon = 1.0e-7;
  u = where(adet > machine_epsilon, u, lident);
  udet = where(adet > machine_epsilon, udet, cone);

  xi = 0.5 * sqrt(udet);        // 4xi^2 = Det [ Sig - Sig^dag  + 1 Tr Sigdag]
  u = 0.5 * u * pow(xi, -1.0);  //  u   = 1/2xi [ Sig - Sig^dag  + 1 Tr Sigdag]

  // Debug test for sanity
  uinv = adj(u);
  b = u * uinv - 1.0;
  assert(norm2(b) < 1.0e-4);

  /*
    Measure: Haar measure dh has d^4a delta(1-|a^2|)
    In polars:
    da = da0 r^2 sin theta dr dtheta dphi delta( 1 - r^2 -a0^2)
    = da0 r^2 sin theta dr dtheta dphi delta( (sqrt(1-a0^) - r)(sqrt(1-a0^) +
    r) )
    = da0 r/2 sin theta dr dtheta dphi delta( (sqrt(1-a0^) - r) )

    Action factor Q(h) dh  = e^-S[h]  dh =  e^{  xi Tr uh} dh    // beta
    enters through xi =  e^{2 xi (h.u)} dh =  e^{2 xi h0u0}.e^{2 xi h1u1}.e^{2
    xi h2u2}.e^{2 xi h3u3} dh

    Therefore for each site, take xi for that site
    i) generate  |a0|<1 with dist
    (1-a0^2)^0.5 e^{2 xi a0 } da0

    Take alpha = 2 xi  = 2 xi [ recall 2 beta/Nc unmod staple norm];
    hence 2.0/Nc factor in Chroma ] A. Generate two uniformly distributed
    pseudo-random numbers R and R', R'', R''' in the unit interval; B. Set X =
    -(ln R)/alpha, X' =-(ln R')/alpha; C. Set C = cos^2(2pi R"), with R"
    another uniform random number in [0,1] ; D. Set A = XC; E. Let d  = X'+A;
    F. If R'''^2 :> 1 - 0.5 d,  go back to A;
    G. Set a0 = 1 - d;

    Note that in step D setting B ~ X - A and using B in place of A in step E
    will generate a second independent a 0 value.
  */

  /////////////////////////////////////////////////////////
  // count the number of sites by picking "1"'s out of hat
  /////////////////////////////////////////////////////////
  Integer hit = 0;
  LatticeReal rtmp(grid);
  rtmp = where(wheremask, rones, rzeros);
  RealD numSites = sum(rtmp);
  RealD numAccepted;
  LatticeInteger Accepted(grid);
  Accepted = Zero();
  LatticeInteger newlyAccepted(grid);

  std::vector<LatticeReal> xr(4, grid);
  std::vector<LatticeReal> a(4, grid);
  LatticeReal d(grid);
  d = Zero();
  LatticeReal alpha(grid);

  //    std::cout<<GridLogMessage<<"xi "<<xi <<std::endl;
  xi = 2.0 * xi;
  alpha = toReal(xi);

  do {
    // A. Generate two uniformly distributed pseudo-random numbers R and R',
    // R'', R''' in the unit interval;
    random(pRNG, xr[0]);
    random(pRNG, xr[1]);
    random(pRNG, xr[2]);
    random(pRNG, xr[3]);

    // B. Set X = - ln R/alpha, X' = -ln R'/alpha
    xr[1] = -log(xr[1]) / alpha;
    xr[2] = -log(xr[2]) / alpha;

    // C. Set C = cos^2(2piR'')
    xr[3] = cos(xr[3] * twopi);
    xr[3] = xr[3] * xr[3];

    LatticeReal xrsq(grid);

    // D. Set A = XC;
    // E. Let d  = X'+A;
    xrsq = xr[2] + xr[1] * xr[3];

    d = where(Accepted, d, xr[2] + xr[1] * xr[3]);

    // F. If R'''^2 :> 1 - 0.5 d,  go back to A;
    LatticeReal thresh(grid);
    thresh = 1.0 - d * 0.5;
    xrsq = xr[0] * xr[0];
    LatticeInteger ione(grid);
    ione = 1;
    LatticeInteger izero(grid);
    izero = Zero();

    newlyAccepted = where(xrsq < thresh, ione, izero);
    Accepted = where(newlyAccepted, newlyAccepted, Accepted);
    Accepted = where(wheremask, Accepted, izero);

    // FIXME need an iSum for integer to avoid overload on return type??
    rtmp = where(Accepted, rones, rzeros);
    numAccepted = sum(rtmp);

    hit++;

  } while ((numAccepted < numSites) && (hit < nheatbath));

  // G. Set a0 = 1 - d;
  a[0] = Zero();
  a[0] = where(wheremask, 1.0 - d, a[0]);

  //////////////////////////////////////////
  //    ii) generate a_i uniform on two sphere radius (1-a0^2)^0.5
  //////////////////////////////////////////

  LatticeReal a123mag(grid);
  a123mag = sqrt(abs(1.0 - a[0] * a[0]));

  LatticeReal cos_theta(grid);
  LatticeReal sin_theta(grid);
  LatticeReal phi(grid);

  random(pRNG, phi);
  phi = phi * twopi;  // uniform in [0,2pi]
  random(pRNG, cos_theta);
  cos_theta = (cos_theta * 2.0) - 1.0;  // uniform in [-1,1]
  sin_theta = sqrt(abs(1.0 - cos_theta * cos_theta));

  a[1] = a123mag * sin_theta * cos(phi);
  a[2] = a123mag * sin_theta * sin(phi);
  a[3] = a123mag * cos_theta;

  ua = toComplex(a[0]) * ident + toComplex(a[1]) * pauli1 +
       toComplex(a[2]) * pauli2 + toComplex(a[3]) * pauli3;

  b = 1.0;
  b = where(wheremask, uinv * ua, b);
  su2Insert(b, V, su2_subgroup);

  // mask the assignment back based on Accptance
  link = where(Accepted, V * link, link);

  //////////////////////////////
  // Debug Checks
  // SU2 check
  LatticeSU2Matrix check(grid);  // rotated matrix after hb
  u = Zero();
  check = ua * adj(ua) - 1.0;
  check = where(Accepted, check, u);
  assert(norm2(check) < 1.0e-4);

  check = b * adj(b) - 1.0;
  check = where(Accepted, check, u);
  assert(norm2(check) < 1.0e-4);

  LatticeMatrix Vcheck(grid);
  Vcheck = Zero();
  Vcheck = where(Accepted, V * adj(V) - 1.0, Vcheck);
  //    std::cout<<GridLogMessage << "SU3 check " <<norm2(Vcheck)<<std::endl;
  assert(norm2(Vcheck) < 1.0e-4);

  // Verify the link stays in SU(3)
  //    std::cout<<GridLogMessage <<"Checking the modified link"<<std::endl;
  Vcheck = link * adj(link) - 1.0;
  assert(norm2(Vcheck) < 1.0e-4);
  /////////////////////////////////
}

template <ONLY_IF_SU>
static void testGenerators(GroupName::SU) {
  Matrix ta;
  Matrix tb;
  std::cout << GridLogMessage
            << "Fundamental - Checking trace ta tb is 0.5 delta_ab"
            << std::endl;
  for (int a = 0; a < AdjointDimension; a++) {
    for (int b = 0; b < AdjointDimension; b++) {
      generator(a, ta);
      generator(b, tb);
      Complex tr = TensorRemove(trace(ta * tb));
      std::cout << GridLogMessage << "(" << a << "," << b << ") =  " << tr
                << std::endl;
      if (a == b) assert(abs(tr - Complex(0.5)) < 1.0e-6);
      if (a != b) assert(abs(tr) < 1.0e-6);
    }
    std::cout << GridLogMessage << std::endl;
  }
  std::cout << GridLogMessage << "Fundamental - Checking if hermitian"
            << std::endl;
  for (int a = 0; a < AdjointDimension; a++) {
    generator(a, ta);
    std::cout << GridLogMessage << a << std::endl;
    assert(norm2(ta - adj(ta)) < 1.0e-6);
  }
  std::cout << GridLogMessage << std::endl;

  std::cout << GridLogMessage << "Fundamental - Checking if traceless"
            << std::endl;
  for (int a = 0; a < AdjointDimension; a++) {
    generator(a, ta);
    Complex tr = TensorRemove(trace(ta));
    std::cout << GridLogMessage << a << " " << std::endl;
    assert(abs(tr) < 1.0e-6);
  }
  std::cout << GridLogMessage << std::endl;
}


template <int N, class vtype>
static Lattice<iScalar<iScalar<iMatrix<vtype, N> > > >
ProjectOnGeneralGroup(const Lattice<iScalar<iScalar<iMatrix<vtype, N> > > > &Umu, GroupName::SU) {
  return ProjectOnGroup(Umu);
}

template <class vtype>
accelerator_inline static iScalar<vtype> ProjectOnGeneralGroup(const iScalar<vtype> &r, GroupName::SU) {
  return ProjectOnGroup(r);
}

template <class vtype, int N>
accelerator_inline static iVector<vtype,N> ProjectOnGeneralGroup(const iVector<vtype,N> &r, GroupName::SU) {
  return ProjectOnGroup(r);
}

template <class vtype,int N, typename std::enable_if< GridTypeMapper<vtype>::TensorLevel == 0 >::type * =nullptr>
accelerator_inline static iMatrix<vtype,N> ProjectOnGeneralGroup(const iMatrix<vtype,N> &arg, GroupName::SU) {
  return ProjectOnGroup(arg);
}

template <typename LatticeMatrixType>
static void taProj(const LatticeMatrixType &in, LatticeMatrixType &out, GroupName::SU) {
  out = Ta(in);
}

/*
 * Fundamental rep gauge xform
 */
template<typename Fundamental,typename GaugeMat>
static void GaugeTransformFundamental( Fundamental &ferm, GaugeMat &g){
  GridBase *grid = ferm._grid;
  conformable(grid,g._grid);
  ferm = g*ferm;
}
/*
 * Adjoint rep gauge xform
 */

template<typename Gimpl>
static void GaugeTransform(typename Gimpl::GaugeField &Umu, typename Gimpl::GaugeLinkField &g){
  GridBase *grid = Umu.Grid();
  conformable(grid,g.Grid());

  typename Gimpl::GaugeLinkField U(grid);
  typename Gimpl::GaugeLinkField ag(grid); ag = adj(g);

  for(int mu=0;mu<Nd;mu++){
    U= PeekIndex<LorentzIndex>(Umu,mu);
    U = g*U*Gimpl::CshiftLink(ag, mu, 1); //BC-aware
    PokeIndex<LorentzIndex>(Umu,U,mu);
  }
}
template<typename Gimpl>
static void GaugeTransform( std::vector<typename Gimpl::GaugeLinkField> &U, typename Gimpl::GaugeLinkField &g){
  GridBase *grid = g.Grid();
  typename Gimpl::GaugeLinkField ag(grid); ag = adj(g);
  for(int mu=0;mu<Nd;mu++){
    U[mu] = g*U[mu]*Gimpl::CshiftLink(ag, mu, 1); //BC-aware
  }
}
template<typename Gimpl>
static void RandomGaugeTransform(GridParallelRNG &pRNG, typename Gimpl::GaugeField &Umu, typename Gimpl::GaugeLinkField &g){
  LieRandomize(pRNG,g,1.0);
  GaugeTransform<Gimpl>(Umu,g);
}

