/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_pvdagm_coherence.cc

    Copyright (C) 2026

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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

// Local-coherence census of the COARSE operator's low singular vectors.
// (Fig 4 of arXiv:2409.03904 transplanted one level down.)
//
// Inputs: the raw fine subspace cache AND the coarse right-singular
// vectors saved by Example_pvdagm_census (CENSUS_EVEC_FILE).
//
// Measurements:
//  CHECK 0  u/v alignment per mode: u_k = A v_k/||A v_k||, align=|<u,v>|.
//           Low-sector non-normality meter AND the "2-for-1 coupon" test
//           (u adds deflation span iff align is small). Also re-measures
//           sigma_k = ||A v_k|| as an evec-file integrity check.
//  CHECK 1  Band overlap: |projection of v_k onto span(psi_c)|^2 where
//           psi_c = coarse images of the 60 RAW fine null vectors.
//           Detached-vs-submerged Ritz band question.
//  CHECK A  Self-coherence: level-2 blocks (BLOCK2) with basis = evecs
//           0..K-1 (block-orthonormalised); completeness of held-out
//           evecs K..NEV-1. Local coherence of the coarse op's own tail.
//  CHECK B  Capture by the available basis: same blocks, basis = psi_c;
//           completeness of every evec. What a 3-level would have had.
//  GS census (inside every block orthonormalisation): per-vector global
//           remainder fraction and MIN-over-blocks remainder fraction.
//           Detects finite-precision rank loss (near-parallel vectors
//           within a block -> normalised noise). Values ~0 flag fiction
//           in the "orthonormal" basis.
//  Idempotency assertion ||(PP+)^2 v - PP+ v||/||PP+ v|| ~ eps on the
//           first probe of each check: catches non-orthonormal-basis and
//           rank-loss corruption in the checks themselves.
//
// Level-2 machinery uses the free block primitives directly on
// std::vector<CoarseVector>, so basis rank K is a RUNTIME parameter
// (no Aggregation template instantiation at level 2).
//
// Env:
//   MASS, SUBSPACE_FILE          as elsewhere (cache MUST exist)
//   CENSUS_EVEC_FILE             REQUIRED: evecs from the census run
//   NEV                          number of evecs to load (default 96)
//   NBASIS2                      K, basis size for CHECK A (default 48)
//   BLOCK                        fine->coarse blocking; MUST match the
//                                census run that made the evecs (4.4.4.4)
//   BLOCK2                       coarse->coarsecoarse blocking (2.2.2.2)

#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>

using namespace std;
using namespace Grid;

RealD mass = 0.00078;
int   Nev     = 96;
int   Nbasis2 = 48;

void ParseEnvironment(void)
{
  if(getenv("MASS"))    mass    = atof(getenv("MASS"));
  if(getenv("NEV"))     Nev     = atoi(getenv("NEV"));
  if(getenv("NBASIS2")) Nbasis2 = atoi(getenv("NBASIS2"));
  std::cout << GridLogMessage << "PARAM: MASS    " << mass    << std::endl;
  std::cout << GridLogMessage << "PARAM: NEV     " << Nev     << std::endl;
  std::cout << GridLogMessage << "PARAM: NBASIS2 " << Nbasis2 << std::endl;
}

template <class Field>
void loadFields(std::vector<Field> &v, std::string const fname){
#ifdef HAVE_LIME
  std::cout << Grid::GridLogMessage << "Loading " << v.size() << " fields from: " << fname << std::endl;
  Grid::emptyUserRecord record;
  Grid::ScidacReader SR;
  SR.open(fname);
  for (int k = 0; k < (int)v.size(); k++)
    SR.readScidacFieldRecord(v[k], record);
  SR.close();
#endif
}

template<class Matrix,class Field>
class PVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  Matrix &_PV;
public:
  PVdagMLinearOperator(Matrix &Mat,Matrix &PV): _Mat(Mat),_PV(PV) {};
  void OpDiag (const Field &in, Field &out) {    assert(0);  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){    assert(0);  };
  void Op     (const Field &in, Field &out){
    Field tmp(in.Grid());
    _Mat.M(in,tmp);
    _PV.Mdag(tmp,out);
  }
  void AdjOp     (const Field &in, Field &out){
    Field tmp(in.Grid());
    _PV.M(in,tmp);
    _Mat.Mdag(tmp,out);
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
    HermOp(in,out);
    ComplexD dot = innerProduct(in,out);
    n1=real(dot); n2=norm2(out);
  }
  void HermOp(const Field &in, Field &out){
    Field tmp(in.Grid());
    Op(in,tmp);
    AdjOp(tmp,out);
  }
};

//////////////////////////////////////////////////////////////////////
// Block Gram-Schmidt with remainder census.
// Mirrors blockOrthonormalize (one GS pass + normalise) but reports,
// per vector: the GLOBAL remainder fraction ||v_k^perp||^2/||v_k||^2
// and the MIN over blocks of the same ratio. A second silent pass
// tightens orthonormality for the completeness measurements.
//////////////////////////////////////////////////////////////////////
template<class Field, class CField>
void blockGSCensus(std::vector<Field> &basis, GridBase *ccGrid, const std::string tag)
{
  int n = basis.size();
  CField ip(ccGrid);
  CField n0(ccGrid), n1(ccGrid);

  std::cout << GridLogMessage << "=== block GS census [" << tag << "] : " << n
            << " vectors, blocks -> " << ccGrid->GlobalDimensions() << std::endl;

  for(int v=0; v<n; v++){
    RealD before = norm2(basis[v]);
    blockInnerProduct(n0, basis[v], basis[v]);
    for(int u=0; u<v; u++){
      blockInnerProductD(ip, basis[u], basis[v]);
      ip = -ip;
      blockZAXPY(basis[v], ip, basis[u], basis[v]);
    }
    RealD after = norm2(basis[v]);
    blockInnerProduct(n1, basis[v], basis[v]);

    // min over blocks of remainder fraction n1/n0
    RealD locmin = 1.0e60;
    {
      typedef typename CField::scalar_object sobj;
      for(int64_t s=0; s<ccGrid->lSites(); s++){
        Coordinate lcoor(ccGrid->_ndimension);
        ccGrid->LocalIndexToLocalCoor(s,lcoor);
        sobj s0, s1;
        peekLocalSite(s0, n0, lcoor);
        peekLocalSite(s1, n1, lcoor);
        RealD r0 = real(TensorRemove(s0));
        RealD r1 = real(TensorRemove(s1));
        if ( r0 > 0.0 ) {
          RealD frac = r1/r0;
          if (frac < locmin) locmin = frac;
        }
      }
      RealD neg = -locmin;
      ccGrid->GlobalMax(neg);
      locmin = -neg;
    }

    std::cout << GridLogMessage << "GS[" << tag << "] vec " << v
              << "  global remainder " << after/before
              << "  min-block remainder " << locmin
              << ( locmin < 1.0e-10 ? "   <== RANK LOSS FLAG" : "" ) << std::endl;

    blockNormalise(ip, basis[v]);
  }
  // Second (silent) pass for numerical hygiene of downstream projections
  for(int v=0; v<n; v++){
    for(int u=0; u<v; u++){
      blockInnerProductD(ip, basis[u], basis[v]);
      ip = -ip;
      blockZAXPY(basis[v], ip, basis[u], basis[v]);
    }
    blockNormalise(ip, basis[v]);
  }
}

//////////////////////////////////////////////////////////////////////
// P P^dag v for a block-orthonormal basis, via the block primitives.
// Runtime basis size; no Aggregation template needed.
//////////////////////////////////////////////////////////////////////
template<class Field, class CField>
void blockPPdag(const std::vector<Field> &q, GridBase *ccGrid, const Field &v, Field &w)
{
  CField ip(ccGrid);
  w = Zero();
  for(int j=0; j<(int)q.size(); j++){
    blockInnerProduct(ip, q[j], v);
    blockZAXPY(w, ip, q[j], w);
  }
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);
  ParseEnvironment();

  const int Ls=24;
  RealD M5=1.8, b=1.5, c=0.5;
  const int nbasis = 60;

  std::vector<int> lat_size {48, 48, 48, 96};

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(lat_size, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  // Fine->coarse blocking: MUST match the census run that made the evecs.
  Coordinate clatt = lat_size;
  Coordinate Block({4,4,4,4});
  if ( getenv("BLOCK") ) {
    GridCmdOptionIntVector(std::string(getenv("BLOCK")),Block);
    GRID_ASSERT(Block.size()==4);
  }
  for(int d=0;d<clatt.size();d++){
    GRID_ASSERT(lat_size[d] % Block[d] == 0);
    clatt[d] = lat_size[d]/Block[d];
  }
  std::cout << GridLogMessage << "Block " << Block << "  coarse lattice " << clatt << std::endl;

  GridCartesian *Coarse4d =  SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian *Coarse5d =  SpaceTimeGrid::makeFiveDimGrid(1,Coarse4d);

  // Coarse->coarsecoarse blocking for the coherence checks
  Coordinate cclatt = clatt;
  Coordinate Block2({2,2,2,2});
  if ( getenv("BLOCK2") ) {
    GridCmdOptionIntVector(std::string(getenv("BLOCK2")),Block2);
    GRID_ASSERT(Block2.size()==4);
  }
  for(int d=0;d<cclatt.size();d++){
    GRID_ASSERT(clatt[d] % Block2[d] == 0);
    cclatt[d] = clatt[d]/Block2[d];
  }
  std::cout << GridLogMessage << "Block2 " << Block2 << "  coarsecoarse lattice " << cclatt << std::endl;

  GridCartesian *CC4d =  SpaceTimeGrid::makeFourDimGrid(cclatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian *CC5d =  SpaceTimeGrid::makeFiveDimGrid(1,CC4d);

  GridParallelRNG RNG5(FGrid);  RNG5.SeedFixedIntegers({5,6,7,8});

  LatticeGaugeField Umu(UGrid);
  std::cout << GridLogMessage << "Reading gauge field" << std::endl;
  FieldMetaData header;
  std::string file("/ccs/home/poare/ckpoint_lat.1000");
  NerscIO::readConfiguration(Umu,header,file);

  MobiusFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c);
  MobiusFermionD Dpv (Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,1.0, M5,b,c);

  typedef PVdagMLinearOperator<MobiusFermionD,LatticeFermionD>       PVdagM_t;
  typedef GeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector                          CoarseVector;
  // Grid index contraction is positional (colour/spin/lorentz order is meaningful),
  // so each MG projection adds one index to the tensor nest rather than reusing a slot:
  // innerProduct(CoarseSiteObj,CoarseSiteObj) returns iScalar<vTComplex>, one level
  // deeper than the fine-level vTComplex (same convention as Example_pvdagm_3level.cc).
  typedef Lattice<iScalar<vTComplex> >                               CoarseCoarseScalar;
  typedef Aggregation<vSpinColourVector,vTComplex,nbasis>            Subspace;

  PVdagM_t PVdagM(Ddwf,Dpv);
  NextToNearestStencilGeometry5D geom(Coarse5d);

  //////////////////////////////////////////////////////////////////////
  // Load the RAW subspace (cache REQUIRED; no generation here)
  //////////////////////////////////////////////////////////////////////
  std::string subspace_file = "/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/subspace_nb"
                              + std::to_string(nbasis) + ".scidac";
  if ( getenv("SUBSPACE_FILE") ) subspace_file = std::string(getenv("SUBSPACE_FILE"));

  uint64_t file_exists = 0;
  if ( UGrid->IsBoss() ) { std::ifstream f(subspace_file); file_exists = f.good() ? 1 : 0; }
  UGrid->GlobalSum(file_exists);
  if ( !file_exists ) {
    std::cout << GridLogMessage << "FATAL: subspace cache not found: " << subspace_file << std::endl;
    GRID_ASSERT(file_exists);
  }

  const int cb = 0;
  Subspace AggregatesGCR(Coarse5d,FGrid,cb);
  loadFields(AggregatesGCR.subspace, subspace_file);

  // RAW copy before CoarsenOperator block-orthogonalises in place
  std::vector<LatticeFermionD> rawNull(nbasis, FGrid);
  for (int k = 0; k < nbasis; k++) rawNull[k] = AggregatesGCR.subspace[k];

  LittleDiracOperator LittleDiracOpPV(geom,FGrid,Coarse5d);
  LittleDiracOpPV.CoarsenOperator(PVdagM, AggregatesGCR);
  NonHermitianLinearOperator<LittleDiracOperator,CoarseVector> LinOpCoarse(LittleDiracOpPV);

  // Coarse images of the raw null vectors
  std::vector<CoarseVector> psi_c(nbasis, Coarse5d);
  for (int k = 0; k < nbasis; k++)
    AggregatesGCR.ProjectToSubspace(psi_c[k], rawNull[k]);

  //////////////////////////////////////////////////////////////////////
  // Load census singular vectors
  //////////////////////////////////////////////////////////////////////
  GRID_ASSERT(getenv("CENSUS_EVEC_FILE"));
  std::string evec_file(getenv("CENSUS_EVEC_FILE"));
  std::vector<CoarseVector> evec(Nev, Coarse5d);
  loadFields(evec, evec_file);

  GRID_ASSERT(Nbasis2 < Nev);

  //////////////////////////////////////////////////////////////////////
  // CHECK 0: sigma re-measurement and u/v alignment per mode
  //////////////////////////////////////////////////////////////////////
  std::cout << GridLogMessage << "=================================================" << std::endl;
  std::cout << GridLogMessage << "CHECK 0: sigma and left/right alignment per mode" << std::endl;
  std::cout << GridLogMessage << "=================================================" << std::endl;
  {
    CoarseVector Av(Coarse5d);
    for(int k=0;k<Nev;k++){
      LinOpCoarse.Op(evec[k],Av);
      RealD nv  = norm2(evec[k]);
      RealD nAv = norm2(Av);
      RealD sigma = std::sqrt(nAv/nv);
      ComplexD uv = innerProduct(Av,evec[k]);
      RealD align = abs(uv)/std::sqrt(nAv*nv);   // |<u|v>|, u = Av/||Av||
      std::cout << GridLogMessage << "COHERENCE: mode " << k
                << "  sigma = " << sigma
                << "  |<u|v>| = " << align << std::endl;
    }
  }

  //////////////////////////////////////////////////////////////////////
  // CHECK 1: band overlap of each evec with span(psi_c)  (GLOBAL span)
  //////////////////////////////////////////////////////////////////////
  std::cout << GridLogMessage << "=================================================" << std::endl;
  std::cout << GridLogMessage << "CHECK 1: evec overlap with GLOBAL span of 60 psi_c" << std::endl;
  std::cout << GridLogMessage << "=================================================" << std::endl;
  {
    // Global (not block) Gram-Schmidt of copies of psi_c
    std::vector<CoarseVector> Q(nbasis, Coarse5d);
    for(int k=0;k<nbasis;k++) Q[k]=psi_c[k];
    for(int k=0;k<nbasis;k++){
      for(int j=0;j<k;j++){
        ComplexD ip = innerProduct(Q[j],Q[k]);
        Q[k] = Q[k] - ip*Q[j];
      }
      RealD nq = norm2(Q[k]);
      GRID_ASSERT(nq>0.0);
      Q[k] = Q[k] * (1.0/std::sqrt(nq));
    }
    for(int k=0;k<Nev;k++){
      RealD ov = 0.0;
      for(int j=0;j<nbasis;j++){
        ComplexD ip = innerProduct(Q[j],evec[k]);
        ov += real(ip*conjugate(ip));
      }
      std::cout << GridLogMessage << "COHERENCE: band overlap mode " << k
                << " = " << ov/norm2(evec[k]) << std::endl;
    }
  }

  //////////////////////////////////////////////////////////////////////
  // CHECK A: self-coherence -- basis = evecs 0..K-1 in BLOCK2 blocks,
  // probes = held-out evecs K..Nev-1
  //////////////////////////////////////////////////////////////////////
  std::cout << GridLogMessage << "=================================================" << std::endl;
  std::cout << GridLogMessage << "CHECK A: tail self-coherence, K=" << Nbasis2 << " basis evecs" << std::endl;
  std::cout << GridLogMessage << "=================================================" << std::endl;
  {
    std::vector<CoarseVector> basisA(Nbasis2, Coarse5d);
    for(int k=0;k<Nbasis2;k++) basisA[k]=evec[k];
    blockGSCensus<CoarseVector,CoarseCoarseScalar>(basisA, CC5d, "A:evecs");

    CoarseVector w(Coarse5d), w2(Coarse5d);
    // Idempotency assertion on first probe
    blockPPdag<CoarseVector,CoarseCoarseScalar>(basisA, CC5d, evec[Nbasis2], w);
    blockPPdag<CoarseVector,CoarseCoarseScalar>(basisA, CC5d, w, w2);
    w2 = w2 - w;
    RealD idem = std::sqrt(norm2(w2)/norm2(w));
    std::cout << GridLogMessage << "COHERENCE: CHECK A idempotency = " << idem
              << "  (expect ~1e-14; O(1) => basis corrupt)" << std::endl;

    RealD sum=0.0;
    for(int k=Nbasis2;k<Nev;k++){
      blockPPdag<CoarseVector,CoarseCoarseScalar>(basisA, CC5d, evec[k], w);
      RealD comp = norm2(w)/norm2(evec[k]);
      sum += comp;
      std::cout << GridLogMessage << "COHERENCE: A completeness mode " << k
                << " = " << comp << std::endl;
    }
    std::cout << GridLogMessage << "COHERENCE: A mean completeness (modes " << Nbasis2
              << ".." << Nev-1 << ") = " << sum/(Nev-Nbasis2) << std::endl;
  }

  //////////////////////////////////////////////////////////////////////
  // CHECK B: capture by the available basis -- basis = 60 psi_c in
  // BLOCK2 blocks, probes = ALL evecs
  //////////////////////////////////////////////////////////////////////
  std::cout << GridLogMessage << "=================================================" << std::endl;
  std::cout << GridLogMessage << "CHECK B: capture by blocked psi_c basis" << std::endl;
  std::cout << GridLogMessage << "=================================================" << std::endl;
  {
    std::vector<CoarseVector> basisB(nbasis, Coarse5d);
    for(int k=0;k<nbasis;k++) basisB[k]=psi_c[k];
    blockGSCensus<CoarseVector,CoarseCoarseScalar>(basisB, CC5d, "B:psi_c");

    CoarseVector w(Coarse5d), w2(Coarse5d);
    blockPPdag<CoarseVector,CoarseCoarseScalar>(basisB, CC5d, evec[0], w);
    blockPPdag<CoarseVector,CoarseCoarseScalar>(basisB, CC5d, w, w2);
    w2 = w2 - w;
    RealD idem = std::sqrt(norm2(w2)/norm2(w));
    std::cout << GridLogMessage << "COHERENCE: CHECK B idempotency = " << idem
              << "  (expect ~1e-14; O(1) => basis corrupt)" << std::endl;

    RealD sumLow=0.0, sumHigh=0.0; int nLow=0, nHigh=0;
    for(int k=0;k<Nev;k++){
      blockPPdag<CoarseVector,CoarseCoarseScalar>(basisB, CC5d, evec[k], w);
      RealD comp = norm2(w)/norm2(evec[k]);
      if (k<nbasis) { sumLow+=comp; nLow++; } else { sumHigh+=comp; nHigh++; }
      std::cout << GridLogMessage << "COHERENCE: B completeness mode " << k
                << " = " << comp << std::endl;
    }
    std::cout << GridLogMessage << "COHERENCE: B mean completeness modes 0.." << nbasis-1
              << " = " << sumLow/nLow << std::endl;
    if (nHigh>0)
      std::cout << GridLogMessage << "COHERENCE: B mean completeness modes " << nbasis
                << ".." << Nev-1 << " = " << sumHigh/nHigh << std::endl;
  }

  std::cout << GridLogMessage << "=================================================" << std::endl;
  std::cout << GridLogMessage << "COHERENCE SUMMARY:" << std::endl;
  std::cout << GridLogMessage << "  CHECK 0 |<u|v>| ~ 1      : low sector effectively normal; 2-for-1 coupon worthless" << std::endl;
  std::cout << GridLogMessage << "  CHECK 1 high for low k   : Ritz band submerged in continuum (not detached)" << std::endl;
  std::cout << GridLogMessage << "  CHECK A ~0.98            : local coherence manifests at coarse level (3-level rescuable)" << std::endl;
  std::cout << GridLogMessage << "  CHECK A ~0.8-0.9         : capture deficit is physical (effective-blocking law)" << std::endl;
  std::cout << GridLogMessage << "  CHECK B << CHECK A       : psi_c basis fails to exploit coherence that exists" << std::endl;
  std::cout << GridLogMessage << "  GS min-block remainder ~0: rank loss; orthogonalisation partly fictional" << std::endl;
  std::cout << GridLogMessage << "=================================================" << std::endl;

  std::cout << GridLogMessage << "Done" << std::endl;
  Grid_finalize();
  return 0;
}
