#include <Grid/Grid.h>

using namespace Grid;

template <typename T>
bool has_correct_group_block_structure(const T& U) {
  std::cout << GridLogMessage << "Checking the structure is " << std::endl;
  std::cout << GridLogMessage << "U  =  (   W    X   )  " << std::endl;
  std::cout << GridLogMessage << "      (  -X^*  W^* )  " << std::endl;
  std::cout << GridLogMessage << std::endl;

  const int nsp = Nc / 2;
  Complex i(0., 1.);
  for (int c1 = 0; c1 < nsp; c1++)  // check on W
  {
    for (int c2 = 0; c2 < nsp; c2++) {
      auto W = PeekIndex<ColourIndex>(U, c1, c2);
      auto Wstar = PeekIndex<ColourIndex>(U, c1 + nsp, c2 + nsp);
      auto Ww = conjugate(Wstar);
      auto amizero = sum(W - Ww);
      auto amizeroo = TensorRemove(amizero);
      assert(amizeroo.real() < 10e-6);
      amizeroo *= i;
      assert(amizeroo.real() < 10e-6);
    }
  }

  for (int c1 = 0; c1 < nsp; c1++) {
    for (int c2 = 0; c2 < nsp; c2++) {
      auto X = PeekIndex<ColourIndex>(U, c1, c2 + nsp);
      auto minusXstar = PeekIndex<ColourIndex>(U, c1 + nsp, c2);
      auto minusXx = conjugate(minusXstar);
      auto amizero = sum(X + minusXx);
      auto amizeroo = TensorRemove(amizero);
      assert(amizeroo.real() < 10e-6);
      amizeroo *= i;
      assert(amizeroo.real() < 10e-6);
    }
  }
  return true;
};

template <typename T>
bool is_element_of_sp2n_group(const T& U) {
  LatticeColourMatrixD aux(U.Grid());
  LatticeColourMatrixD identity(U.Grid());
  identity = 1.0;
  LatticeColourMatrixD Omega(U.Grid());
  Sp<Nc>::Omega(Omega);

  std::cout << GridLogMessage << "Check matrix is non-zero " << std::endl;
  assert(norm2(U) > 1e-8);

  std::cout << GridLogMessage << "Unitary check" << std::endl;
  aux = U * adj(U) - identity;
  std::cout << GridLogMessage << "U adjU - 1 = " << norm2(aux) << std::endl;
  assert(norm2(aux) < 1e-8);

  aux = Omega - (U * Omega * transpose(U));
  std::cout << GridLogMessage << "Omega - U Omega transpose(U) = " << norm2(aux)
            << std::endl;
  assert(norm2(aux) < 1e-8);

  std::cout << GridLogMessage
            << "|Det| = " << norm2(Determinant(U)) / U.Grid()->gSites()
            << std::endl;
  assert(norm2(Determinant(U)) / U.Grid()->gSites() - 1 < 1e-8);

  return has_correct_group_block_structure(U);
}

int main (int argc, char **argv)
{
    Grid_init(&argc,&argv);
    
    Coordinate latt_size   = GridDefaultLatt();
    Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();

    GridCartesian             Grid(latt_size,simd_layout,mpi_layout);
    GridRedBlackCartesian     RBGrid(&Grid);
    
    LatticeGaugeField Umu(&Grid);
    LatticeColourMatrixD U(&Grid);
    
    std::vector<int> pseeds({1,2,3,4,5});
    std::vector<int> sseeds({6,7,8,9,10});
    GridParallelRNG  pRNG(&Grid); pRNG.SeedFixedIntegers(pseeds);
    GridSerialRNG    sRNG;       sRNG.SeedFixedIntegers(sseeds);
    
    std::cout << GridLogMessage << "Checking Cold Configuration " << std::endl;
    Sp<Nc>::ColdConfiguration(pRNG,Umu);
    U = PeekIndex<LorentzIndex>(Umu,1);
    assert(is_element_of_sp2n_group(U));
    
    std::cout << GridLogMessage << "Checking Hot Configuration" << std::endl;
    Sp<Nc>::HotConfiguration(pRNG,Umu);
    U = PeekIndex<LorentzIndex>(Umu,1);
    assert(is_element_of_sp2n_group(U));
    
    std::cout << GridLogMessage << "Checking Tepid Configuration" << std::endl;
    Sp<Nc>::TepidConfiguration(pRNG,Umu);
    U = PeekIndex<LorentzIndex>(Umu,1);
    assert(is_element_of_sp2n_group(U));
    
    Grid_finalize();


    }

