#include <Grid/Grid.h>

#define verbose 1

using namespace Grid;

static void antisymm_base() {

  const int this_nc = 6;
  const int this_n = this_nc/2;
  const int this_irrep_dim = Sp_TwoIndex<this_nc, AntiSymmetric>::Dimension;
  const int this_algebra_dim = Sp<this_nc>::AlgebraDimension;
  typedef Sp_TwoIndex<this_nc, AntiSymmetric>::iGroupMatrix<Complex> Matrix;
  typedef Sp_TwoIndex<this_nc, AntiSymmetric>::iGroupTwoIndexMatrix<Complex> ASMatrix;
  
  Matrix Omega;
  Matrix eij_a;
  Matrix eij_b;
  Matrix eij_c;
  Matrix e_sum;
  Omega = Zero();
  for (int i = 0; i < this_n; i++)
  {
    Omega()()(i, this_n + i) = 1.;
    Omega()()(this_n + i, i) = -1;
  }

  RealD realA;
  RealD realB;
  std::cout << GridLogMessage << "2as dimension is " << this_irrep_dim << std::endl;
  std::cout << GridLogMessage << "algebra dimension is " << this_algebra_dim << std::endl;
  realA = Sp_TwoIndex<this_nc, AntiSymmetric>::Dimension + Sp_TwoIndex<this_nc, Symmetric>::Dimension;
  assert ( realA == this_nc * this_nc - 1); // Nc x Nc = dim(2indxS) + dim(2indxAS) + dim(singlet)

  std::cout << GridLogMessage << "checking base is antisymmetric " << std::endl;
  for (int a=0; a < this_irrep_dim; a++)
  {
      Sp_TwoIndex<this_nc, AntiSymmetric>::base(a, eij_c);
      e_sum = eij_c + transpose(eij_c);
      std::cout << GridLogMessage << "e_ab + e_ab^T " << norm2(e_sum) << std::endl;
      assert(norm2(e_sum) < 1e-8);
        
  }
  std::cout << GridLogMessage << "Checking Tr (e^(ab) Omega ) = 0 and Tr (e^(ab) e^(cd) = delta^((ab)(cd)) ) " << std::endl;
  for (int a=0; a < Sp_TwoIndex<this_nc, AntiSymmetric>::Dimension; a++) {
      Sp_TwoIndex<this_nc, AntiSymmetric>::base(a, eij_a);
      realA = norm2(trace(Omega*eij_a));
      std::cout << GridLogMessage << "Omega trace for (ab) = " << a << std::endl;
      assert(realA == 0);
      for (int b=0; b < Sp_TwoIndex<this_nc, AntiSymmetric>::Dimension; b++) {
          Sp_TwoIndex<this_nc, AntiSymmetric>::base(b, eij_b);
          auto d_ab = TensorRemove(trace(eij_a * eij_b));
#if verbose
          std::cout << GridLogMessage << "Tr( e_{ab=" << a << "} e_{cd=" << b << "} ) = " << d_ab << std::endl;
#endif
          std::cout << GridLogMessage << "Orthonormality for (ab) = " << a << std::endl;
          if (a==b) {
              assert(real(d_ab)+1 < 1e-8);
              assert(imag(d_ab) < 1e-8);
          } else {
              assert(real(d_ab) < 1e-8);
              assert(imag(d_ab) < 1e-8);
          }
      }
  }

  int sum = 0;
  int sum_im = 0;
  Vector<Matrix> ta_fund(this_algebra_dim);
  Vector<Matrix> eij(this_irrep_dim);
  Matrix tmp_l;
  Matrix tmp_r;
  for (int n = 0; n < this_algebra_dim; n++)
  {
      Sp<this_nc>::generator(n, ta_fund[n]);
  }
  for (int a = 0; a < this_irrep_dim; a++)
  {
      Sp_TwoIndex<this_nc, AntiSymmetric>::base(a, eij[a]);
  }
  for (int gen_id = 0; gen_id < this_algebra_dim; gen_id++)
  {
      Complex iTr;
      sum = 0;
      sum_im = 0;
      std::cout << GridLogMessage <<  "generator number " << gen_id << std::endl;
      for (int a = 0; a < this_irrep_dim; a++)
      {
      
          tmp_l = adj(eij[a])*ta_fund[gen_id]*eij[a];
          tmp_r = adj(eij[a])*eij[a]*transpose(ta_fund[gen_id]);
#if verbose
          std::cout << GridLogMessage << " as_indx = " << a << " eDag T_F e = " << std::endl << tmp_l << std::endl;
          std::cout << GridLogMessage << " as_indx = " << a << " eDag e T_F^T = " << std::endl << tmp_r << std::endl;
#endif
          std::cout << GridLogMessage << " as_indx = " << a << " Tr(sum) = " << TensorRemove(trace(tmp_l+tmp_r)) << std::endl;
          sum += real(TensorRemove(trace(tmp_l+tmp_r)));
          sum_im += imag(TensorRemove(trace(tmp_l+tmp_r)));
      }
      std::cout << GridLogMessage << "re-evaluated trace of the generator " << gen_id << " is " << sum << " " << sum_im << std::endl;
      assert ( sum < 1e-8) ;
      assert ( sum_im < 1e-8) ;
  }
}

static void symm_base() {

  const int this_nc = 6;
  const int this_n = this_nc/2;
  const int this_irrep_dim = Sp_TwoIndex<this_nc, Symmetric>::Dimension;
  const int this_algebra_dim = Sp<this_nc>::AlgebraDimension;
  typedef Sp_TwoIndex<this_nc, Symmetric>::iGroupMatrix<Complex> Matrix;
  typedef Sp_TwoIndex<this_nc, Symmetric>::iGroupTwoIndexMatrix<Complex> SMatrix;
  
  Matrix Omega;
  Matrix eij_a;
  Matrix eij_b;
  Matrix eij_c;
  Matrix e_sum;
  Omega = Zero();
  for (int i = 0; i < this_n; i++)
  {
    Omega()()(i, this_n + i) = 1.;
    Omega()()(this_n + i, i) = -1;
  }

  RealD realA;
  RealD realB;
  std::cout << GridLogMessage << "symm dimension is " << this_irrep_dim << std::endl;
  std::cout << GridLogMessage << "algebra dimension is " << this_algebra_dim << std::endl;
  realA = Sp_TwoIndex<this_nc, AntiSymmetric>::Dimension + Sp_TwoIndex<this_nc, Symmetric>::Dimension;
  assert ( realA == this_nc * this_nc - 1); // Nc x Nc = dim(2indxS) + dim(2indxAS) + dim(singlet)

  std::cout << GridLogMessage << "checking base is symmetric " << std::endl;
  for (int a=0; a < this_irrep_dim; a++)
  {
      Sp_TwoIndex<this_nc, Symmetric>::base(a, eij_c);
      e_sum = eij_c - transpose(eij_c);
      std::cout << GridLogMessage << "e_ab - e_ab^T " << norm2(e_sum) << std::endl;
      assert(norm2(e_sum) < 1e-8);
  }
  std::cout << GridLogMessage << "Checking Tr (e^(ab) Omega ) = 0 and Tr (e^(ab) e^(cd) = delta^((ab)(cd)) ) " << std::endl;
  for (int a=0; a < Sp_TwoIndex<this_nc, Symmetric>::Dimension; a++) {
      Sp_TwoIndex<this_nc, Symmetric>::base(a, eij_a);
      realA = norm2(trace(Omega*eij_a));
      std::cout << GridLogMessage << "Omega trace for (ab) = " << a << std::endl;
      assert(realA == 0);
      for (int b=0; b < Sp_TwoIndex<this_nc, Symmetric>::Dimension; b++) {
          Sp_TwoIndex<this_nc, Symmetric>::base(b, eij_b);
          auto d_ab = TensorRemove(trace(eij_a * eij_b));
#if verbose
          std::cout << GridLogMessage << "Tr( e_{ab=" << a << "} e_{cd=" << b << "} ) = " << d_ab << std::endl;
#endif
          std::cout << GridLogMessage << "Orthonormality for (ab) = " << a << std::endl;
          if (a==b) {
              assert(real(d_ab)-1 < 1e-8);
              assert(imag(d_ab) < 1e-8);
          } else {
              assert(real(d_ab) < 1e-8);
              assert(imag(d_ab) < 1e-8);
          }
      }
  }

  int sum = 0;
  int sum_im = 0;
  Vector<Matrix> ta_fund(this_algebra_dim);
  Vector<Matrix> eij(this_irrep_dim);
  Matrix tmp_l;
  Matrix tmp_r;
  for (int n = 0; n < this_algebra_dim; n++)
  {
      Sp<this_nc>::generator(n, ta_fund[n]);
  }
  for (int a = 0; a < this_irrep_dim; a++)
  {
      Sp_TwoIndex<this_nc, Symmetric>::base(a, eij[a]);
  }
  for (int gen_id = 0; gen_id < this_algebra_dim; gen_id++)
  {
      Complex iTr;
      sum = 0;
      sum_im = 0;
      std::cout << GridLogMessage <<  "generator number " << gen_id << std::endl;
      for (int a = 0; a < this_irrep_dim; a++)
      {
      
          tmp_l = adj(eij[a])*ta_fund[gen_id]*eij[a];
          tmp_r = adj(eij[a])*eij[a]*transpose(ta_fund[gen_id]);
#if verbose
          std::cout << GridLogMessage << " as_indx = " << a << " eDag T_F e = " << std::endl << tmp_l << std::endl;
          std::cout << GridLogMessage << " as_indx = " << a << " eDag e T_F^T = " << std::endl << tmp_r << std::endl;
#endif
          std::cout << GridLogMessage << " as_indx = " << a << " Tr(sum) = " << TensorRemove(trace(tmp_l+tmp_r)) << std::endl;
          sum += real(TensorRemove(trace(tmp_l+tmp_r)));
          sum_im += imag(TensorRemove(trace(tmp_l+tmp_r)));
      }
      std::cout << GridLogMessage << "re-evaluated trace of the generator " << gen_id << " is " << sum << " " << sum_im << std::endl;
      assert ( sum < 1e-8) ;
      assert ( sum_im < 1e-8) ;
  }
}


int main(int argc, char** argv) {
  std::cout << GridLogMessage << "Checking AntiSymmetric base " << std::endl;
  antisymm_base();
  std::cout << GridLogMessage << "*************** " << std::endl;
  std::cout << GridLogMessage << "Checking Symmetric base " << std::endl;
  symm_base();
}

