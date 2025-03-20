#include <Grid/Grid.h>

#define verbose 0

using namespace Grid;

template<int this_nc>
static void check_dimensions() {
    
    const int this_n = this_nc/2;
    const int this_algebra_dim = Sp<this_nc>::AlgebraDimension;
    
    RealD realA;
    std::cout << GridLogMessage << "Nc = " << this_n << " 2as dimension is " << Sp_TwoIndex<this_nc, AntiSymmetric>::Dimension << std::endl;
    std::cout << GridLogMessage << "Nc = " << this_n << " 2s dimension is " << Sp_TwoIndex<this_nc, Symmetric>::Dimension << std::endl;
    std::cout << GridLogMessage << "Nc = " << this_n << " algebra dimension is " << this_algebra_dim << std::endl;
    realA = Sp_TwoIndex<this_nc, AntiSymmetric>::Dimension + Sp_TwoIndex<this_nc, Symmetric>::Dimension;
    std::cout << GridLogMessage << "Checking dim(2AS) + dim(AS) + 1 = Nc * Nc " << this_algebra_dim << std::endl;
    assert ( realA == this_nc * this_nc - 1); // Nc x Nc = dim(2indxS) + dim(2indxAS) + dim(singlet)
}

template<int this_nc, TwoIndexSymmetry S>
static void run_symmetry_checks() {
    typedef typename Sp_TwoIndex<this_nc, S>::template iGroupMatrix<Complex> Matrix;
    const int this_n = this_nc/2;
    const int this_irrep_dim = Sp_TwoIndex<this_nc, S>::Dimension;
    const int this_algebra_dim = Sp<this_nc>::AlgebraDimension;
    Matrix eij_c;
    Matrix e_sum;
    RealD realS = S;
    
    std::cout << GridLogMessage << "checking base has symmetry " << S << std::endl;
    for (int a=0; a < this_irrep_dim; a++)
    {
        Sp_TwoIndex<this_nc, S>::base(a, eij_c);
        e_sum = eij_c - realS * transpose(eij_c);
        std::cout << GridLogMessage << "e_ab - (" << S << " * e_ab^T ) = " << norm2(e_sum) << std::endl;
        assert(norm2(e_sum) < 1e-8);
          
    }
}

template<int this_nc, TwoIndexSymmetry S>
static void run_traces_checks() {
    typedef typename Sp_TwoIndex<this_nc, S>::template iGroupMatrix<Complex> Matrix;
    const int this_n = this_nc/2;
    const int this_irrep_dim = Sp_TwoIndex<this_nc, S>::Dimension;
    const int this_algebra_dim = Sp<this_nc>::AlgebraDimension;
    Matrix eij_a;
    Matrix eij_b;
    Matrix Omega;
    Sp<this_nc>::Omega(Omega);
    RealD realS = S;
    RealD realA;
    
    std::cout << GridLogMessage << "Checking Tr (e^(ab) Omega ) = 0 and Tr (e^(ab) e^(cd) = delta^((ab)(cd)) ) " << std::endl;
    for (int a=0; a < Sp_TwoIndex<this_nc, S>::Dimension; a++) {
        Sp_TwoIndex<this_nc, S>::base(a, eij_a);
        realA = norm2(trace(Omega*eij_a));
        std::cout << GridLogMessage << "Checkig Omega-trace for e_{ab=" << a << "} " << std::endl;
        //std::cout << GridLogMessage << "Tr ( Omega e_{ab=" << a << "} ) = " << realA << std::endl;
        assert(realA < 1e-8);
        for (int b=0; b < Sp_TwoIndex<this_nc, S>::Dimension; b++) {
            Sp_TwoIndex<this_nc, S>::base(b, eij_b);
            auto d_ab = TensorRemove(trace(eij_a * eij_b));
    #if verbose
            std::cout << GridLogMessage << "Tr( e_{ab=" << a << "} e_{cd=" << b << "} ) = " << d_ab << std::endl;
    #endif
            std::cout << GridLogMessage << "Checking orthonormality for e_{ab = " << a << "} " << std::endl;
            if (a==b) {
                assert(real(d_ab) - realS < 1e-8);
            } else {
                assert(real(d_ab) < 1e-8);
            }
            assert(imag(d_ab) < 1e-8);
            assert(imag(d_ab) < 1e-8);
        }
    }
    
}

template<int this_nc, TwoIndexSymmetry S>
static void run_generators_checks() {
    const int this_n = this_nc/2;
    const int this_irrep_dim = Sp_TwoIndex<this_nc, S>::Dimension;
    const int this_algebra_dim = Sp<this_nc>::AlgebraDimension;
    typedef typename Sp_TwoIndex<this_nc, S>::template iGroupMatrix<Complex> Matrix;
    int sum = 0;
    int sum_im = 0;
    std::vector<Matrix> ta_fund(this_algebra_dim);
    std::vector<Matrix> eij(this_irrep_dim);
    Matrix tmp_l;
    Matrix tmp_r;
    for (int n = 0; n < this_algebra_dim; n++)
    {
        Sp<this_nc>::generator(n, ta_fund[n]);  // generators in the fundamental
    }
      for (int a = 0; a < this_irrep_dim; a++)
    {
        Sp_TwoIndex<this_nc, S>::base(a, eij[a]);   // base functions e_ij^a for upgrading gauge links from fund to 2-index
    }
    for (int gen_id = 0; gen_id < this_algebra_dim; gen_id++)
    {
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
            //std::cout << GridLogMessage << " as_indx = " << a << " Tr(eDag T_F e + eDag e T_F^T) = " << TensorRemove(trace(tmp_l+tmp_r)) << std::endl;
            sum += real(TensorRemove(trace(tmp_l+tmp_r)));
            sum_im += imag(TensorRemove(trace(tmp_l+tmp_r)));
        }
        std::cout << GridLogMessage << "re-evaluated trace of the generator " << gen_id << " is " << sum << " " << sum_im << std::endl;
        assert ( sum < 1e-8) ;
        assert ( sum_im < 1e-8) ;
    }
    
}
    
template<int this_nc, TwoIndexSymmetry S>
static void run_base_checks() {
    std::cout << GridLogMessage << " ****** " << std::endl;
    std::cout << GridLogMessage << "Running checks for Nc = " << this_nc << " TwoIndex Symmetry = " << S << std::endl;
    run_symmetry_checks<this_nc, S>();
    run_traces_checks<this_nc, S>();
    run_generators_checks<this_nc, S>();
}

int main(int argc, char** argv) {
    check_dimensions<2>();
    check_dimensions<4>();
    check_dimensions<6>();
    check_dimensions<8>();
    
    run_base_checks<2, Symmetric>();    // For Nc=2 the AS is the singlet
    run_base_checks<4, Symmetric>();
    run_base_checks<4, AntiSymmetric>();
    run_base_checks<6, Symmetric>();
    run_base_checks<6, AntiSymmetric>();
    run_base_checks<8, Symmetric>();
    run_base_checks<8, AntiSymmetric>();
}
