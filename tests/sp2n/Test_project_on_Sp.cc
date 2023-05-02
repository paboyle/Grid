#include <Grid/Grid.h>

using namespace Grid;

int SpGroupQuiz (const LatticeColourMatrixD Uin)
{
    double vol = Uin.Grid()->gSites();
    const int nsp = Nc / 2;
    LatticeColourMatrixD Omega(Uin.Grid());
    Sp<Nc>::Omega(Omega);
    LatticeColourMatrixD aux(Uin.Grid());
    LatticeColourMatrixD identity(Uin.Grid());
    Complex i(0., 1.);
    identity = 1;
    
    std::cout << GridLogMessage << "Check matrix is non-zero " << std::endl;
    assert(norm2(Uin) > 1e-8);
    aux = Uin*adj(Uin) - identity;
    std::cout << GridLogMessage << "U adjU - 1 = " << norm2(aux)/vol << std::endl;
    assert(norm2(aux) < 1e-8);
    aux = Omega - (Uin * Omega * transpose(Uin));
    std::cout << GridLogMessage << "Omega - U Omega transpose(U) = " << norm2(aux)/vol << std::endl;
    assert(norm2(aux) < 1e-8);
    //Sp<Nc>::OmegaInvariance(Uin)
    std::cout << GridLogMessage << "Checking the structure is " << std::endl;
    std::cout << GridLogMessage << "U  =  (   W    X   )  " << std::endl;
    std::cout << GridLogMessage << "      (  -X^*  W^* )  " << std::endl;
    for (int c1 = 0; c1 < nsp; c1++) //check on W
    {
        for (int c2 = 0; c2 < nsp; c2++)
        {
            auto W = PeekIndex<ColourIndex>(Uin,c1,c2);
            auto Wstar =  PeekIndex<ColourIndex>(Uin,c1+nsp,c2+nsp);
            auto Ww = conjugate( Wstar );
            auto amizero = sum(W - Ww);
            auto amizeroo = TensorRemove(amizero);
            assert(  amizeroo.real() < 10e-6 );
            amizeroo *= i;
            assert(  amizeroo.real() < 10e-6 );
        }
    }
    for (int c1 = 0; c1 < nsp ; c1++)
    {
        for (int c2 = 0; c2 < nsp; c2++)
        {
            auto X = PeekIndex<ColourIndex>(Uin,c1,c2+nsp);
            auto minusXstar = PeekIndex<ColourIndex>(Uin,c1+nsp,c2);
            auto minusXx = conjugate(minusXstar);
            auto amizero = sum (X + minusXx);
            auto amizeroo = TensorRemove(amizero);
            assert(  amizeroo.real() < 10e-6 );
            amizeroo *= i;
            assert(  amizeroo.real() < 10e-6 );
        }
    }
    std::cout << GridLogMessage << "|Det| = " << norm2( Determinant(Uin) ) / vol << std::endl;
    assert( norm2( Determinant(Uin) ) / vol - 1 < 1e-8);
    return 0;
}

int SpAntiHermitianAlgebraQuiz (const LatticeColourMatrixD Uin)
{
    double vol = Uin.Grid()->gSites();
    const int nsp = Nc / 2;
    LatticeColourMatrixD Omega(Uin.Grid());
    Sp<Nc>::Omega(Omega);
    LatticeColourMatrixD aux(Uin.Grid());
    LatticeColourMatrixD identity(Uin.Grid());
    Complex i(0., 1.);
    identity = 1;
    
    std::cout << GridLogMessage << "Check matrix is non-zero " << std::endl;
    assert(norm2(Uin) > 1e-8);
    aux = Uin - adj(Uin);
    std::cout << GridLogMessage << "SpTa ::: T - Tda = " << norm2(aux)/vol << std::endl;
    aux = Uin + adj(Uin);
    std::cout << GridLogMessage << "SpTa ::: T + Tda = " << norm2(aux)/vol << std::endl;
    assert( norm2(aux) - 1 < 1e-8);
    std::cout << GridLogMessage << "Check that Omega T Omega + conj(T) = 0 " << std::endl;
    aux = Omega*Uin*Omega + conjugate(Uin);
    assert( norm2(aux) < 1e-8);
    std::cout << GridLogMessage << "Checking the structure is " << std::endl;
    std::cout << GridLogMessage << "U  =  (   W    X   )  " << std::endl;
    std::cout << GridLogMessage << "      (  -X^*  W^* )  " << std::endl;
    for (int c1 = 0; c1 < nsp; c1++) //check on W
    {
        for (int c2 = 0; c2 < nsp; c2++)
        {
            auto W = PeekIndex<ColourIndex>(Uin,c1,c2);
            auto Wstar =  PeekIndex<ColourIndex>(Uin,c1+nsp,c2+nsp);
            auto Ww = conjugate( Wstar );
            auto amizero = sum(W - Ww);
            auto amizeroo = TensorRemove(amizero);
            assert(  amizeroo.real() < 10e-6 );
            amizeroo *= i;
            assert(  amizeroo.real() < 10e-6 );
        }
    }
    for (int c1 = 0; c1 < nsp ; c1++)
    {
        for (int c2 = 0; c2 < nsp; c2++)
        {
            auto X = PeekIndex<ColourIndex>(Uin,c1,c2+nsp);
            auto minusXstar = PeekIndex<ColourIndex>(Uin,c1+nsp,c2);
            auto minusXx = conjugate(minusXstar);
            auto amizero = sum (X + minusXx);
            auto amizeroo = TensorRemove(amizero);
            assert(  amizeroo.real() < 10e-6 );
            amizeroo *= i;
            assert(  amizeroo.real() < 10e-6 );
        }
    }
    return 0;
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
    LatticeColourMatrixD aux(&Grid);
    LatticeColourMatrixD identity(&Grid);
    LatticeColourMatrixD Omega(&Grid);
    Sp<Nc>::Omega(Omega);

    identity = 1.0;
    RealD epsilon = 0.1;
    Complex i(0., 1.);
    double vol = Umu.Grid()->gSites();
    const int nsp = Nc / 2;
    
    std::vector<int> pseeds({1,2,3,4,5});
    GridParallelRNG  pRNG(&Grid); pRNG.SeedFixedIntegers(pseeds);
    
    SU<Nc>::HotConfiguration(pRNG,Umu);
    U = PeekIndex<LorentzIndex>(Umu,0);
    
    std::cout << GridLogMessage << "Starting with random SUn matrix " << std::endl;
    aux = U*adj(U) - identity;
    std::cout <<GridLogMessage << "U adjU - 1 = " << norm2(aux) << std::endl;
    assert ( norm2(aux) < 1e-8 );
    if (Nc != 2)    // Sp2 = SU2
    {
        std::cout << GridLogMessage <<  "Checking matrix does NOT leave Omega invariant, to avoid a trivial test " << std::endl;
        aux = Omega - (U * Omega * transpose(U));
        assert ( norm2(aux) > 1e-8 );
    }
    
    U = U + epsilon*identity - i*identity;
    std::cout << GridLogMessage << "Unitary matrix deformed " << std::endl;
    auto det =  sum( Determinant(U) );
    std::cout << GridLogMessage << "Re(Det) = " << real(det) / vol << std::endl;
    det = det*i;
    std::cout << GridLogMessage << "Im(Det) = " << real(det) / vol << std::endl;
    aux = U*adj(U) - identity;
    std::cout << GridLogMessage << "now U adjU - 1 = " << norm2(aux) << std::endl;
    
    //  ProjectOnSpGroup
    std::cout << GridLogMessage << "Testing ProjectOnSpGroup" << std::endl;
    std::cout << GridLogMessage << "Apply ProjectOnSpGroup to deformed matrix" << std::endl;
    U = ProjectOnSpGroup(U);
    std::cout << GridLogMessage << "Run checks:" << std::endl;
    SpGroupQuiz(U);
    det =  sum( Determinant(U) );
    std::cout << GridLogMessage << "Re(Det) after ProjectOnSpGroup (nothing to assert) = " << real(det) / vol << std::endl;
    det = det*i;
    std::cout << GridLogMessage << "Im(Det) after ProjectOnSpGroup (nothing to assert) = " << real(det) / vol << std::endl;
    
    //  ProjectOnGaugeGroup
    SU<Nc>::HotConfiguration(pRNG,Umu); //refresh
    U = PeekIndex<LorentzIndex>(Umu,0);
    std::cout << GridLogMessage << "Testing ProjectOnGaugeGroup" << std::endl;
    U = U + epsilon*identity - i*identity;
    std::cout << GridLogMessage << "Apply ProjectOnGaugeGroup to deformed matrix" << std::endl;
    Sp<Nc>::ProjectOnGaugeGroup(U);
    std::cout << GridLogMessage <<"Run checks:" << std::endl;
    SpGroupQuiz(U);
    det =  sum( Determinant(U) );
    std::cout << GridLogMessage << "Re(Det) after ProjectOnGaugeGroup (nothing to assert) = " << real(det) / vol << std::endl;
    det = det*i;
    std::cout << GridLogMessage << "Im(Det) after ProjectOnGaugeGroup (nothing to assert) = " << real(det) / vol << std::endl;
    
    //  ProjectGn
    SU<Nc>::HotConfiguration(pRNG,Umu); //refresh
    U = PeekIndex<LorentzIndex>(Umu,0);
    std::cout << GridLogMessage << "Testing ProjectGn" << std::endl;
    U = U + epsilon*identity - i*identity;
    std::cout << GridLogMessage << "Apply ProjectGn to deformed matrix" << std::endl;
    Sp<Nc>::ProjectGn(U);
    std::cout << GridLogMessage << "Run checks:" << std::endl;
    SpGroupQuiz(U);
    det =  sum( Determinant(U) );
    std::cout << GridLogMessage << "Re(Det) after ProjectGn (constrained) = " << real(det) / vol << std::endl;
    assert ( (real(det) / vol) - 1 < 1e-8 );
    det = det*i;
    std::cout << GridLogMessage << "Im(Det) after ProjectGn (constrained) = " << real(det) / vol << std::endl;
    assert ( real(det) / vol < 1e-8 );
    
    //  SpTa
    SU<Nc>::HotConfiguration(pRNG,Umu);//refresh
    U = PeekIndex<LorentzIndex>(Umu,0);
    std::cout << GridLogMessage << "Testing SpTa" << std::endl;
    U = U + epsilon*identity - i*identity;
    std::cout << GridLogMessage << "Apply SpTa to deformed matrix" << std::endl;
    U = SpTa(U);
    std::cout << GridLogMessage << "Run checks:" << std::endl;
    SpAntiHermitianAlgebraQuiz(U);
        
    Grid_finalize();

}
