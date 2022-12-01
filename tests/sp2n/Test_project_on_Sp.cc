#include <Grid/Grid.h>

using namespace Grid;

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
    
    const int nsp = Nc / 2;
    
    identity = 1.0;
    RealD epsilon = 0.01;
    Complex i(0., 1.);
    RealD u = 0.;
    double vol = Umu.Grid()->gSites();
    
    std::vector<int> pseeds({1,2,3,4,5});
    std::vector<int> sseeds({6,7,8,9,10});
    GridParallelRNG  pRNG(&Grid); pRNG.SeedFixedIntegers(pseeds);
    GridSerialRNG    sRNG;       sRNG.SeedFixedIntegers(sseeds);
    
    SU<Nc>::HotConfiguration(pRNG,Umu);
    U = PeekIndex<LorentzIndex>(Umu,2);
    
    aux = U*adj(U) - identity;
    std::cout <<GridLogMessage << std::endl;
    std::cout << GridLogMessage << "Starting with random SUn matrix " << std::endl;
    std::cout << GridLogMessage << "Unitary check " << std::endl;
    std::cout <<GridLogMessage << "U adjU - 1 = " << norm2(aux) << std::endl;
    assert ( norm2(aux) < 1e-8 );
    std::cout <<GridLogMessage << std::endl;
    if (Nc != 2)
    {
        std::cout << GridLogMessage << "This matrix should not leave Omega invariant, expect a warning" << std::endl;
    }
    Sp<Nc>::OmegaInvariance(U);
    std::cout <<GridLogMessage << std::endl;
    
    U = U + epsilon*identity;
    aux = U*adj(U) - identity;
    
    std::cout << GridLogMessage << "Unitary matrix deformed " << std::endl;
    std::cout << GridLogMessage << "now U adjU - 1 = " << norm2(aux) << std::endl;
    std::cout <<GridLogMessage << std::endl;

    std::cout <<GridLogMessage << std::endl;
    std::cout << GridLogMessage << "Projecting on Sp2n " << std::endl;

    U = ProjectOnSpGroup(U);
    //U = ProjectOnGroup(U);
    aux = U*adj(U) - identity;
    std::cout <<GridLogMessage << std::endl;
    std::cout <<GridLogMessage << std::endl;
    std::cout << GridLogMessage << "Unitary check after Sp(2n) projection " << std::endl;
    std::cout << GridLogMessage << "U adjU - 1 = " << norm2(aux) << std::endl;
    assert( norm2(aux) < 1e-8);
    std::cout <<GridLogMessage << std::endl;
    std::cout <<GridLogMessage << std::endl;
    
    // checks on determinant
    std::cout << GridLogMessage << "Det after Projection on Sp2n = " << norm2( Determinant(U) ) / vol << std::endl;
    std::cout <<GridLogMessage << std::endl;
    std::cout <<GridLogMessage << std::endl;
    
    // actual sp2n check
    std::cout << GridLogMessage << "Checking invariance after projection "<< std::endl;
    Sp<Nc>::OmegaInvariance(U);
    
    // checks on elements
    
    std::cout <<GridLogMessage << std::endl;
    std::cout <<GridLogMessage << std::endl;
    std::cout << GridLogMessage << "Checking the structure is " << std::endl;
    std::cout << GridLogMessage << "U  =  (   W    X   )  " << std::endl;
    std::cout << GridLogMessage << "      (  -X^*  W^* )  " << std::endl;
    std::cout <<GridLogMessage << std::endl;
    for (int c1 = 0; c1 < nsp; c1++) //check on W
    {
        for (int c2 = 0; c2 < nsp; c2++)
        {
            auto W = PeekIndex<ColourIndex>(U,c1,c2);
            auto Wstar =  PeekIndex<ColourIndex>(U,c1+nsp,c2+nsp);
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
            auto X = PeekIndex<ColourIndex>(U,c1,c2+nsp);
            auto minusXstar = PeekIndex<ColourIndex>(U,c1+nsp,c2);
            auto minusXx = conjugate(minusXstar);
            auto amizero = sum (X + minusXx);
            auto amizeroo = TensorRemove(amizero);
            assert(  amizeroo.real() < 10e-6 );
            amizeroo *= i;
            assert(  amizeroo.real() < 10e-6 );
        }
    }
    
    std::cout << GridLogMessage << "ok" << std::endl;
    
    // an explicit check for sp2
    /*
    if (Nc == 2)
    {
        assert(Nc==2);
        ColourMatrix A;
        A = Zero();

        Complex a(25041994., 12.);
        Complex b(39., 0.22);
        Complex d(10000., -2222.3333);
    
        A()()(0,0) = a;
        A()()(0,1) = b;
        A()()(1,0) = i;
        A()()(1,1) = d;
        std::cout <<GridLogMessage << std::endl;
        std::cout <<GridLogMessage << std::endl;
        std::cout << GridLogMessage << "An explicit check for Sp2" << std::endl;
        std::cout <<GridLogMessage << std::endl;
        std::cout << GridLogMessage << "Building a non unitary matrix by hand with funny entries " << std::endl;
        std::cout << GridLogMessage << "A = " << A << std::endl;
        std::cout << GridLogMessage << "Projecting on Sp2 " << std::endl;
        A = ProjectOnSpGroup(A);
        std::cout << GridLogMessage << "now A = " << A << std::endl;
        std::cout << GridLogMessage << "A(0,0) - conjA(1,1) = " << A()()(0,0) - adj ( A()()(1,1) )<< std::endl;
        std::cout << GridLogMessage << "A(0,1) + conjA(1,0) = " << A()()(0,1) + adj ( A()()(1,0) )<< std::endl;
    }*/
    
    Grid_finalize();


}
