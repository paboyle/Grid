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
    
    // Will test resimplectification (from ProjectOnGaugeGroup, ProjectOnSpGroup) and projection on the algebra (from ProjectSp2nAlgebra)
    
    const int nsp = Nc / 2;
    
    identity = 1.0;
    RealD epsilon = 0.01;
    Complex i(0., 1.);
    RealD u = 0.;
    double vol = Umu.Grid()->gSites();
    
    std::vector<int> pseeds({1,2,3,4,5});
    GridParallelRNG  pRNG(&Grid); pRNG.SeedFixedIntegers(pseeds);
    
    SU<Nc>::HotConfiguration(pRNG,Umu);
    U = PeekIndex<LorentzIndex>(Umu,0);
    
    aux = U*adj(U) - identity;
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
    std::cout << GridLogMessage << "Simplectify" << std::endl;

    U = ProjectOnSpGroup(U);

    aux = U*adj(U) - identity;
    std::cout <<GridLogMessage << std::endl;
    std::cout <<GridLogMessage << std::endl;
    std::cout << GridLogMessage << "Unitary check after simplectification " << std::endl;
    std::cout << GridLogMessage << "U adjU - 1 = " << norm2(aux) << std::endl;
    assert( norm2(aux) < 1e-8);
    
    // checks on determinant
    std::cout << GridLogMessage << "Det after simplectification = " << norm2( Determinant(U) ) / vol << std::endl;
    assert( norm2(aux) - 1 < 1e-8);
    
    // actual sp2n check
    std::cout << GridLogMessage << "Checking invariance after simplectification. Will not kill if it fails, but the next check will " << std::endl;
    Sp<Nc>::OmegaInvariance(U);
    
    // checks on elements
    
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
    
    std::cout << GridLogMessage << "Check the ProjectOnGaugeGroup function to simplectify" << std::endl;
    U = U + 2932.111*identity;
    std::cout << GridLogMessage << "Resimplectify deformed matrix" << std::endl;
    Sp<Nc>::ProjectOnGaugeGroup(U);
    aux = U*adj(U) - identity;
    std::cout << GridLogMessage << "U adjU - 1 = " << norm2(aux) << std::endl;
    assert( norm2(aux) < 1e-8);
    Sp<Nc>::OmegaInvariance(U);
    
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
    
    std::cout << GridLogMessage << "Now check projection on the algebra" << std::endl;
    
    U = PeekIndex<LorentzIndex>(Umu,1);
    U = U + 666.*identity;
    aux = U*adj(U) - identity;
    std::cout << GridLogMessage << "Matrix deformed " << std::endl;
    std::cout << GridLogMessage << "U adjU - 1 = " << norm2(aux) << std::endl;
    std::cout << GridLogMessage << "Project on sp2n algebra" << std::endl;
    U = ProjectSp2nAlgebra(U);
    aux = U*adj(U) - identity;
    std::cout << GridLogMessage << "U adjU - 1 = " << norm2(aux) << std::endl;
    assert( norm2(aux) < 1e-8);
    std::cout << GridLogMessage << "Check that Omega U Omega = conj(U)" << std::endl;
    
    LatticeColourMatrixD Omega(&Grid);
    
    Sp<Nc>::Omega(Omega);
    aux = Omega*U*Omega - conjugate(U);
    std::cout << GridLogMessage << "Omega U Omega - conj(U) = " << norm2(aux) << std::endl;
    assert( norm2(aux) < 1e-8);
    
    
    std::cout << GridLogMessage << "Checking the structure is " << std::endl;
    std::cout << GridLogMessage << "U  =  (   W    X   )  " << std::endl;
    std::cout << GridLogMessage << "      (  X^*  -W^* )  " << std::endl;
    std::cout <<GridLogMessage << std::endl;
    for (int c1 = 0; c1 < nsp; c1++) //check on W
    {
        for (int c2 = 0; c2 < nsp; c2++)
        {
            auto W = PeekIndex<ColourIndex>(U,c1,c2);
            auto Wstar =  PeekIndex<ColourIndex>(U,c1+nsp,c2+nsp);
            auto Ww = conjugate( Wstar );
            auto amizero = sum(W + Ww);
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
            auto amizero = sum (X - minusXx);
            auto amizeroo = TensorRemove(amizero);
            assert(  amizeroo.real() < 10e-6 );
            amizeroo *= i;
            assert(  amizeroo.real() < 10e-6 );
        }
    }
    
    
    
    Grid_finalize();


}
