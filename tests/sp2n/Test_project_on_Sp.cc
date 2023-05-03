#include <Grid/Grid.h>

using namespace Grid;

template <typename T>
bool has_correct_group_block_structure(const T& U){
    std::cout << GridLogMessage << "Checking the structure is " << std::endl;
    std::cout << GridLogMessage << "U  =  (   W    X   )  " << std::endl;
    std::cout << GridLogMessage << "      (  -X^*  W^* )  " << std::endl;
    std::cout <<GridLogMessage << std::endl;

    const int nsp = Nc / 2;
    Complex i(0., 1.);
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
    return true;
};

template <typename T>
bool is_element_of_sp2n_group(T U) {// does explicitly take a copy in order to not spoil the matrix for further use
 
    LatticeColourMatrixD aux(U.Grid());
    LatticeColourMatrixD identity(U.Grid());
    identity = 1.0;
    
    std::cout << GridLogMessage << "Unitary check" << std::endl;
    aux = U*adj(U) - identity;
    std::cout << GridLogMessage << "U adjU - 1 = " << norm2(aux) << std::endl;
    assert( norm2(aux) < 1e-8);
    
    std::cout << GridLogMessage << "Checking Omega invariance" << std::endl;
    Sp<Nc>::OmegaInvariance(U);     // no assertion here, but the next check will kill us if we are not simplectic
    
  return has_correct_group_block_structure(U);
}

template<typename T>
void test_group_projections(T U) {
    RealD Delta = 666.;
    LatticeColourMatrixD identity(U.Grid());
    identity = 1.0;

    std::cout <<GridLogMessage << std::endl;
    std::cout << GridLogMessage << "#   #   #   #" << std::endl;
    std::cout << GridLogMessage << "Group" << std::endl;
    std::cout << GridLogMessage << "#   #   #   #" << std::endl;
    std::cout <<GridLogMessage << std::endl;

    std::string name="ProjectOnSpGroup";
    std::cout << GridLogMessage << "Testing "<< name << std::endl;
    std::cout << GridLogMessage << "Apply to deformed matrix" << std::endl;

    U = U + Delta*identity;
    U = ProjectOnSpGroup(U);
    assert(is_element_of_sp2n_group(U));
   
    name = "ProjectOnGaugeGroup";
    std::cout << GridLogMessage << "Testing "<< name << std::endl;
    std::cout << GridLogMessage << "Apply to deformed matrix" << std::endl;

    U = U + Delta*identity;
    Sp<Nc>::ProjectOnGaugeGroup(U);
    assert(is_element_of_sp2n_group(U));
   
    name = "ProjectGn";
    std::cout << GridLogMessage << "Testing "<< name << std::endl;
    std::cout << GridLogMessage << "Apply to deformed matrix" << std::endl;

    U = U + Delta*identity;
    Sp<Nc>::ProjectGn(U);
    assert(is_element_of_sp2n_group(U));
}



int main (int argc, char **argv)
{
    Grid_init(&argc,&argv);
    
    
    Coordinate latt_size   = GridDefaultLatt();
    Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();

    GridCartesian             Grid(latt_size,simd_layout,mpi_layout);
    
    
    LatticeGaugeField Umu(&Grid);
    LatticeColourMatrixD U(&Grid);
    LatticeColourMatrixD Up(&Grid);
    LatticeColourMatrixD aux(&Grid);
    LatticeColourMatrixD identity(&Grid);
    
    // Will test resimplectification-related functionalities (from ProjectOnGaugeGroup, ProjectOnSpGroup, ProjectGn) and projection on the algebra (from ProjectSp2nAlgebra)
    // we work with matrices with positive determinant so detU = 1 even if in principle ProjectOnGaugeGroup and ProjectOnSpGroup allow for detU=-1
    // so the checks will be the same for the three functions
    // NB only ProjectGn is the proper simplectification function
    
    
    const int nsp = Nc / 2;
    
    identity = 1.0;
    RealD epsilon = 0.01;
    RealD Delta = 666.;
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
    
    test_group_projections(U);

    std::cout <<GridLogMessage << std::endl;
    std::cout << GridLogMessage << "#   #   #   #" << std::endl;
    std::cout << GridLogMessage << "Algebra" << std::endl;
    std::cout << GridLogMessage << "#   #   #   #" << std::endl;
    std::cout <<GridLogMessage << std::endl;
    std::cout << GridLogMessage << "Testing SpTa" << std::endl;
    
    U = PeekIndex<LorentzIndex>(Umu,1);
    U = U + Delta*identity;
    std::cout << GridLogMessage << "Matrix deformed " << std::endl;
    std::cout << GridLogMessage << "Apply SpTa to deformed matrix" << std::endl;
    U = SpTa(U);
    
    aux = U - adj(U);
    std::cout << GridLogMessage << "SpTa ::: T - Tda = " << norm2(aux) << std::endl;
    aux = U + adj(U);
    std::cout << GridLogMessage << "SpTa ::: T + Tda = " << norm2(aux) << std::endl;
    
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
    
    //test Ta
    /*
    U = U + 666.*identity;
    Up = Ta(U);
    aux = Up - adj(Up);
    std::cout << GridLogMessage << "TA !!! T - Tda = " << norm2(aux) << std::endl;
    aux = Up + adj(Up);
    std::cout << GridLogMessage << "TA !!! T + Tda = " << norm2(aux) << std::endl;*/
    
    // test taProj
    std::cout << GridLogMessage << "Testing taProj" << std::endl;
    U = U + Delta*identity;
    std::cout << GridLogMessage << "Matrix deformed " << std::endl;
    std::cout << GridLogMessage << "Apply taProj to deformed matrix" << std::endl;
    Sp<Nc>::taProj(U, Up);
    aux = Up - adj(Up);
    std::cout << GridLogMessage << "taProj ::: T - Tda = " << norm2(aux) << std::endl;
    aux = Up + adj(Up);
    std::cout << GridLogMessage << "taProj ::: T + Tda = " << norm2(aux) << std::endl;
    
    std::cout << GridLogMessage << "Check that Omega U Omega = conj(U)" << std::endl;
    Sp<Nc>::Omega(Omega);
    aux = Omega*Up*Omega - conjugate(Up);
    std::cout << GridLogMessage << "Omega U Omega - conj(U) = " << norm2(aux) << std::endl;
    assert( norm2(aux) < 1e-8);
    // before it was
    aux = Omega*U*Omega - conjugate(U);
    std::cout << GridLogMessage << " before taProj Omega U Omega - conj(U) = " << norm2(aux) << std::endl;

    U = Up;
    
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
