#include <Grid.h>

#include <qcd/utils/CovariantCshift.h>
#include <qcd/utils/WilsonLoops.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

class suN {
public:

  static int generators(int ncolour)   { return ncolour*ncolour-1; }
  static int su2subgroups(int ncolour) { return (ncolour*(ncolour-1))/2; }
  
  template<typename CComplex,int N> using suNmatrix = iScalar<iScalar<iMatrix<CComplex,N> > > ;

  ////////////////////////////////////////////////////////////////////////
  // There are N^2-1 generators for SU(N).
  //
  // We take a traceless hermitian generator basis as follows
  //
  // * Normalisation: trace ta tb = 1/2 delta_ab
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
  //    - for each row col pair there is a (sigma_x) and a (sigma_y) like generator
  //
  //
  //   t^a_ij = { in 0.. Nc(Nc-1)/2 -1} =>  delta_{i,i1} delta_{j,i2} +  delta_{i,i1} delta_{j,i2}  
  //   t^a_ij = { in Nc(Nc-1)/2 ... Nc^(Nc-1) -1} =>  i delta_{i,i1} delta_{j,i2} - i delta_{i,i1} delta_{j,i2}  
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
  ////////////////////////////////////////////////////////////////////////
  template<class CComplex,int Ncolour> 
  static void suNgenerator(int lieIndex,suNmatrix<CComplex,Ncolour> &ta){
    // map lie index to which type of generator
    int diagIndex;
    int su2Index;
    int sigxy;
    int NNm1 =  Ncolour*(Ncolour-1);
    if ( lieIndex>= NNm1 ) {
      diagIndex = lieIndex -NNm1;
      suNgeneratorDiagonal(diagIndex,ta);
      return;
    }
    sigxy   = lieIndex&0x1;
    su2Index= lieIndex>>1;
    if ( sigxy ) suNgeneratorSigmaY(su2Index,ta);
    else         suNgeneratorSigmaX(su2Index,ta);
  }
  template<class CComplex,int Ncolour> 
  static void suNgeneratorSigmaX(int su2Index,suNmatrix<CComplex,Ncolour> &ta){
    ta=zero;
    int i1,i2;
    su2SubGroupIndex<Ncolour>(i1,i2,su2Index);
    ta()()(i1,i2)=1.0;
    ta()()(i2,i1)=1.0;
    ta= ta *0.5;
  }
  template<class CComplex,int Ncolour> 
  static void suNgeneratorSigmaY(int su2Index,suNmatrix<CComplex,Ncolour> &ta){
    ta=zero;
    Complex i(0.0,1.0);
    int i1,i2;
    su2SubGroupIndex<Ncolour>(i1,i2,su2Index);
    ta()()(i1,i2)=-i;
    ta()()(i2,i1)= i;
    ta= ta *0.5;
  }
  template<class CComplex,int Ncolour> 
  static void suNgeneratorDiagonal(int diagIndex,suNmatrix<CComplex,Ncolour> &ta){
    ta=zero;
    int trsq=0;
    int last=diagIndex+1;
    for(int i=0;i<=diagIndex;i++){
      ta()()(i,i) = 1.0;
      trsq++;
    }
    ta()()(last,last) = -last;
    trsq+=last*last;
    RealD nrm = 1.0/std::sqrt(2.0*trsq);
    ta = ta *nrm;
  }
  ////////////////////////////////////////////////////////////////////////
  // Map a 
  //
  ////////////////////////////////////////////////////////////////////////
  template<int Ncolour>
  static void su2SubGroupIndex(int &i1,int &i2,int su2_index){

    assert( (su2_index>=0) && (su2_index< (Ncolour*(Ncolour-1))/2) );

    int spare=su2_index;
    for(i1=0;spare >= (Ncolour-1-i1);i1++ ){
      spare = spare - (Ncolour-1-i1);  // remove the Nc-1-i1 terms                                 
    }
    i2=i1+1+spare;
  }
  template<class CComplex,int Ncolour>
  static void su2Extract(std::vector<LatticeComplex> &r,const Lattice<suNmatrix<CComplex,Ncolour> > &source, int su2_index)
  {
    assert(r.size() == 4); // store in 4 real parts
    
    for(int i=0;i<4;i++){
      conformable(r[i],source);
    }
    
    int i1,i2;    
    su2SubGroupIndex<Ncolour>(i1,i2,su2_index);
    
    /* Compute the b(k) of A_SU(2) = b0 + i sum_k bk sigma_k */ 
    r[0] = real(source()()(i1,i1) + source()()(i2,i2));
    r[1] = imag(source()()(i1,i2) + source()()(i2,i1));
    r[2] = real(source()()(i1,i2) - source()()(i2,i1));
    r[3] = imag(source()()(i1,i1) - source()()(i2,i2));
  }


  template<int Ncolour> static void printGenerators(void)
  {
    for(int gen=0;gen<suN::generators(Ncolour);gen++){
      suN::suNmatrix<Complex,Ncolour> ta;
      suN::suNgenerator(gen,ta);
      std::cout<< "Nc = "<<Ncolour<<" t_"<<gen<<std::endl;
      std::cout<<ta<<std::endl;
    }
  }

  template<int Ncolour> static void testGenerators(void){
    suNmatrix<Complex,Ncolour> ta;
    suNmatrix<Complex,Ncolour> tb;
    std::cout<<"Checking trace ta tb is 0.5 delta_ab"<<std::endl;
    for(int a=0;a<generators(Ncolour);a++){
      for(int b=0;b<generators(Ncolour);b++){
	suNgenerator(a,ta);
	suNgenerator(b,tb);
	Complex tr =TensorRemove(trace(ta*tb)); 
	std::cout<<tr<<" ";
	if(a==b) assert(abs(tr-Complex(0.5))<1.0e-6);
	if(a!=b) assert(abs(tr)<1.0e-6);
      }
      std::cout<<std::endl;
    }
    std::cout<<"Checking hermitian"<<std::endl;
    for(int a=0;a<generators(Ncolour);a++){
      suNgenerator(a,ta);
      std::cout<<a<<" ";
      assert(norm2(ta-adj(ta))<1.0e-6);
    }    
    std::cout<<std::endl;

    std::cout<<"Checking traceless"<<std::endl;
    for(int a=0;a<generators(Ncolour);a++){
      suNgenerator(a,ta);
      Complex tr =TensorRemove(trace(ta)); 
      std::cout<<a<<" ";
      assert(abs(tr)<1.0e-6);
    }    
    std::cout<<std::endl;
  }

};


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);


  std::vector<int> simd_layout = GridDefaultSimd(4,vComplexF::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  std::vector<int> latt_size  ({4,4,4,4});

  GridCartesian     Fine(latt_size,simd_layout,mpi_layout);

  LatticeGaugeField Umu(&Fine);

  std::cout<<"*********************************************"<<std::endl;
  std::cout<<"* Generators for SU(2)"<<std::endl;
  std::cout<<"*********************************************"<<std::endl;
  suN::printGenerators<2>();
  suN::testGenerators<2>();
  std::cout<<"*********************************************"<<std::endl;
  std::cout<<"* Generators for SU(3)"<<std::endl;
  std::cout<<"*********************************************"<<std::endl;
  suN::printGenerators<3>();
  suN::testGenerators<3>();
  std::cout<<"*********************************************"<<std::endl;
  std::cout<<"* Generators for SU(4)"<<std::endl;
  std::cout<<"*********************************************"<<std::endl;
  suN::printGenerators<4>();
  suN::testGenerators<4>();
  std::cout<<"*********************************************"<<std::endl;
  std::cout<<"* Generators for SU(5)"<<std::endl;
  std::cout<<"*********************************************"<<std::endl;
  suN::printGenerators<5>();
  suN::testGenerators<5>();

  Grid_finalize();
}


