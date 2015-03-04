#include "Grid.h"
#include "Grid_vRealD.h"
#include "Grid_vRealF.h"
#include "Grid_vComplexD.h"
#include "Grid_vComplexF.h"
#include "Grid_Cartesian.h"
#include "Grid_Lattice.h"

using namespace std;
using namespace dpo;
using namespace dpo::QCD;


int main (int argc, char ** argv)
{

#ifdef KNL
  struct sigaction sa,osa;
  sigemptyset (&sa.sa_mask);
  sa.sa_sigaction= sa_action;
  sa.sa_flags    = SA_ONSTACK|SA_SIGINFO;

  sigaction(SIGSEGV,&sa,NULL);
  sigaction(SIGTRAP,&sa,NULL);
#endif

    std::vector<int> latt_size(4);
    std::vector<int> simd_layout(4);

 for(int omp=32;omp<237;omp*=2){

#ifdef OMP
   omp_set_num_threads(omp);
#endif 

  for(int lat=16;lat<=32;lat+=8){
    latt_size[0] = lat;
    latt_size[1] = lat;
    latt_size[2] = lat;
    latt_size[3] = lat;
    double volume = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
    
    simd_layout[0] = 1;
    simd_layout[1] = 2;
    simd_layout[2] = 2;
    simd_layout[3] = 2;
    
    GridCartesian Fine(latt_size,simd_layout);
    GridRedBlackCartesian rbFine(latt_size,simd_layout);

    
    LatticeColourMatrix Foo(&Fine);
    LatticeColourMatrix Bar(&Fine);

    LatticeSpinColourMatrix scFoo(&Fine);
    LatticeSpinColourMatrix scBar(&Fine);

    LatticeColourMatrix Shifted(&Fine);
    LatticeColourMatrix ShiftedCheck(&Fine);
    LatticeColourMatrix rShifted(&rbFine);
    LatticeColourMatrix bShifted(&rbFine);
   
    LatticeColourMatrix rFoo(&rbFine);
    LatticeColourMatrix bFoo(&rbFine);
    
    LatticeColourMatrix FooBar(&Fine);
    LatticeSpinColourMatrix scFooBar(&Fine);
    
    LatticeColourVector     cVec(&Fine);
    LatticeSpinVector       sVec(&Fine);
    LatticeSpinColourVector scVec(&Fine);

    LatticeColourMatrix     cMat(&Fine);
    LatticeSpinMatrix       sMat(&Fine);
    LatticeSpinColourMatrix scMat(&Fine);
    
    LatticeComplex scalar(&Fine);

    SpinMatrix GammaFive;
    iSpinMatrix<vComplex> iGammaFive;
    ColourMatrix cmat;
    
    random(Foo);
    gaussian(Bar);
    random(scFoo);
    random(scBar);

    random(cMat);
    random(sMat);
    random(scMat);
    random(cVec);
    random(sVec);
    random(scVec);

    fflush(stdout);
    cVec = cMat * cVec;  // LatticeColourVector     = LatticeColourMatrix     * LatticeColourVector
    sVec = sMat * sVec;  // LatticeSpinVector       = LatticeSpinMatrix       * LatticeSpinVector
    scVec= scMat * scVec;// LatticeSpinColourVector = LatticeSpinColourMatrix * LatticeSpinColourVector
    scVec= cMat * scVec; // LatticeSpinColourVector = LatticeColourMatrix     * LatticeSpinColourVector
    scVec= sMat * scVec; // LatticeSpinColourVector = LatticeSpinMatrix       * LatticeSpinColourVector
    
    cMat = outerProduct(cVec,cVec);
    scalar = localInnerProduct(cVec,cVec);
    
    scMat = sMat*scMat;  // LatticeSpinColourMatrix = LatticeSpinMatrix       * LatticeSpinColourMatrix

    // Non-lattice (const objects) * Lattice
    ColourMatrix cm;
    SpinColourMatrix scm;
    
    scMat = cMat*scMat;
    scm = cm * scm;         // SpinColourMatrix  = ColourMatrix     * SpinColourMatrix
    scm = scm *cm;          // SpinColourMatrix  = SpinColourMartix * ColourMatrix
    scm = GammaFive * scm ; // SpinColourMatrix  = SpinMatrix       * SpinColourMatrix
    scm = scm* GammaFive  ; // SpinColourMatrix  = SpinColourMatrix * SpinMatrix
    
    sMat = adj(sMat);       // LatticeSpinMatrix adjoint
    sMat = iGammaFive*sMat; // SpinMatrix * LatticeSpinMatrix
    sMat = GammaFive*sMat;  // SpinMatrix * LatticeSpinMatrix
    scMat= adj(scMat);
    cMat= adj(cMat);
    cm=adj(cm);
    scm=adj(scm);
    
//    Foo = Foo+scalar; // LatticeColourMatrix+Scalar
//    Foo = Foo*scalar; // LatticeColourMatrix*Scalar
//    Foo = Foo-scalar; // LatticeColourMatrix-Scalar
//    Foo = scalar*Foo; // Scalar*LatticeColourMatrix
//    Foo = scalar+Foo; // Scalar+LatticeColourMatrix
//    Foo = scalar-Foo; // Scalar-LatticeColourMatrix
    
    LatticeComplex trscMat(&Fine);
    trscMat = trace(scMat); // Trace
    
    FooBar = Bar;
    int shift=1;
    int dir=0;
 
    random(Foo);
    pickCheckerboard(CbRed,rFoo,Foo);    // Pick out red or black checkerboards
    pickCheckerboard(CbBlack,bFoo,Foo);

    Shifted  = Cshift(Foo,dir,shift);    // Shift everything
    bShifted = Cshift(rFoo,dir,shift);   // Shift red
    rShifted = Cshift(bFoo,dir,shift);   // Shift black
    
    setCheckerboard(ShiftedCheck,bShifted); // Put them all together
    setCheckerboard(ShiftedCheck,rShifted); // and check the results (later)

    // Lattice SU(3) x SU(3)
    FooBar = Foo * Bar;
    
    // Lattice 12x12 GEMM
    scFooBar = scFoo * scBar;
    
    // Benchmark some simple operations LatticeSU3 * Lattice SU3.
    double t0,t1,flops;
    double bytes;
    int ncall=100;
    int Nc = dpo::QCD::Nc;

    flops = ncall*1.0*volume*(8*Nc*Nc*Nc);
    bytes = ncall*1.0*volume*Nc*Nc    *2*3*sizeof(dpo::Real);
    printf("%f flop and %f bytes\n",flops,bytes/ncall);
        FooBar = Foo * Bar;
    t0=usecond();
    for(int i=0;i<ncall;i++){
        mult(FooBar,Foo,Bar); // this is better
    }
    t1=usecond();
#ifdef OMP
    printf("mult NumThread %d , Lattice size %d , %f us per call\n",omp_get_max_threads(),lat,(t1-t0)/ncall);
#endif
    printf("mult NumThread %d , Lattice size %d , %f Mflop/s\n",omp,lat,flops/(t1-t0));
    printf("mult NumThread %d , Lattice size %d , %f MB/s\n",omp,lat,bytes/(t1-t0));

/*
        mult(FooBar,Foo,Bar);
	//        FooBar = Foo * Bar;
    t0=usecond();
    for(int i=0;i<ncall;i++){
      //      mult(FooBar,Foo,Cshift(Bar,1,-2));
      mult(FooBar,Foo,Bar);
	//        FooBar = Foo * Bar; // this is bad
    }
    t1=usecond();
    
    printf("A: NumThread %d , Lattice size %d , %f us per call\n",omp,lat,(t1-t0)/ncall);
    printf("A: NumThread %d , Lattice size %d , %f Mflop/s\n",omp,lat,flops/(t1-t0));
    printf("A: NumThread %d , Lattice size %d , %f MB/s\n",omp,lat,bytes/(t1-t0));
*/

    pickCheckerboard(0,rFoo,FooBar);
    pickCheckerboard(1,bFoo,FooBar);
    
    setCheckerboard(FooBar,rFoo);
    setCheckerboard(FooBar,bFoo);
    
    double nrm=0;
    
    std::vector<int> coor(4);
    for(coor[0]=0;coor[0]<latt_size[0];coor[0]++){
    for(coor[1]=0;coor[1]<latt_size[1];coor[1]++){
    for(coor[2]=0;coor[2]<latt_size[2];coor[2]++){
    for(coor[3]=0;coor[3]<latt_size[3];coor[3]++){
 
        std::complex<dpo::Real> diff;
                    
        std::vector<int> shiftcoor = coor;
        shiftcoor[dir]=(shiftcoor[dir]-shift+latt_size[dir])%latt_size[dir];
                    
        ColourMatrix foo;
        ColourMatrix bar;
        ColourMatrix shifted1;
        ColourMatrix shifted2;
        ColourMatrix shifted3;
        ColourMatrix foobar1;
        ColourMatrix foobar2;
        ColourMatrix mdiff,amdiff;
                    
        peekSite(shifted1,Shifted,coor);
        peekSite(shifted2,Foo,shiftcoor);
        peekSite(shifted3,ShiftedCheck,coor);
       
        peekSite(foo,Foo,coor);
        
        mdiff = shifted1-shifted2;
        amdiff=adj(mdiff);
        ColourMatrix prod = amdiff*mdiff;
        TReal Ttr=real(trace(prod));
        double nn=Ttr._internal._internal;
        if ( nn > 0 )
            cout<<"Shift real trace fail "<<coor[0]<<coor[1]<<coor[2]<<coor[3] <<endl;
        
        for(int r=0;r<3;r++){
        for(int c=0;c<3;c++){
            diff =shifted1._internal._internal[r][c]-shifted2._internal._internal[r][c];
            nn=real(conj(diff)*diff);
            if ( nn > 0 )
                cout<<"Shift fail "<<coor[0]<<coor[1]<<coor[2]<<coor[3] <<" "
                    <<shifted1._internal._internal[r][c]<<" "<<shifted2._internal._internal[r][c]
                    << " "<< foo._internal._internal[r][c]<<endl;
            else if(0)
                cout<<"Shift pass "<<coor[0]<<coor[1]<<coor[2]<<coor[3] <<" "
                    <<shifted1._internal._internal[r][c]<<" "<<shifted2._internal._internal[r][c]
                    << " "<< foo._internal._internal[r][c]<<endl;
        }}
        
        for(int r=0;r<3;r++){
        for(int c=0;c<3;c++){
            diff =shifted3._internal._internal[r][c]-shifted2._internal._internal[r][c];
            nn=real(conj(diff)*diff);
            if ( nn > 0 )
                cout<<"Shift rb fail "<<coor[0]<<coor[1]<<coor[2]<<coor[3] <<" "
                <<shifted3._internal._internal[r][c]<<" "<<shifted2._internal._internal[r][c]
                << " "<< foo._internal._internal[r][c]<<endl;
            else if(0)
                cout<<"Shift rb pass "<<coor[0]<<coor[1]<<coor[2]<<coor[3] <<" "
                <<shifted3._internal._internal[r][c]<<" "<<shifted2._internal._internal[r][c]
                << " "<< foo._internal._internal[r][c]<<endl;
        }}
        peekSite(bar,Bar,coor);
                    
        peekSite(foobar1,FooBar,coor);
        foobar2 = foo*bar;
        for(int r=0;r<Nc;r++){
        for(int c=0;c<Nc;c++){
            diff =foobar2._internal._internal[r][c]-foobar1._internal._internal[r][c];
            nrm = nrm + real(conj(diff)*diff);
        }}
    }}}}

    std::cout << "LatticeColorMatrix * LatticeColorMatrix nrm diff = "<<nrm<<std::endl;

   } // loop for lat
 } // loop for omp

}
