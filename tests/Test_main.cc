    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_main.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#include "Grid.h"


using namespace std;
using namespace Grid;
using namespace Grid::QCD;

/*
 Grid_main.cc(232): error: no suitable user-defined conversion from "Grid::iScalar<Grid::iMatrix<Grid::iScalar<Grid::Complex>, 4>>" to "const Grid::iScalar<Grid::iScalar<Grid::iMatrix<Grid::Complex, 3>>>" exists
c_m = peekIdiot<SpinColourMatrix>(scm,1,2);
*/
template<class vobj> auto peekIdiot(const vobj &rhs,int i,int j) -> decltype(peekIndex<2>(rhs,0,0))
{
  return peekIndex<2>(rhs,i,j);
}
template<class vobj> auto peekDumKopf(const vobj &rhs,int i,int j) -> decltype(peekIndex<3>(rhs,0,0))
{
  return peekIndex<3>(rhs,i,j);
}
template<class vobj> auto peekDumKopf(const vobj &rhs,int i) -> decltype(peekIndex<3>(rhs,0))
{
  return peekIndex<3>(rhs,i);
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();

  latt_size.resize(4);

#ifdef AVX512
 for(int omp=128;omp<236;omp+=16){
#else
 for(int omp=1;omp<2;omp*=20){
#endif

#ifdef OMP
   omp_set_num_threads(omp);
#endif 

  for(int lat=8;lat<=16;lat+=40){

    std::cout << "Lat "<<lat<<std::endl;

    latt_size[0] = lat;
    latt_size[1] = lat;
    latt_size[2] = lat;
    latt_size[3] = lat;
    double volume = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
    
    GridCartesian           Fine(latt_size,simd_layout,mpi_layout);
    GridRedBlackCartesian rbFine(latt_size,simd_layout,mpi_layout);
    GridParallelRNG       FineRNG(&Fine);
    GridSerialRNG       SerialRNG;
    GridSerialRNG       SerialRNG1;

    FineRNG.SeedRandomDevice();
    SerialRNG.SeedRandomDevice();

    std::cout <<"SerialRNG" << SerialRNG._generators[0] <<std::endl;

    std::vector<typename GridSerialRNG::RngStateType> saved;
    SerialRNG.GetState(saved,0);
    SerialRNG1.SetState(saved,0);

    RealD dd1,dd2;

    std::cout << "Testing RNG state save restore"<<std::endl;
    for(int i=0;i<10;i++){
      random(SerialRNG,dd1);
      random(SerialRNG1,dd2);
      std::cout << "i "<<i<<" "<<dd1<< " " <<dd2<<std::endl;
    }
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
    
    LatticeLorentzColourMatrix lcMat(&Fine);


    LatticeComplex scalar(&Fine);
    LatticeReal    rscalar(&Fine);
    LatticeReal    iscalar(&Fine);

    SpinMatrix GammaFive;
    iSpinMatrix<vComplex> iGammaFive;
    ColourMatrix cmat;
    
    random(FineRNG,Foo);
    gaussian(FineRNG,Bar);
    random(FineRNG,scFoo);
    random(FineRNG,scBar);

    random(FineRNG,cMat);
    random(FineRNG,sMat);
    random(FineRNG,scMat);
    random(FineRNG,lcMat);
    random(FineRNG,cVec);
    random(FineRNG,sVec);
    random(FineRNG,scVec);


    fflush(stdout);
    
    TComplex tr = trace(cmat);

     
    cVec = cMat * cVec;  // LatticeColourVector     = LatticeColourMatrix     * LatticeColourVector
    sVec = sMat * sVec;  // LatticeSpinVector       = LatticeSpinMatrix       * LatticeSpinVector
    scVec= scMat * scVec;// LatticeSpinColourVector = LatticeSpinColourMatrix * LatticeSpinColourVector
    scVec= cMat * scVec; // LatticeSpinColourVector = LatticeColourMatrix     * LatticeSpinColourVector
    scVec= sMat * scVec; // LatticeSpinColourVector = LatticeSpinMatrix       * LatticeSpinColourVector
    
    cMat = outerProduct(cVec,cVec);
    scalar = localInnerProduct(cVec,cVec);

    cMat = Ta(cMat);  //traceless antihermitian


    scalar += scalar;
    scalar -= scalar;
    scalar *= scalar;
    add(scalar,scalar,scalar);
    sub(scalar,scalar,scalar);
    mult(scalar,scalar,scalar);

    mac(scalar,scalar,scalar);
    scalar = scalar+scalar;
    scalar = scalar-scalar;
    scalar = scalar*scalar;

    scalar=outerProduct(scalar,scalar);

    scalar=adj(scalar);

    //    rscalar=real(scalar);
    //    iscalar=imag(scalar);
    //    scalar =cmplx(rscalar,iscalar);
    PokeIndex<ColourIndex>(cVec,scalar,1);


    scalar=transpose(scalar);
    scalar=TransposeIndex<ColourIndex>(scalar);
    scalar=TraceIndex<SpinIndex>(scalar);
    scalar=PeekIndex<ColourIndex>(cVec,0);

    scalar=trace(scalar);
    scalar=localInnerProduct(cVec,cVec);
    scalar=localNorm2(cVec);
    
//     -=,+=,*=,()
//     add,+,sub,-,mult,mac,*
//     adj,conjugate
//     real,imag
//     transpose,transposeIndex  
//     trace,traceIndex
//     peekIndex
//     innerProduct,outerProduct,
//     localNorm2
//     localInnerProduct
    
  
    scMat = sMat*scMat;  // LatticeSpinColourMatrix = LatticeSpinMatrix       * LatticeSpinColourMatrix

    
    ///////////////////////
    // Non-lattice (const objects) * Lattice
    ColourMatrix cm;
    SpinColourMatrix scm;
    vSpinColourMatrix vscm;
    Complex cplx(1.0);
    Integer myint=1;
    double mydouble=1.0;

    //    vSpinColourMatrix vscm;
    scMat = cMat*scMat;
    scm = cm * scm;         // SpinColourMatrix  = ColourMatrix     * SpinColourMatrix
    scm = scm *cm;          // SpinColourMatrix  = SpinColourMartix * ColourMatrix
    scm = GammaFive * scm ; // SpinColourMatrix  = SpinMatrix       * SpinColourMatrix
    scm = scm* GammaFive  ; // SpinColourMatrix  = SpinColourMatrix * SpinMatrix

    scm = scm*cplx;
    vscm = vscm*cplx;
    scMat = scMat*cplx;

    scm = cplx*scm;
    vscm = cplx*vscm;
    scMat = cplx*scMat;
    scm = myint*scm;
    vscm = myint*vscm;
    scMat = scMat*myint;
    
    scm = scm*mydouble;
    vscm = vscm*mydouble;
    scMat = scMat*mydouble;
    scMat = mydouble*scMat;
    cMat = mydouble*cMat;
  
    sMat = adj(sMat);       // LatticeSpinMatrix adjoint
    sMat = iGammaFive*sMat; // SpinMatrix * LatticeSpinMatrix
    sMat = GammaFive*sMat;  // SpinMatrix * LatticeSpinMatrix
    scMat= adj(scMat);
    cMat= adj(cMat);
    cm=adj(cm);
    scm=adj(scm);
    scm=transpose(scm);
    scm=transposeIndex<1>(scm);
    
   
    random(SerialRNG, cm);
    std::cout<<GridLogMessage << cm << std::endl;

    cm = Ta(cm);
    TComplex tracecm= trace(cm);      
    std::cout<<GridLogMessage << cm << std::endl;


    cm = Exponentiate(cm, 2.0, 12);
    std::cout<<GridLogMessage << cm << "  " << std::endl;
    Complex det = Determinant(cm);
    std::cout<<GridLogMessage << "determinant: " << det <<  std::endl;
    std::cout<<GridLogMessage << "norm: " << norm2(cm) <<  std::endl;

    cm = ProjectOnGroup(cm);
    std::cout<<GridLogMessage << cm << "  " << std::endl;
    std::cout<<GridLogMessage << "norm: " << norm2(cm) <<  std::endl;
    cm = ProjectOnGroup(cm);
    std::cout<<GridLogMessage << cm << "  " << std::endl;
    std::cout<<GridLogMessage << "norm: " << norm2(cm) <<  std::endl;


    //    det = Determinant(cm);
    //    std::cout<<GridLogMessage << "determinant: " << det <<  std::endl;


//    Foo = Foo+scalar; // LatticeColourMatrix+Scalar
//    Foo = Foo*scalar; // LatticeColourMatrix*Scalar
//    Foo = Foo-scalar; // LatticeColourMatrix-Scalar
//    Foo = scalar*Foo; // Scalar*LatticeColourMatrix
//    Foo = scalar+Foo; // Scalar+LatticeColourMatrix
//    Foo = scalar-Foo; // Scalar-LatticeColourMatrix
    
    LatticeComplex trscMat(&Fine);
    trscMat = trace(scMat); // Trace

    // Exponentiate test
    std::vector<int> mysite {0,0,0,0};
    random(FineRNG,cMat);
    cMat = Ta(cMat);
    peekSite(cm, cMat, mysite);
    std::cout<<GridLogMessage << cm << "  " << std::endl;
    cm = Exponentiate(cm, 1.0, 12);
    std::cout<<GridLogMessage << cm << "  " << std::endl;
    std::cout<<GridLogMessage << "norm: " << norm2(cm) <<  std::endl;


    std::cout<<GridLogMessage << "norm cMmat : " << norm2(cMat) <<  std::endl;
    cMat = expMat(cMat, ComplexD(1.0, 0.0));
    std::cout<<GridLogMessage << "norm expMat: " << norm2(cMat) <<  std::endl;
    peekSite(cm, cMat, mysite);
    std::cout<<GridLogMessage << cm << "  " << std::endl;
    std::cout<<GridLogMessage << "determinant: " << Determinant(cm) <<  std::endl;
    std::cout<<GridLogMessage << "norm: " << norm2(cm) <<  std::endl;


    // LatticeComplex trlcMat(&Fine);
    // trlcMat = trace(lcMat); // Trace involving iVector - now generates error
    

    { // Peek-ology and Poke-ology, with a little app-ology
      Complex      c;
      ColourMatrix c_m;   
      SpinMatrix   s_m;   
      SpinColourMatrix sc_m; 

      s_m = TensorIndexRecursion<ColourIndex>::traceIndex(sc_m); // Map to traceColour
      c_m = TensorIndexRecursion<SpinIndex>::traceIndex(sc_m); // map to traceSpin

      c   = TensorIndexRecursion<SpinIndex>::traceIndex(s_m); 
      c   = TensorIndexRecursion<ColourIndex>::traceIndex(c_m);
      
      s_m = TensorIndexRecursion<ColourIndex>::peekIndex(scm,0,0);
      c_m = TensorIndexRecursion<SpinIndex>::peekIndex(scm,1,2);
      //      c_m = peekSpin<SpinColourMatrix>(scm,1,2);
      //      c_m = peekIdiot<SpinColourMatrix>(scm,1,2);

      printf("c. Level %d\n",c_m.TensorLevel);
      printf("c. Level %d\n",c_m().TensorLevel);
      printf("c. Level %d\n",c_m()().TensorLevel);
      
      c_m()()    = scm()(0,0); //ColourComponents of CM <= ColourComponents of SpinColourMatrix
      scm()(1,1) = cm()();  //ColourComponents of CM <= ColourComponents of SpinColourMatrix
      c          = scm()(1,1)(1,2);
      scm()(1,1)(2,1) = c;

      //      pokeIndex<ColourIndex> (c_m,c,0,0);
    }

    FooBar = Bar;
    /*
    { 
      std::vector<int> coor(4);
      for(int d=0;d<4;d++) coor[d] = 0;
      peekSite(cmat,Foo,coor);
      Foo = zero;
      pokeSite(cmat,Foo,coor);
    }
    random(Foo);
    */
    lex_sites(Foo);


    Integer mm[4];
    mm[0]=1;
    mm[1]=Fine._rdimensions[0];
    mm[2]=Fine._ldimensions[0]*Fine._ldimensions[1];
    mm[3]=Fine._ldimensions[0]*Fine._ldimensions[1]*Fine._ldimensions[2];

    LatticeInteger lex(&Fine);
    lex=zero;
    for(int d=0;d<4;d++){
      LatticeInteger coor(&Fine);
      LatticeCoordinate(coor,d);
      lex = lex + coor*mm[d];
      
    }



    //    Bar = zero;
    //    Bar = where(lex<Integer(10),Foo,Bar);


    cout << "peeking sites..\n";
    {
      std::vector<int> coor(4);
      for(coor[3]=0;coor[3]<latt_size[3]/mpi_layout[3];coor[3]++){
      for(coor[2]=0;coor[2]<latt_size[2]/mpi_layout[2];coor[2]++){
      for(coor[1]=0;coor[1]<latt_size[1]/mpi_layout[1];coor[1]++){
      for(coor[0]=0;coor[0]<latt_size[0]/mpi_layout[0];coor[0]++){
        ColourMatrix bar;
        peekSite(bar,Bar,coor);
        for(int r=0;r<3;r++){
        for(int c=0;c<3;c++){
	  //	  cout<<"bar "<<coor[0]<<coor[1]<<coor[2]<<coor[3] <<" "<<bar()()(r,c)<<std::endl;
	}}
      }}}}
    } 

    //setCheckerboard(ShiftedCheck,rFoo); 
    //setCheckerboard(ShiftedCheck,bFoo); 


    // Lattice SU(3) x SU(3)
    Fine.Barrier();
    FooBar = Foo * Bar;

    // Lattice 12x12 GEMM
    scFooBar = scFoo * scBar;

    // Benchmark some simple operations LatticeSU3 * Lattice SU3.
    double t0,t1,flops;
    double bytes;
    int ncall=5000;
    int Nc = Grid::QCD::Nc;

    LatticeGaugeField U(&Fine);
    //    LatticeColourMatrix Uy = peekLorentz(U,1);
    //    LatticeColourMatrix Uy = peekDumKopf(U,1);

    flops = ncall*1.0*volume*(8*Nc*Nc*Nc);
    bytes = ncall*1.0*volume*Nc*Nc    *2*3*sizeof(Grid::Real);
    if ( Fine.IsBoss() ) {
      printf("%f flop and %f bytes\n",flops,bytes/ncall);
    }
    FooBar = Foo * Bar;
    Fine.Barrier();
    t0=usecond();
    for(int i=0;i<ncall;i++){
      Fine.Barrier();
      mult(FooBar,Foo,Bar); // this is better
    }

    t1=usecond();
    Fine.Barrier();
    if ( Fine.IsBoss() ) {
#ifdef OMP
      printf("mult NumThread %d , Lattice size %d , %f us per call\n",omp_get_max_threads(),lat,(t1-t0)/ncall);
#endif
      printf("mult NumThread %d , Lattice size %d , %f us per call\n",omp,lat,(t1-t0)/ncall);
      printf("mult NumThread %d , Lattice size %d , %f Mflop/s\n",omp,lat,flops/(t1-t0));
      printf("mult NumThread %d , Lattice size %d , %f MB/s\n",omp,lat,bytes/(t1-t0));
    }
    mult(FooBar,Foo,Bar);
    FooBar = Foo * Bar;

    bytes = ncall*1.0*volume*Nc*Nc    *2*5*sizeof(Grid::Real);
    Fine.Barrier();
    t0=usecond();
    for(int i=0;i<ncall;i++){
      Fine.Barrier();
      mult(FooBar,Foo,Cshift(Bar,1,-1));
      //mult(FooBar,Foo,Bar);
      //FooBar = Foo * Bar; // this is bad
    }
    t1=usecond();
    Fine.Barrier();

    FooBar = Foo * Bar;
    
    if ( Fine.IsBoss() ) {
      printf("Cshift Mult: NumThread %d , Lattice size %d , %f us per call\n",omp,lat,(t1-t0)/ncall);
      printf("Cshift Mult: NumThread %d , Lattice size %d , %f Mflop/s\n",omp,lat,flops/(t1-t0));
      printf("Cshift Mult: NumThread %d , Lattice size %d , %f MB/s\n",omp,lat,bytes/(t1-t0));
    }
    //    pickCheckerboard(0,rFoo,FooBar);
    //    pickCheckerboard(1,bFoo,FooBar);
    //    setCheckerboard(FooBar,rFoo);
    //    setCheckerboard(FooBar,bFoo);
    
    double nrm=0;


    LatticeColourMatrix deriv(&Fine);
    double half=0.5;
    deriv = 0.5*Cshift(Foo,0,1) - 0.5*Cshift(Foo,0,-1);


    for(int dir=0;dir<4;dir++){
      for(int shift=0;shift<latt_size[dir];shift++){

	pickCheckerboard(0,rFoo,Foo);    // Pick out red or black checkerboards
	pickCheckerboard(1,bFoo,Foo);
    
	if ( Fine.IsBoss() ) {
	  std::cout<<GridLogMessage << "Shifting both parities by "<< shift <<" direction "<< dir <<std::endl;
	}
	Shifted  = Cshift(Foo,dir,shift);    // Shift everything

	bShifted = Cshift(rFoo,dir,shift);   // Shift red->black
	rShifted = Cshift(bFoo,dir,shift);   // Shift black->red
    
	ShiftedCheck=zero;
	setCheckerboard(ShiftedCheck,bShifted); // Put them all together
	setCheckerboard(ShiftedCheck,rShifted); // and check the results (later)
    
	// Check results
	std::vector<int> coor(4);
	for(coor[3]=0;coor[3]<latt_size[3]/mpi_layout[3];coor[3]++){
	for(coor[2]=0;coor[2]<latt_size[2]/mpi_layout[2];coor[2]++){
	for(coor[1]=0;coor[1]<latt_size[1]/mpi_layout[1];coor[1]++){
	for(coor[0]=0;coor[0]<latt_size[0]/mpi_layout[0];coor[0]++){

	 
 
        std::complex<Grid::Real> diff;
                    
        std::vector<int> shiftcoor = coor;
        shiftcoor[dir]=(shiftcoor[dir]+shift+latt_size[dir])%(latt_size[dir]/mpi_layout[dir]);

	std::vector<int> rl(4);
	for(int dd=0;dd<4;dd++){
	  rl[dd] = latt_size[dd]/simd_layout[dd]/mpi_layout[dd];
	}
	int lex =  coor[0]%rl[0]
	  + (coor[1]%rl[1])*rl[0]
	  + (coor[2]%rl[2])*rl[0]*rl[1]
	  + (coor[3]%rl[3])*rl[0]*rl[1]*rl[2];
	lex += 
	  +1000*(coor[0]/rl[0])
	  +1000*(coor[1]/rl[1])*simd_layout[0]
	  +1000*(coor[2]/rl[2])*simd_layout[0]*simd_layout[1]
	  +1000*(coor[3]/rl[3])*simd_layout[0]*simd_layout[1]*simd_layout[2];

	int lex_coor = shiftcoor[0]%rl[0]
	  + (shiftcoor[1]%rl[1])*rl[0]
	  + (shiftcoor[2]%rl[2])*rl[0]*rl[1]
	  + (shiftcoor[3]%rl[3])*rl[0]*rl[1]*rl[2];
	lex_coor += 
	  +1000*(shiftcoor[0]/rl[0])
	  +1000*(shiftcoor[1]/rl[1])*simd_layout[0]
	  +1000*(shiftcoor[2]/rl[2])*simd_layout[0]*simd_layout[1]
	  +1000*(shiftcoor[3]/rl[3])*simd_layout[0]*simd_layout[1]*simd_layout[2];

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
	Complex trprod = trace(prod);
        Real Ttr=real(trprod);
        double nn=Ttr;
        if ( nn > 0 )
            cout<<"Shift real trace fail "<<coor[0]<<coor[1]<<coor[2]<<coor[3] <<endl;
     

        for(int r=0;r<3;r++){
        for(int c=0;c<3;c++){
            diff =shifted1()()(r,c)-shifted2()()(r,c);
            nn=real(conjugate(diff)*diff);
            if ( nn > 0 )
                cout<<"Shift fail (shifted1/shifted2-ref) "<<coor[0]<<coor[1]<<coor[2]<<coor[3] <<" "
                    <<shifted1()()(r,c)<<" "<<shifted2()()(r,c)
                    << " "<< foo()()(r,c)<< " lex expect " << lex_coor << " lex "<<lex<<endl;
            else if(0)
                cout<<"Shift pass 1vs2 "<<coor[0]<<coor[1]<<coor[2]<<coor[3] <<" "
                    <<shifted1()()(r,c)<<" "<<shifted2()()(r,c)
                    << " "<< foo()()(r,c)<< " lex expect " << lex_coor << " lex "<<lex<<endl;
        }}
        
        for(int r=0;r<3;r++){
        for(int c=0;c<3;c++){
            diff =shifted3()()(r,c)-shifted2()()(r,c);
            nn=real(conjugate(diff)*diff);
            if ( nn > 0 )
                cout<<"Shift rb fail (shifted3/shifted2-ref) "<<coor[0]<<coor[1]<<coor[2]<<coor[3] <<" "
                <<shifted3()()(r,c)<<" "<<shifted2()()(r,c)
                << " "<< foo()()(r,c)<< " lex expect " << lex_coor << " lex "<<lex<<endl;
            else if(0)
                cout<<"Shift rb pass 3vs2 "<<coor[0]<<coor[1]<<coor[2]<<coor[3] <<" "
                <<shifted3()()(r,c)<<" "<<shifted2()()(r,c)
                << " "<< foo()()(r,c)<< " lex expect " << lex_coor << " lex "<<lex<<endl;
        }}
        peekSite(bar,Bar,coor);
                    
        peekSite(foobar1,FooBar,coor);
        foobar2 = foo*bar;
        for(int r=0;r<Nc;r++){
        for(int c=0;c<Nc;c++){
            diff =foobar2()()(r,c)-foobar1()()(r,c);
            nrm = nrm + real(conjugate(diff)*diff);
        }}
    }}}}
	if( Fine.IsBoss() ){
	  std::cout<<GridLogMessage << "LatticeColorMatrix * LatticeColorMatrix nrm diff = "<<nrm<<std::endl;
	}
      }}

   } // loop for lat
 } // loop for omp


 std::cout<<GridLogMessage << sizeof(vComplexF) << std::endl;
 
 Grid_finalize();
}
