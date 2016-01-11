    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/utils/SUn.h

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
#ifndef QCD_UTIL_SUN_H
#define QCD_UTIL_SUN_H

namespace Grid {
  namespace QCD {

template<int ncolour>
class SU {
public:

  static int generators(void)   { return ncolour*ncolour-1; }
  static int su2subgroups(void) { return (ncolour*(ncolour-1))/2; }

  template<typename vtype> using iSUnMatrix              = iScalar<iScalar<iMatrix<vtype, ncolour> > > ;
  template<typename vtype> using iSU2Matrix              = iScalar<iScalar<iMatrix<vtype, 2> > > ;

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Types can be accessed as SU<2>::Matrix , SU<2>::vSUnMatrix, SU<2>::LatticeMatrix etc...
  //////////////////////////////////////////////////////////////////////////////////////////////////
  typedef iSUnMatrix<Complex>         Matrix;
  typedef iSUnMatrix<ComplexF>        MatrixF;
  typedef iSUnMatrix<ComplexD>        MatrixD;

  typedef iSUnMatrix<vComplex>       vMatrix;
  typedef iSUnMatrix<vComplexF>      vMatrixF;
  typedef iSUnMatrix<vComplexD>      vMatrixD;

  typedef Lattice<vMatrix>     LatticeMatrix;
  typedef Lattice<vMatrixF>    LatticeMatrixF;
  typedef Lattice<vMatrixD>    LatticeMatrixD;

  typedef iSU2Matrix<Complex>         SU2Matrix;
  typedef iSU2Matrix<ComplexF>        SU2MatrixF;
  typedef iSU2Matrix<ComplexD>        SU2MatrixD;

  typedef iSU2Matrix<vComplex>       vSU2Matrix;
  typedef iSU2Matrix<vComplexF>      vSU2MatrixF;
  typedef iSU2Matrix<vComplexD>      vSU2MatrixD;

  typedef Lattice<vSU2Matrix>     LatticeSU2Matrix;
  typedef Lattice<vSU2MatrixF>    LatticeSU2MatrixF;
  typedef Lattice<vSU2MatrixD>    LatticeSU2MatrixD;


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
  template<class cplx>
  static void generator(int lieIndex,iSUnMatrix<cplx> &ta){
    // map lie index to which type of generator
    int diagIndex;
    int su2Index;
    int sigxy;
    int NNm1 =  ncolour*(ncolour-1);
    if ( lieIndex>= NNm1 ) {
      diagIndex = lieIndex -NNm1;
      generatorDiagonal(diagIndex,ta);
      return;
    }
    sigxy   = lieIndex&0x1;
    su2Index= lieIndex>>1;
    if ( sigxy ) generatorSigmaY(su2Index,ta);
    else         generatorSigmaX(su2Index,ta);
  }
  template<class cplx>
  static void generatorSigmaX(int su2Index,iSUnMatrix<cplx> &ta){
    ta=zero;
    int i1,i2;
    su2SubGroupIndex(i1,i2,su2Index);
    ta()()(i1,i2)=1.0;
    ta()()(i2,i1)=1.0;
    ta= ta *0.5;
  }
  template<class cplx>
  static void generatorSigmaY(int su2Index,iSUnMatrix<cplx> &ta){
    ta=zero;
    cplx i(0.0,1.0);
    int i1,i2;
    su2SubGroupIndex(i1,i2,su2Index);
    ta()()(i1,i2)=-i;
    ta()()(i2,i1)= i;
    ta= ta *0.5;
  }
  template<class cplx>
  static void generatorDiagonal(int diagIndex,iSUnMatrix<cplx> &ta){
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
  // Map a su2 subgroup number to the pair of rows that are non zero
  ////////////////////////////////////////////////////////////////////////
  static void su2SubGroupIndex(int &i1,int &i2,int su2_index){

    assert( (su2_index>=0) && (su2_index< (ncolour*(ncolour-1))/2) );

    int spare=su2_index;
    for(i1=0;spare >= (ncolour-1-i1);i1++ ){
      spare = spare - (ncolour-1-i1);  // remove the Nc-1-i1 terms  
    }
    i2=i1+1+spare;
  }

  //////////////////////////////////////////////////////////////////////////////////////////
  // Pull out a subgroup and project on to real coeffs x pauli basis
  //////////////////////////////////////////////////////////////////////////////////////////
  template<class vcplx>
  static void su2Extract(    Lattice<iSinglet<vcplx> > &Determinant,
			     Lattice<iSU2Matrix<vcplx> > &subgroup, 
			     const Lattice<iSUnMatrix<vcplx> > &source, 
			     int su2_index)
  {
    GridBase *grid(source._grid);
    conformable(subgroup,source);
    conformable(subgroup,Determinant);
    int i0,i1;    
    su2SubGroupIndex(i0,i1,su2_index);

PARALLEL_FOR_LOOP
    for(int ss=0;ss<grid->oSites();ss++){
      subgroup._odata[ss]()()(0,0) = source._odata[ss]()()(i0,i0);
      subgroup._odata[ss]()()(0,1) = source._odata[ss]()()(i0,i1);
      subgroup._odata[ss]()()(1,0) = source._odata[ss]()()(i1,i0);
      subgroup._odata[ss]()()(1,1) = source._odata[ss]()()(i1,i1);

      iSU2Matrix<vcplx> Sigma = subgroup._odata[ss];

      Sigma = Sigma-adj(Sigma)+trace(adj(Sigma));

      subgroup._odata[ss] = Sigma;

      // this should be purely real
      Determinant._odata[ss] = Sigma()()(0,0)*Sigma()()(1,1)
   	                     - Sigma()()(0,1)*Sigma()()(1,0);

    }

  }

  //////////////////////////////////////////////////////////////////////////////////////////
  // Set matrix to one and insert a pauli subgroup
  //////////////////////////////////////////////////////////////////////////////////////////
  template<class vcplx>
  static void su2Insert( const Lattice<iSU2Matrix<vcplx> > &subgroup, 
			       Lattice<iSUnMatrix<vcplx> > &dest,
			int su2_index)
  {
    GridBase *grid(dest._grid);
    conformable(subgroup,dest);
    int i0,i1;    
    su2SubGroupIndex(i0,i1,su2_index);

    dest = 1.0; // start out with identity
PARALLEL_FOR_LOOP
    for(int ss=0;ss<grid->oSites();ss++){
      dest._odata[ss]()()(i0,i0) = subgroup._odata[ss]()()(0,0);
      dest._odata[ss]()()(i0,i1) = subgroup._odata[ss]()()(0,1);
      dest._odata[ss]()()(i1,i0) = subgroup._odata[ss]()()(1,0);
      dest._odata[ss]()()(i1,i1) = subgroup._odata[ss]()()(1,1);
    }
  }


  ///////////////////////////////////////////////
  // Generate e^{ Re Tr Staple Link} dlink 
  //
  // *** Note Staple should be appropriate linear compbination between all staples.
  // *** If already by beta pass coefficient 1.0.
  // *** This routine applies the additional 1/Nc factor that comes after trace in action.
  //
  ///////////////////////////////////////////////
  static void SubGroupHeatBath(      GridSerialRNG       &sRNG,
				     GridParallelRNG     &pRNG,
				     RealD beta,//coeff multiplying staple in action (with no 1/Nc)
				     LatticeMatrix &link, 
				     const LatticeMatrix &barestaple, // multiplied by action coeffs so th
				     int su2_subgroup,
				     int nheatbath,
				     LatticeInteger &wheremask)
  {
    GridBase *grid = link._grid;

    int ntrials=0;
    int nfails=0;
    const RealD twopi=2.0*M_PI;

    LatticeMatrix staple(grid); 

    staple = barestaple * (beta/ncolour);

    LatticeMatrix    V(grid);    V = link*staple;

    // Subgroup manipulation in the lie algebra space
    LatticeSU2Matrix    u(grid);   // Kennedy pendleton "u" real projected normalised Sigma
    LatticeSU2Matrix uinv(grid);
    LatticeSU2Matrix   ua(grid);   // a in pauli form
    LatticeSU2Matrix    b(grid);   // rotated matrix after hb

    // Some handy constant fields
    LatticeComplex ones (grid); ones = 1.0;
    LatticeComplex zeros(grid); zeros=zero;
    LatticeReal rones (grid); rones = 1.0;
    LatticeReal rzeros(grid); rzeros=zero;
    LatticeComplex udet(grid); // determinant of real(staple)
    LatticeInteger mask_true (grid); mask_true = 1;
    LatticeInteger mask_false(grid); mask_false= 0;

  /*
PLB 156 P393 (1985) (Kennedy and Pendleton)

Note: absorb "beta" into the def of sigma compared to KP paper; staple
passed to this routine has "beta" already multiplied in

Action linear in links h and of form:

    beta S = beta  Sum_p (1 - 1/Nc Re Tr Plaq )

Writing Sigma = 1/Nc (beta Sigma') where sum over staples is "Sigma' "

     beta S = const - beta/Nc Re Tr h Sigma'
            = const - Re Tr h Sigma
      
Decompose h and Sigma into (1, sigma_j) ; h_i real, h^2=1, Sigma_i complex arbitrary.

    Tr h Sigma = h_i Sigma_j Tr (sigma_i sigma_j)  = h_i Sigma_j 2 delta_ij
 Re Tr h Sigma = 2 h_j Re Sigma_j

Normalised re Sigma_j = xi u_j 

With u_j a unit vector and U can be in SU(2);

Re Tr h Sigma = 2 h_j Re Sigma_j = 2 xi (h.u)

4xi^2 = Det [ Sig - Sig^dag  + 1 Tr Sigdag]
 u   = 1/2xi [ Sig - Sig^dag  + 1 Tr Sigdag]

 xi = sqrt(Det)/2;

Write a= u h in SU(2); a has pauli decomp a_j;

Note: Product b' xi is unvariant because scaling Sigma leaves
      normalised vector "u" fixed; Can rescale Sigma so b' = 1.
  */

    ////////////////////////////////////////////////////////
    // Real part of Pauli decomposition
    // Note a subgroup can project to zero in cold start
    ////////////////////////////////////////////////////////
    su2Extract(udet,u,V,su2_subgroup); 

    //////////////////////////////////////////////////////
    // Normalising this vector if possible; else identity
    //////////////////////////////////////////////////////
    LatticeComplex xi(grid);  

    LatticeSU2Matrix lident(grid);
    
    SU2Matrix ident = Complex(1.0); 
    SU2Matrix pauli1; SU<2>::generator(0,pauli1);
    SU2Matrix pauli2; SU<2>::generator(1,pauli2);
    SU2Matrix pauli3; SU<2>::generator(2,pauli3);
    pauli1 = timesI(pauli1)*2.0;
    pauli2 = timesI(pauli2)*2.0;
    pauli3 = timesI(pauli3)*2.0;

    LatticeComplex cone(grid);
    LatticeReal adet(grid);
    adet = abs(toReal(udet));
    lident=Complex(1.0);
    cone  =Complex(1.0);
    Real machine_epsilon=1.0e-7;
    u   = where(adet>machine_epsilon,u,lident);
    udet= where(adet>machine_epsilon,udet,cone);

    xi  = 0.5*sqrt(udet);     //4xi^2 = Det [ Sig - Sig^dag  + 1 Tr Sigdag]
    u   = 0.5*u*pow(xi,-1.0); //  u   = 1/2xi [ Sig - Sig^dag  + 1 Tr Sigdag]

    // Debug test for sanity
    uinv=adj(u);
    b=u*uinv-1.0;
    assert(norm2(b)<1.0e-4);

  /*
Measure: Haar measure dh has d^4a delta(1-|a^2|)
In polars:
  da = da0 r^2 sin theta dr dtheta dphi delta( 1 - r^2 -a0^2) 
     = da0 r^2 sin theta dr dtheta dphi delta( (sqrt(1-a0^) - r)(sqrt(1-a0^) + r) )
     = da0 r/2 sin theta dr dtheta dphi delta( (sqrt(1-a0^) - r) )

Action factor Q(h) dh  = e^-S[h]  dh =  e^{  xi Tr uh} dh    // beta enters through xi
                                     =  e^{2 xi (h.u)} dh
                                     =  e^{2 xi h0u0}.e^{2 xi h1u1}.e^{2 xi h2u2}.e^{2 xi h3u3} dh

Therefore for each site, take xi for that site
i) generate  |a0|<1 with dist 
   (1-a0^2)^0.5 e^{2 xi a0 } da0

Take alpha = 2 xi  = 2 xi [ recall 2 beta/Nc unmod staple norm]; hence 2.0/Nc factor in Chroma ]
A. Generate two uniformly distributed pseudo-random numbers R and R', R'', R''' in the unit interval; 
B. Set X = -(ln R)/alpha, X' =-(ln R')/alpha;
C. Set C = cos^2(2pi R"), with R" another uniform random number in [0,1] ;
D. Set A = XC;
E. Let d  = X'+A;
F. If R'''^2 :> 1 - 0.5 d,  go back to A;
G. Set a0 = 1 - d;

Note that in step D setting B ~ X - A and using B in place of A in step E will generate a second independent a 0 value.
  */

    /////////////////////////////////////////////////////////
    // count the number of sites by picking "1"'s out of hat
    /////////////////////////////////////////////////////////
    Integer hit=0;
    LatticeReal rtmp(grid);
    rtmp=where(wheremask,rones,rzeros);
    RealD numSites = sum(rtmp);
    RealD numAccepted;
    LatticeInteger      Accepted(grid); Accepted = zero;
    LatticeInteger newlyAccepted(grid);
   
    std::vector<LatticeReal> xr(4,grid);
    std::vector<LatticeReal>  a(4,grid);
    LatticeReal d(grid); d=zero;
    LatticeReal alpha(grid);

    //    std::cout<<GridLogMessage<<"xi "<<xi <<std::endl;
    alpha = toReal(2.0*xi);

    do { 

      // A. Generate two uniformly distributed pseudo-random numbers R and R', R'', R''' in the unit interval; 
      random(pRNG,xr[0]); 
      random(pRNG,xr[1]); 
      random(pRNG,xr[2]); 
      random(pRNG,xr[3]); 

      // B. Set X = - ln R/alpha, X' = -ln R'/alpha
      xr[1] = -log(xr[1])/alpha;  
      xr[2] = -log(xr[2])/alpha;

      // C. Set C = cos^2(2piR'')
      xr[3] = cos(xr[3]*twopi);
      xr[3] = xr[3]*xr[3];

      LatticeReal xrsq(grid);

      //D. Set A = XC;
      //E. Let d  = X'+A;
      xrsq = xr[2]+xr[1]*xr[3];

      d = where(Accepted,d,xr[2]+xr[1]*xr[3]);
	      
      //F. If R'''^2 :> 1 - 0.5 d,  go back to A;
      LatticeReal thresh(grid); thresh = 1.0-d*0.5;
      xrsq = xr[0]*xr[0];
      LatticeInteger ione(grid); ione = 1;
      LatticeInteger izero(grid); izero=zero;

      newlyAccepted = where ( xrsq < thresh,ione,izero);
      Accepted = where ( newlyAccepted, newlyAccepted,Accepted);
      Accepted = where ( wheremask, Accepted,izero);

      // FIXME need an iSum for integer to avoid overload on return type??
      rtmp=where(Accepted,rones,rzeros);
      numAccepted = sum(rtmp);

      hit++;

    } while ( (numAccepted < numSites) && ( hit < nheatbath) );

    // G. Set a0 = 1 - d;
    a[0]=zero;
    a[0]=where(wheremask,1.0-d,a[0]);

    //////////////////////////////////////////
    //    ii) generate a_i uniform on two sphere radius (1-a0^2)^0.5
    //////////////////////////////////////////

    LatticeReal a123mag(grid);
    a123mag = sqrt(abs(1.0-a[0]*a[0]));

    LatticeReal cos_theta (grid);
    LatticeReal sin_theta (grid);
    LatticeReal       phi (grid);

    random(pRNG,phi);            phi = phi * twopi;       // uniform in [0,2pi]
    random(pRNG,cos_theta); cos_theta=(cos_theta*2.0)-1.0;  // uniform in [-1,1]
    sin_theta = sqrt(abs(1.0-cos_theta*cos_theta));

    a[1]    = a123mag * sin_theta * cos(phi);
    a[2]    = a123mag * sin_theta * sin(phi);
    a[3]    = a123mag * cos_theta;
    
    ua = toComplex(a[0])*ident 
       + toComplex(a[1])*pauli1 
       + toComplex(a[2])*pauli2 
       + toComplex(a[3])*pauli3;

    b = 1.0;
    b = where(wheremask,uinv*ua,b);
    su2Insert(b,V,su2_subgroup);

    //mask the assignment back based on Accptance
    link = where(Accepted,V * link,link);

    //////////////////////////////
    // Debug Checks
    // SU2 check
    LatticeSU2Matrix    check(grid);   // rotated matrix after hb
    u = zero;
    check = ua * adj(ua) - 1.0;
    check = where(Accepted,check,u);
    assert(norm2(check)<1.0e-4);

    check = b*adj(b)-1.0;
    check = where(Accepted,check,u);
    assert(norm2(check)<1.0e-4);
    
    LatticeMatrix Vcheck(grid);
    Vcheck = zero;
    Vcheck = where(Accepted,V*adj(V) - 1.0,Vcheck);
    //    std::cout<<GridLogMessage << "SU3 check " <<norm2(Vcheck)<<std::endl;
    assert(norm2(Vcheck)<1.0e-4);
    
    // Verify the link stays in SU(3)
    //    std::cout<<GridLogMessage <<"Checking the modified link"<<std::endl;
    Vcheck = link*adj(link) - 1.0;
    assert(norm2(Vcheck)<1.0e-4);
    /////////////////////////////////
  }

  static void printGenerators(void)
  {
    for(int gen=0;gen<generators();gen++){
      Matrix ta;
      generator(gen,ta);
      std::cout<<GridLogMessage<< "Nc = "<<ncolour<<" t_"<<gen<<std::endl;
      std::cout<<GridLogMessage<<ta<<std::endl;
    }
  }

  static void testGenerators(void){
    Matrix ta;
    Matrix tb;
    std::cout<<GridLogMessage<<"Checking trace ta tb is 0.5 delta_ab"<<std::endl;
    for(int a=0;a<generators();a++){
      for(int b=0;b<generators();b++){
	generator(a,ta);
	generator(b,tb);
	Complex tr =TensorRemove(trace(ta*tb)); 
	std::cout<<GridLogMessage<<tr<<" ";
	if(a==b) assert(abs(tr-Complex(0.5))<1.0e-6);
	if(a!=b) assert(abs(tr)<1.0e-6);
      }
      std::cout<<GridLogMessage<<std::endl;
    }
    std::cout<<GridLogMessage<<"Checking hermitian"<<std::endl;
    for(int a=0;a<generators();a++){
      generator(a,ta);
      std::cout<<GridLogMessage<<a<<" ";
      assert(norm2(ta-adj(ta))<1.0e-6);
    }    
    std::cout<<GridLogMessage<<std::endl;

    std::cout<<GridLogMessage<<"Checking traceless"<<std::endl;
    for(int a=0;a<generators();a++){
      generator(a,ta);
      Complex tr =TensorRemove(trace(ta)); 
      std::cout<<GridLogMessage<<a<<" ";
      assert(abs(tr)<1.0e-6);
    }    
    std::cout<<GridLogMessage<<std::endl;
  }

  // reunitarise??
  static void LieRandomize(GridParallelRNG     &pRNG,LatticeMatrix &out,double scale=1.0){
    GridBase *grid = out._grid;

    LatticeComplex ca (grid);
    LatticeMatrix  lie(grid);
    LatticeMatrix  la (grid);
    Complex ci(0.0,scale);
    Complex cone(1.0,0.0);
    Matrix ta;

    lie=zero;
    for(int a=0;a<generators();a++){

      random(pRNG,ca); 

      ca = (ca+conjugate(ca))*0.5;
      ca = ca - 0.5;

      generator(a,ta);

      la=ci*ca*ta;

      lie = lie+la; // e^{i la ta}
    }
    taExp(lie,out);
  }

  static void GaussianLieAlgebraMatrix(GridParallelRNG     &pRNG,LatticeMatrix &out,double scale=1.0){
    GridBase *grid = out._grid;
    LatticeReal ca (grid);
    LatticeMatrix  la (grid);
    Complex ci(0.0,scale);
    Matrix ta;

    out=zero;
    for(int a=0;a<generators();a++){
      gaussian(pRNG,ca); 
      generator(a,ta);
      la=toComplex(ca)*ci*ta;
      out += la; 
    }

  }


  static void HotConfiguration(GridParallelRNG &pRNG,LatticeGaugeField &out){
    LatticeMatrix Umu(out._grid);
    for(int mu=0;mu<Nd;mu++){
      LieRandomize(pRNG,Umu,1.0);
      PokeIndex<LorentzIndex>(out,Umu,mu);
    }
  }
  static void TepidConfiguration(GridParallelRNG &pRNG,LatticeGaugeField &out){
    LatticeMatrix Umu(out._grid);
    for(int mu=0;mu<Nd;mu++){
      LieRandomize(pRNG,Umu,0.01);
      pokeLorentz(out,Umu,mu);
    }
  }
  static void ColdConfiguration(GridParallelRNG &pRNG,LatticeGaugeField &out){
    LatticeMatrix Umu(out._grid);
    Umu=1.0;
    for(int mu=0;mu<Nd;mu++){
      pokeLorentz(out,Umu,mu);
    }
  }

  static void taProj( const LatticeMatrix &in,  LatticeMatrix &out){
    out = Ta(in);
  }
  static void taExp( const LatticeMatrix &x,  LatticeMatrix &ex){ 

    LatticeMatrix xn(x._grid);
    RealD nfac = 1.0;

    xn = x;
    ex =xn+Complex(1.0); // 1+x

    // Do a 12th order exponentiation
    for(int i=2; i <= 12; ++i)
    {
      nfac = nfac/RealD(i); //1/2, 1/2.3 ...
      xn   = xn * x; // x2, x3,x4....
      ex   = ex+ xn*nfac;// x2/2!, x3/3!....
    }
  }

};

 typedef SU<2> SU2;
 typedef SU<3> SU3;
 typedef SU<4> SU4;
 typedef SU<5> SU5;

  }
}
#endif
