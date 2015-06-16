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
  template<class vreal,class vcplx>
  static void su2Extract(std::vector<Lattice<iSinglet  <vreal> > > &r,
			       const Lattice<iSUnMatrix<vcplx> > &source, 
			 int su2_index)
  {
    GridBase *grid(source._grid);

    assert(r.size() == 4); // store in 4 real parts
    
    for(int i=0;i<4;i++){
      conformable(r[i],source);
    }
    
    int i1,i2;    
    su2SubGroupIndex(i1,i2,su2_index);

    /* Compute the b(k) of A_SU(2) = b0 + i sum_k bk sigma_k */ 
    
    //    r[0] = toReal(real(peekColour(source,i1,i1)) + real(peekColour(source,i2,i2)));
    //    r[1] = toReal(imag(peekColour(source,i1,i2)) + imag(peekColour(source,i2,i1)));
    //    r[2] = toReal(real(peekColour(source,i1,i2)) - real(peekColour(source,i2,i1)));
    //    r[3] = toReal(imag(peekColour(source,i1,i1)) - imag(peekColour(source,i2,i2)));
    r[0] = toReal(real(peekColour(source,i1,i1)) + real(peekColour(source,i2,i2)));
    r[1] = toReal(imag(peekColour(source,i1,i2)) + imag(peekColour(source,i2,i1)));
    r[2] = toReal(real(peekColour(source,i1,i2)) - real(peekColour(source,i2,i1)));
    r[3] = toReal(imag(peekColour(source,i1,i1)) - imag(peekColour(source,i2,i2)));
    
  }

  template<class vreal,class vcplx>
  static void su2Insert(const std::vector<Lattice<iSinglet<vreal> > > &r,
			Lattice<iSUnMatrix<vcplx> > &dest,
			int su2_index)
  {
    typedef typename Lattice<iSUnMatrix<vcplx> >::scalar_type cplx;
    typedef Lattice<iSinglet<vcplx> > Lcomplex;
    GridBase * grid = dest._grid;

    assert(r.size() == 4); // store in 4 real parts

    Lcomplex tmp(grid);

    std::vector<Lcomplex > cr(4,grid);
    for(int i=0;i<r.size();i++){
      conformable(r[i],dest);
      cr[i]=toComplex(r[i]);
    }

    int i1,i2;
    su2SubGroupIndex(i1,i2,su2_index);

    cplx one   (1.0,0.0);
    cplx cplx_i(0.0,1.0);

    tmp =   cr[0]*one+ cr[3]*cplx_i;   pokeColour(dest,tmp,i1,i2); 
    tmp =   cr[2]*one+ cr[1]*cplx_i;   pokeColour(dest,tmp,i1,i2);
    tmp =  -cr[2]*one+ cr[1]*cplx_i;   pokeColour(dest,tmp,i2,i1);
    tmp =   cr[0]*one- cr[3]*cplx_i;   pokeColour(dest,tmp,i2,i2);
  }

  static void SubGroupHeatBath(      GridSerialRNG       &sRNG,
				     GridParallelRNG     &pRNG,
			             RealD beta,
				     LatticeMatrix &link, 
			       const LatticeMatrix &staple, 
				     int su2_subgroup,
				     int nheatbath,
				     int& ntrials,
				     int& nfails,
				     LatticeInteger &wheremask)
  {
    GridBase *grid = link._grid;

    LatticeMatrix V(grid);
    V = link*staple;

    std::vector<LatticeReal> r(4,grid);
    std::vector<LatticeReal> a(4,grid);
    su2Extract(r,V,su2_subgroup); // HERE

    LatticeReal r_l(grid);
    r_l  = r[0]*r[0]+r[1]*r[1]+r[2]*r[2]+r[3]*r[3];
    r_l = sqrt(r_l);

    LatticeReal ftmp(grid);
    LatticeReal ftmp1(grid);
    LatticeReal ftmp2(grid);
    LatticeReal one (grid);    one = 1.0;
    LatticeReal zz  (grid);    zz  = zero;
    LatticeReal recip(grid);   recip = one/r_l;
    
    Real machine_epsilon= 1.0e-14;
    
    ftmp = where(r_l>machine_epsilon,recip,one);
    a[0] = where(r_l>machine_epsilon,   r[0] * ftmp , one);
    a[1] = where(r_l>machine_epsilon, -(r[1] * ftmp), zz);
    a[2] = where(r_l>machine_epsilon, -(r[2] * ftmp), zz);
    a[3] = where(r_l>machine_epsilon, -(r[3] * ftmp), zz);

    r_l *=  beta / ncolour;
    
    ftmp1 = where(wheremask,one,zz);
    Real num_sites = TensorRemove(sum(ftmp1));

    Integer itrials = (Integer)num_sites;
    ntrials = 0;
    nfails = 0;

    LatticeInteger lupdate(grid);
    LatticeInteger lbtmp(grid);
    LatticeInteger lbtmp2(grid); lbtmp2=zero;

    int n_done = 0;
    int nhb = 0;

    r[0] = a[0];
    lupdate = 1;

    LatticeReal ones (grid); ones = 1.0;
    LatticeReal zeros(grid); zeros=zero;

    const RealD twopi=2.0*M_PI;
    while ( nhb < nheatbath && n_done < num_sites ) {

       ntrials += itrials;

       random(pRNG,r[1]);
       std::cout<<"RANDOM SPARSE FLOAT r[1]"<<std::endl;
       std::cout<<r[1]<<std::endl;

       random(pRNG,r[2]);
       random(pRNG,ftmp);

       r[1] = log(r[1]);
       r[2] = log(r[2]);

       ftmp = ftmp*twopi;
       r[3] = cos(ftmp);
       r[3] = r[3]*r[3];

       r[1] += r[2] * r[3];
       r[2]  = r[1] / r_l;

       random(pRNG,ftmp);
       r[1] = ftmp*ftmp;

       {
	 LatticeInteger mask_true (grid); mask_true = 1;
	 LatticeInteger mask_false(grid); mask_false= 0;
	 LatticeReal    thresh(grid); thresh = (1.0 + 0.5*r[2]);
	 lbtmp = where(r[1] <= thresh,mask_true,mask_false);
       }
       lbtmp2= lbtmp && lupdate;       
       r[0]  = where(lbtmp2, 1.0+r[2], r[0]);

       ftmp1 = where(lbtmp2,ones,zeros);
       RealD sitesum = sum(ftmp1);
       Integer itmp = sitesum;

       n_done += itmp;
       itrials -= itmp;
       nfails += itrials;
       lbtmp   = !lbtmp;
       lupdate = lupdate & lbtmp;
       ++nhb;
    }
    // Now create r[1], r[2] and r[3] according to the spherical measure 
    // Take absolute value to guard against round-off 
    random(pRNG,ftmp1);
    r[2] = 1.0 - 2.0*ftmp1;
 
    ftmp1 = abs(1.0 - r[0]*r[0]);
    r[3]   = -(sqrt(ftmp1) * r[2]);

    // Take absolute value to guard against round-off 
    r_l    = sqrt(abs(ftmp1 - r[3]*r[3]));
 
    random(pRNG,ftmp1);
    ftmp1 *= twopi;
    r[1]    = r_l * cos(ftmp1);
    r[2]    = r_l * sin(ftmp1);

    // Update matrix is B = R * A, with B, R and A given by b_i, r_i and a_i 
    std::vector<LatticeReal> b(4,grid);
    b[0] = r[0]*a[0] - r[1]*a[1] - r[2]*a[2] - r[3]*a[3];
    b[1] = r[0]*a[1] + r[1]*a[0] - r[2]*a[3] + r[3]*a[2];
    b[2] = r[0]*a[2] + r[2]*a[0] - r[3]*a[1] + r[1]*a[3];
    b[3] = r[0]*a[3] + r[3]*a[0] - r[1]*a[2] + r[2]*a[1];

     //
     // Now fill an SU(3) matrix V with the SU(2) submatrix su2_index
     // parametrized by a_k in the sigma matrix basis.
     //
    su2Insert(b,V,su2_subgroup);

    // U = V*U
    LatticeMatrix tmp(grid);
    tmp = V * link;

    //mask the assignment back
    link = where(wheremask,tmp,link);

  }

  static void printGenerators(void)
  {
    for(int gen=0;gen<generators();gen++){
      Matrix ta;
      generator(gen,ta);
      std::cout<< "Nc = "<<ncolour<<" t_"<<gen<<std::endl;
      std::cout<<ta<<std::endl;
    }
  }

  static void testGenerators(void){
    Matrix ta;
    Matrix tb;
    std::cout<<"Checking trace ta tb is 0.5 delta_ab"<<std::endl;
    for(int a=0;a<generators();a++){
      for(int b=0;b<generators();b++){
	generator(a,ta);
	generator(b,tb);
	Complex tr =TensorRemove(trace(ta*tb)); 
	std::cout<<tr<<" ";
	if(a==b) assert(abs(tr-Complex(0.5))<1.0e-6);
	if(a!=b) assert(abs(tr)<1.0e-6);
      }
      std::cout<<std::endl;
    }
    std::cout<<"Checking hermitian"<<std::endl;
    for(int a=0;a<generators();a++){
      generator(a,ta);
      std::cout<<a<<" ";
      assert(norm2(ta-adj(ta))<1.0e-6);
    }    
    std::cout<<std::endl;

    std::cout<<"Checking traceless"<<std::endl;
    for(int a=0;a<generators();a++){
      generator(a,ta);
      Complex tr =TensorRemove(trace(ta)); 
      std::cout<<a<<" ";
      assert(abs(tr)<1.0e-6);
    }    
    std::cout<<std::endl;
  }

  // reunitarise??

  static void taProj( const LatticeMatrix &in,  LatticeMatrix &out){
    out = Ta(in);
  }
  static void taExp( const LatticeMatrix &x,  LatticeMatrix &ex){ 
    LatticeMatrix xn   = x;

    RealD nfac = 1.0;
    ex = 1+x; // 1+x

    // Do a 12th order exponentiation
    for(int i= 2; i <= 12; ++i)
    {
      nfac = nfac/i; 
      xn   = xn * x; // x2, x3,x4....
      ex  += xn*nfac;// x2/2!, x3/3!....
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
