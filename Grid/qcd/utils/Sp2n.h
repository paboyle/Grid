
#ifndef QCD_UTIL_Sp2n_H
#define QCD_UTIL_Sp2n_H

NAMESPACE_BEGIN(Grid);

// Sp(2N)
// ncolour = N

// need to be careful with n and 2n
// I am defining the Sp class for Sp(2n) to be such that the template variable ncolour
// is n inside 2n and the typedef at the end of the file should eliminate possible confusion.

// the other routines, like projectOnSp2n, N will be the dimension of the actual number of colors for consistency with the sun routines

template <int ncolour>
class Sp {
public:
    static const int Dimension = ncolour*2;
    static const int AlgebraDimension = ncolour*(2*ncolour +1);
    static int su2subgroups(void) { return (ncolour * (ncolour - 1)) / 2; }
    
    
    template <typename vtype>
    using iSp2nMatrix = iScalar<iScalar<iMatrix<vtype, Dimension> > >;
    template <typename vtype>
    using iSU2Matrix = iScalar<iScalar<iMatrix<vtype, 2> > >;
    template <typename vtype>
    using iSp2nAlgebraVector = iScalar<iScalar<iVector<vtype, AlgebraDimension> > >;
    
    typedef iSp2nMatrix<Complex> Matrix;
    typedef iSp2nMatrix<ComplexF> MatrixF;
    typedef iSp2nMatrix<ComplexD> MatrixD;
    
    typedef iSp2nMatrix<vComplex> vMatrix;
    typedef iSp2nMatrix<vComplexF> vMatrixF;
    typedef iSp2nMatrix<vComplexD> vMatrixD;
    
    typedef iSp2nAlgebraVector<Complex> AlgebraVector;
    typedef iSp2nAlgebraVector<ComplexF> AlgebraVectorF;
    typedef iSp2nAlgebraVector<ComplexD> AlgebraVectorD;
    
    typedef iSp2nAlgebraVector<vComplex> vAlgebraVector;
    typedef iSp2nAlgebraVector<vComplexF> vAlgebraVectorF;
    typedef iSp2nAlgebraVector<vComplexD> vAlgebraVectorD;
    
    typedef Lattice<vMatrix> LatticeMatrix;
    typedef Lattice<vMatrixF> LatticeMatrixF;
    typedef Lattice<vMatrixD> LatticeMatrixD;
    
    typedef Lattice<vAlgebraVector> LatticeAlgebraVector;
    typedef Lattice<vAlgebraVectorF> LatticeAlgebraVectorF;
    typedef Lattice<vAlgebraVectorD> LatticeAlgebraVectorD;
    
    
    
    // Sp(2N) has N(2N+1) = 2N^2+N generators
    //
    // normalise the generators such that
    // Trace ( Ta Tb) = 1/2 delta_ab
    //
    // N generators in the cartan, 2N^2 off
    // off diagonal:
    //     there are 6 types named a,b,c,d and w,z
    //     abcd are N(N-1)/2 each while wz are N each
    
    
    
    
    template <class cplx>
    static void generator(int lieIndex, iSp2nMatrix<cplx> &ta) {
        // map lie index into type of generators: diagonal, abcd type, wz type

        int diagIndex;
        int aIndex, bIndex, cIndex, dIndex;
        int wIndex, zIndex; // a,b,c,d are N(N-1)/2 and w,z are N
        int mod = ncolour * (ncolour-1) * 0.5;
        int offdiag = 2*ncolour*ncolour; // number of generators not in the cartan subalgebra
        int wmod = 4*mod;
        int zmod = wmod+ncolour;
        if (lieIndex >= offdiag) {
            diagIndex = lieIndex - offdiag; // 0, ... ,N-1
            //std::cout << GridLogMessage << "diag type " << std::endl;
            generatorDiagtype(diagIndex, ta);
            return;
        }
        if ( (lieIndex >= wmod) && (lieIndex < zmod) ) {
            //std::cout << GridLogMessage << "w type " << std::endl;
            wIndex = lieIndex- wmod;    // 0, ... ,N-1
            generatorWtype(wIndex,ta);
            return;
        }
        if ( (lieIndex >= zmod) && (lieIndex < offdiag) ) {
            //std::cout << GridLogMessage << "z type " << std::endl;
            //std::cout << GridLogMessage << "lie index " << lieIndex << std::endl;
            //std::cout << GridLogMessage << "z mod " << zmod << std::endl;
            zIndex = lieIndex - zmod;   // 0, ... ,N-1
            generatorZtype(zIndex,ta);
            return;
        }
        if (lieIndex < mod) { // atype 0, ... , N(N-1)/2=mod
            //std::cout << GridLogMessage << "a type " << std::endl;
            aIndex = lieIndex;
            //std::cout << GridLogMessage << "a indx " << aIndex << std::endl;
            generatorAtype(aIndex, ta);
            return;
        }
        if ( (lieIndex >= mod) && lieIndex < 2*mod) { // btype mod, ... , 2mod-1
            //std::cout << GridLogMessage << "b type " << std::endl;
            bIndex = lieIndex - mod;
            generatorBtype(bIndex, ta);
            return;
        }
        if ( (lieIndex >= 2*mod) && lieIndex < 3*mod) { // ctype 2mod, ... , 3mod-1
            //std::cout << GridLogMessage << "c type " << std::endl;
            cIndex = lieIndex - 2*mod;
            generatorCtype(cIndex, ta);
            return;
        }
        if ( (lieIndex >= 3*mod) && lieIndex < wmod) { // ctype 3mod, ... , 4mod-1 = wmod-1
             //std::cout << GridLogMessage << "d type " << std::endl;
             dIndex = lieIndex - 3*mod;
             generatorDtype(dIndex, ta);
            return;
         }
            
    } //end of generator
    
    template <class cplx>
    static void generatorDiagtype(int diagIndex, iSp2nMatrix<cplx> &ta) {
        
        // ta(i,i) = - ta(i+N,i+N) = 1/2 for each i index of the cartan subalgebra
        
        ta = Zero();
        RealD nrm = 1.0 / 2;

        ta()()(diagIndex,diagIndex) = nrm;
        ta()()(diagIndex+ncolour,diagIndex+ncolour) = -nrm;
    }
    
    template <class cplx>
    static void generatorAtype(int aIndex, iSp2nMatrix<cplx> &ta) {
        
        // ta(i,j) = ta(j,i) = -ta(i+N,j+N) = -ta(j+N,i+N) = 1 / 2 sqrt(2)
        // with i<j and i=0,...,N-2
        // follows that j=i+1, ... , N
        int i1, i2;
        ta = Zero();
        RealD nrm = 1 / (2 * std::sqrt(2) );
        
        su2SubGroupIndex(i1, i2, aIndex);
        ta()()(i1,i2) = 1;
        ta()()(i2,i1) = 1;
        ta()()(i1+ncolour,i2+ncolour) = -1;
        ta()()(i2+ncolour,i1+ncolour) = -1;

        ta = ta * nrm;
    }
    
    template <class cplx>
    static void generatorBtype(int bIndex, iSp2nMatrix<cplx> &ta) {
        
        // ta(i,j) = -ta(j,i) = ta(i+N,j+N) = -ta(j+N,i+N) = i / 1/ 2 sqrt(2)
        // with i<j and i=0,...,N-2
        // follows that j=i+1, ... , N-1
        
        int i1, i2;
        ta = Zero();
        cplx i(0.0, 1.0);
        RealD nrm = 1 / (2 * std::sqrt(2));
        su2SubGroupIndex(i1, i2, bIndex);
        

        ta()()(i1,i2) = i;
        ta()()(i2,i1) = -i;
        ta()()(i1+ncolour,i2+ncolour) = i;
        ta()()(i2+ncolour,i1+ncolour) = -i;

        ta = ta * nrm;
    }
    
    template <class cplx>
    static void generatorCtype(int cIndex, iSp2nMatrix<cplx> &ta) {
        
        // ta(i,j+N) = ta(j,i+N) = ta(i+N,j) = ta(j+N,i) = 1 / 2 sqrt(2)
        
        
        int i1, i2;
        ta = Zero();
        RealD nrm = 1 / (2 * std::sqrt(2) );
        su2SubGroupIndex(i1, i2, cIndex);
        
        ta()()(i1,i2+ncolour) = 1;
        ta()()(i2,i1+ncolour) = 1;
        ta()()(i1+ncolour,i2) = 1;
        ta()()(i2+ncolour,i1) = 1;

        ta = ta * nrm;
    }
    
    template <class cplx>
    static void generatorDtype(int dIndex, iSp2nMatrix<cplx> &ta) {
        
        // ta(i,j+N) = ta(j,i+N) = -ta(i+N,j) = -ta(j+N,i) = i /  2 sqrt(2)
        
        int i1, i2;
        ta = Zero();
        cplx i(0.0, 1.0);
        RealD nrm = 1 / (2 * std::sqrt(2)  );
        su2SubGroupIndex(i1, i2, dIndex);

        ta()()(i1,i2+ncolour) = i;
        ta()()(i2,i1+ncolour) = i;
        ta()()(i1+ncolour,i2) = -i;
        ta()()(i2+ncolour,i1) = -i;

        ta = ta * nrm;
    }
    
    template <class cplx>
     static void generatorWtype(int wIndex, iSp2nMatrix<cplx> &ta) {
         
         // ta(i,i+N) =  ta(i+N,i) = 1/2
         
         ta = Zero();
         RealD nrm = 1.0 / 2; //check

         ta()()(wIndex,wIndex+ncolour) = 1;
         ta()()(wIndex+ncolour,wIndex) = 1;
         
         ta = ta * nrm;
     }
    
    template <class cplx>
     static void generatorZtype(int zIndex, iSp2nMatrix<cplx> &ta) {
         
         // ta(i,i+N) = - ta(i+N,i) = i/2
         
         ta = Zero();
         RealD nrm = 1.0 / 2; //check
         cplx i(0.0, 1.0);
         ta()()(zIndex,zIndex+ncolour) = i;
         ta()()(zIndex+ncolour,zIndex) = -i;
         
         ta = ta * nrm;
     }
    
    
    ////////////////////////////////////////////////////////////////////////
    // Map a su2 subgroup number to the pair of rows that are non zero
    ////////////////////////////////////////////////////////////////////////
    static void su2SubGroupIndex(int &i1, int &i2, int su2_index) {
      assert((su2_index >= 0) && (su2_index < (ncolour * (ncolour - 1)) / 2));

      int spare = su2_index;
      for (i1 = 0; spare >= (ncolour - 1 - i1); i1++) {
        spare = spare - (ncolour - 1 - i1);  // remove the Nc-1-i1 terms
      }
      i2 = i1 + 1 + spare;
    }
    
    
    
    
    
    static void printGenerators(void) {
        for (int gen = 0; gen < AlgebraDimension; gen++) {
            Matrix ta;
            generator(gen, ta);
            std::cout << GridLogMessage << "Nc (2n) = " << 2*ncolour << std::endl;
            std::cout << GridLogMessage << " t_" << gen << std::endl;
            std::cout << GridLogMessage << ta << std::endl;
        }
    }
    
    
    
    static void testGenerators(void) {
        Matrix ta;
        Matrix tb;
        std::cout << GridLogMessage << "Fundamental - Checking trace ta tb is 0.5 delta_ab " << std::endl;
        for (int a = 0; a < AlgebraDimension; a++) {
            for (int b = 0; b < AlgebraDimension; b++) {
                generator(a,ta);
                generator(b,tb);
                Complex tr = TensorRemove(trace( ta * tb) );
                std::cout << GridLogMessage << "(" << a << "," << b << ") =  " << tr
                << std::endl;
                if (a == b) assert(abs(tr - Complex(0.5)) < 1.0e-6);
                if (a != b) assert(abs(tr) < 1.0e-6);

            }
        }
        std::cout << GridLogMessage << std::endl;
        std::cout << GridLogMessage << "Fundamental - Checking if hermitian" << std::endl;
        for (int a = 0; a < AlgebraDimension; a++) {
            generator(a,ta);
            std::cout << GridLogMessage << a << std::endl;
            assert(norm2(ta - adj(ta)) < 1.0e-6);
        }
        std::cout << GridLogMessage << std::endl;
        std::cout << GridLogMessage << "Fundamental - Checking if traceless" << std::endl;
        for (int a = 0; a < AlgebraDimension; a++) {
          generator(a, ta);
          Complex tr = TensorRemove(trace(ta));
          std::cout << GridLogMessage << a << std::endl;
          assert(abs(tr) < 1.0e-6);
        }
        
    }
    
    
    static void GaussianFundamentalLieAlgebraMatrix(GridParallelRNG &pRNG, //same as sun
                                                    LatticeMatrix &out,
                                                    Real scale = 1.0) {
      GridBase *grid = out.Grid();
      LatticeReal ca(grid);
      LatticeMatrix la(grid);
      Complex ci(0.0, scale);
      Matrix ta;

      out = Zero();
      for (int a = 0; a < AlgebraDimension; a++) {
        gaussian(pRNG, ca);
        generator(a, ta);
        la = toComplex(ca) * ta;
        out += la;
      }
      out *= ci;
    }
    
    
    
    template <typename LatticeMatrixType>
    static void taExp(const LatticeMatrixType &x, LatticeMatrixType &ex) { // same as sun
      typedef typename LatticeMatrixType::scalar_type ComplexType;

      LatticeMatrixType xn(x.Grid());
      RealD nfac = 1.0;

      xn = x;
      ex = xn + ComplexType(1.0);  // 1+x

      // Do a 12th order exponentiation
      for (int i = 2; i <= 12; ++i) {
        nfac = nfac / RealD(i);  // 1/2, 1/2.3 ...
        xn = xn * x;             // x2, x3,x4....
        ex = ex + xn * nfac;     // x2/2!, x3/3!....
      }
    }
    
 
    
    
/*    template <typename GaugeField>
    static void HotSpConfiguration(GridParallelRNG &pRNG, GaugeField &out) {
      typedef typename GaugeField::vector_type vector_type;
      typedef iSp2nMatrix<vector_type> vMatrixType;
      typedef Lattice<vMatrixType> LatticeMatrixType;

      LatticeMatrixType Umu(out.Grid());
      for (int mu = 0; mu < Nd; mu++) {
        LieRandomize(pRNG, Umu, 1.0);   //def
        PokeIndex<LorentzIndex>(out, Umu, mu);
      }
    }
    template<typename GaugeField>
    static void TepidSpConfiguration(GridParallelRNG &pRNG,GaugeField &out){
      typedef typename GaugeField::vector_type vector_type;
      typedef iSp2nMatrix<vector_type> vMatrixType;
      typedef Lattice<vMatrixType> LatticeMatrixType;

      LatticeMatrixType Umu(out.Grid());
      for(int mu=0;mu<Nd;mu++){
        LieRandomize(pRNG,Umu,0.01);    //def
        PokeIndex<LorentzIndex>(out,Umu,mu);
      }
    } */
    
    template<typename GaugeField>
    static void ColdConfiguration(GaugeField &out){
      typedef typename GaugeField::vector_type vector_type;
      typedef iSp2nMatrix<vector_type> vMatrixType;
      typedef Lattice<vMatrixType> LatticeMatrixType;

      LatticeMatrixType Umu(out.Grid());
      Umu=1.0;
      for(int mu=0;mu<Nd;mu++){
        PokeIndex<LorentzIndex>(out,Umu,mu);
      }
    }
    template<typename GaugeField>
    static void ColdConfiguration(GridParallelRNG &pRNG,GaugeField &out){
      ColdConfiguration(out);
    }
    
    
    
}; // end of class Sp
    
    


template<int N> //same of su
 LatticeComplexD SpDeterminant(const Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > &Umu)
 {
   GridBase *grid=Umu.Grid();
   auto lvol = grid->lSites();
   LatticeComplexD ret(grid);

   autoView(Umu_v,Umu,CpuRead);
   autoView(ret_v,ret,CpuWrite);
   thread_for(site,lvol,{
     Eigen::MatrixXcd EigenU = Eigen::MatrixXcd::Zero(N,N);
     Coordinate lcoor;
     grid->LocalIndexToLocalCoor(site, lcoor);
     iScalar<iScalar<iMatrix<ComplexD, N> > > Us;
     peekLocalSite(Us, Umu_v, lcoor);
     for(int i=0;i<N;i++){
       for(int j=0;j<N;j++){
     EigenU(i,j) = Us()()(i,j);
       }}
     ComplexD det = EigenU.determinant();
     pokeLocalSite(det,ret_v,lcoor);
   });
   return ret;
 }
 
 template<int N>
 static void ProjectSp2n(Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > &Umu)
 {
   Umu      = ProjectOnSpGroup(Umu);
   auto det = SpDeterminant(Umu); // ok ?

   det = conjugate(det);

   for(int i=0;i<N;i++){
     auto element = PeekIndex<ColourIndex>(Umu,N-1,i);
     element = element * det;
     PokeIndex<ColourIndex>(Umu,element,Nc-1,i);
   }
 }
 template<int N>
 static void ProjectSp2n(Lattice<iVector<iScalar<iMatrix<vComplexD, N> >,Nd> > &U)
 {
   GridBase *grid=U.Grid();
   for(int mu=0;mu<Nd;mu++){
     auto Umu = PeekIndex<LorentzIndex>(U,mu);
     Umu      = ProjectOnSpGroup(Umu);
     ProjectSp2n(Umu);
     PokeIndex<LorentzIndex>(U,Umu,mu);
   }
 }

typedef Sp<1> Sp2;
typedef Sp<2> Sp4;
typedef Sp<3> Sp6;
typedef Sp<4> Sp8;

NAMESPACE_END(Grid);
#endif
