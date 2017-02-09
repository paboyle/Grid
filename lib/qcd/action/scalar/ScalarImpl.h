#ifndef SCALAR_IMPL
#define SCALAR_IMPL


namespace Grid {
  //namespace QCD {

  template <class S>
  class ScalarImplTypes {
  public:
    typedef S Simd;
    
    template <typename vtype>
    using iImplField = iScalar<iScalar<iScalar<vtype> > >;
    
    typedef iImplField<Simd> SiteField;
    
    
    typedef Lattice<SiteField> Field;
    
    static inline void generate_momenta(Field& P, GridParallelRNG& pRNG){
      gaussian(pRNG, P);
    }
    
    static inline Field projectForce(Field& P){return P;}
    
    static inline void update_field(Field& P, Field& U, double ep){
      U += P*ep;
    }
    
    static inline RealD FieldSquareNorm(Field& U){
      return (- sum(trace(U*U))/2.0);
    }
    
    static inline void HotConfiguration(GridParallelRNG &pRNG, Field &U) {
      gaussian(pRNG, U);
    }
    
    static inline void TepidConfiguration(GridParallelRNG &pRNG, Field &U) {
      gaussian(pRNG, U);
    }
    
    static inline void ColdConfiguration(GridParallelRNG &pRNG, Field &U) {
      U = 1.0;
    }
    
  };

  template <class S, unsigned int N>
  class ScalarMatrixImplTypes {
  public:
    typedef S Simd;
    
    template <typename vtype>
    using iImplField = iScalar<iScalar<iMatrix<vtype, N> > >;
    
    typedef iImplField<Simd> SiteField;
    
    
    typedef Lattice<SiteField> Field;
    
    static inline void generate_momenta(Field& P, GridParallelRNG& pRNG){
      gaussian(pRNG, P);
    }
    
    static inline Field projectForce(Field& P){return P;}
    
    static inline void update_field(Field& P, Field& U, double ep){
      U += P*ep;
    }
    
    static inline RealD FieldSquareNorm(Field& U){
      return (TensorRemove(- sum(trace(U*U))*0.5).real());
    }
    
    static inline void HotConfiguration(GridParallelRNG &pRNG, Field &U) {
      gaussian(pRNG, U);
    }
    
    static inline void TepidConfiguration(GridParallelRNG &pRNG, Field &U) {
      gaussian(pRNG, U);
    }
    
    static inline void ColdConfiguration(GridParallelRNG &pRNG, Field &U) {
      U = 1.0;
    }
    
  };


  
  
  typedef ScalarImplTypes<vReal> ScalarImplR;
  typedef ScalarImplTypes<vRealF> ScalarImplF;
  typedef ScalarImplTypes<vRealD> ScalarImplD;
  
  //} 
} 

#endif
