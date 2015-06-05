#ifndef G5_HERMITIAN_LINOP
#define G5_HERMITIAN_LINOP
namespace Grid {
  namespace QCD {
////////////////////////////////////////////////////////////////////
// Wrap an already herm matrix
////////////////////////////////////////////////////////////////////
template<class Matrix,class Field>
class Gamma5HermitianLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
public:
  Gamma5HermitianLinearOperator(Matrix &Mat): _Mat(Mat){};
  void Op     (const Field &in, Field &out){
    _Mat.M(in,out);
  }
  void AdjOp     (const Field &in, Field &out){
    _Mat.M(in,out);
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){

    HermOp(in,out);
    
    ComplexD dot;
    dot= innerProduct(in,out);
    n1=real(dot);
    
    dot = innerProduct(out,out);
    n2=real(dot);
  }
  void HermOp(const Field &in, Field &out){
    Field tmp(in._grid);
    _Mat.M(in,tmp);
    G5R5(out,tmp);
  }
};

}}
#endif
