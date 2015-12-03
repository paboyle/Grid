#ifndef G5_HERMITIAN_LINOP
#define G5_HERMITIAN_LINOP

namespace Grid {
  namespace QCD {

////////////////////////////////////////////////////////////////////
// Wrap an already herm matrix
////////////////////////////////////////////////////////////////////
template<class Matrix,class Field>
class Gamma5R5HermitianLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
public:
  Gamma5R5HermitianLinearOperator(Matrix &Mat): _Mat(Mat){};
  void Op     (const Field &in, Field &out){
    HermOp(in,out);
  }
  void AdjOp     (const Field &in, Field &out){
    HermOp(in,out);
  }
  void OpDiag (const Field &in, Field &out) {
    Field tmp(in._grid);
    _Mat.Mdiag(in,tmp);
    G5R5(out,tmp);
  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {
    Field tmp(in._grid);
    _Mat.Mdir(in,tmp,dir,disp);
    G5R5(out,tmp);
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


template<class Matrix,class Field>
class Gamma5HermitianLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  Gamma g5;
public:
    Gamma5HermitianLinearOperator(Matrix &Mat): _Mat(Mat), g5(Gamma::Gamma5) {};
  void Op     (const Field &in, Field &out){
    HermOp(in,out);
  }
  void AdjOp     (const Field &in, Field &out){
    HermOp(in,out);
  }
  void OpDiag (const Field &in, Field &out) {
    Field tmp(in._grid);
    _Mat.Mdiag(in,tmp);
    out=g5*tmp;
  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {
    Field tmp(in._grid);
    _Mat.Mdir(in,tmp,dir,disp);
    out=g5*tmp;
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
    out=g5*tmp;
  }
};


}}
#endif
