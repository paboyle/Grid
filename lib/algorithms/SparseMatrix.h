#ifndef  GRID_ALGORITHM_SPARSE_MATRIX_H
#define  GRID_ALGORITHM_SPARSE_MATRIX_H

#include <Grid.h>

namespace Grid {

  /////////////////////////////////////////////////////////////////////////////////////////////
  // Interface defining what I expect of a general sparse matrix, such as a Fermion action
  /////////////////////////////////////////////////////////////////////////////////////////////
    template<class Field> class SparseMatrixBase {
    public:
      GridBase *_grid;
      // Full checkerboar operations
      virtual RealD M    (const Field &in, Field &out)=0;
      virtual RealD Mdag (const Field &in, Field &out)=0;
      virtual void  MdagM(const Field &in, Field &out,RealD &ni,RealD &no) {
	Field tmp (in._grid);
	ni=M(in,tmp);
	no=Mdag(tmp,out);
      }
      SparseMatrixBase(GridBase *grid) : _grid(grid) {};
    };

  /////////////////////////////////////////////////////////////////////////////////////////////
  // Interface augmented by a red black sparse matrix, such as a Fermion action
  /////////////////////////////////////////////////////////////////////////////////////////////
    template<class Field> class CheckerBoardedSparseMatrixBase : public SparseMatrixBase<Field> {
    public:
      GridBase *_cbgrid;
      // half checkerboard operaions
      virtual  void Meooe    (const Field &in, Field &out)=0;
      virtual  void Mooee    (const Field &in, Field &out)=0;
      virtual  void MooeeInv (const Field &in, Field &out)=0;

      virtual  void MeooeDag    (const Field &in, Field &out)=0;
      virtual  void MooeeDag    (const Field &in, Field &out)=0;
      virtual  void MooeeInvDag (const Field &in, Field &out)=0;

      // Schur decomp operators
      virtual  RealD Mpc      (const Field &in, Field &out) {
	Field tmp(in._grid);

	Meooe(in,tmp);
	MooeeInv(tmp,out);
	Meooe(out,tmp);

	Mooee(in,out);
	return axpy_norm(out,-1.0,tmp,out);
      }
      virtual  RealD MpcDag   (const Field &in, Field &out){
	Field tmp(in._grid);

	MeooeDag(in,tmp);
	MooeeInvDag(tmp,out);
	MeooeDag(out,tmp);

	MooeeDag(in,out);
	return axpy_norm(out,-1.0,tmp,out);
      }
      virtual void MpcDagMpc(const Field &in, Field &out,RealD &ni,RealD &no) {
	Field tmp(in._grid);
	ni=Mpc(in,tmp);
	no=MpcDag(tmp,out);
	//	std::cout<<"MpcDagMpc "<<ni<<" "<<no<<std::endl;
      }
      CheckerBoardedSparseMatrixBase(GridBase *grid,GridBase *cbgrid) : SparseMatrixBase<Field>(grid), _cbgrid(cbgrid) {};
    };

}

#endif
