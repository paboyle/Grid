#ifndef  GRID_ALGORITHM_SPARSE_MATRIX_H
#define  GRID_ALGORITHM_SPARSE_MATRIX_H

#include <Grid.h>

namespace Grid {

  /////////////////////////////////////////////////////////////////////////////////////////////
  // Interface defining what I expect of a general sparse matrix, such as a Fermion action
  /////////////////////////////////////////////////////////////////////////////////////////////
    template<class Field> class SparseMatrixBase {
    public:
      virtual GridBase *Grid(void) =0;
      // Full checkerboar operations
      virtual RealD M    (const Field &in, Field &out)=0;
      virtual RealD Mdag (const Field &in, Field &out)=0;
      virtual void  MdagM(const Field &in, Field &out,RealD &ni,RealD &no) {
	Field tmp (in._grid);
	ni=M(in,tmp);
	no=Mdag(tmp,out);
      }
      virtual  void Mdiag    (const Field &in, Field &out)=0;
      virtual  void Mdir     (const Field &in, Field &out,int dir, int disp)=0;
    };

  /////////////////////////////////////////////////////////////////////////////////////////////
  // Interface augmented by a red black sparse matrix, such as a Fermion action
  /////////////////////////////////////////////////////////////////////////////////////////////
    template<class Field> class CheckerBoardedSparseMatrixBase : public SparseMatrixBase<Field> {
    public:
      virtual GridBase *RedBlackGrid(void)=0;
      // half checkerboard operaions
      virtual  void Meooe    (const Field &in, Field &out)=0;
      virtual  void Mooee    (const Field &in, Field &out)=0;
      virtual  void MooeeInv (const Field &in, Field &out)=0;

      virtual  void MeooeDag    (const Field &in, Field &out)=0;
      virtual  void MooeeDag    (const Field &in, Field &out)=0;
      virtual  void MooeeInvDag (const Field &in, Field &out)=0;

    };

}

#endif
