#ifndef GRID_SCHUR_RED_BLACK_H
#define GRID_SCHUR_RED_BLACK_H

  /*
   * Red black Schur decomposition
   *
   *  M = (Mee Meo) =  (1             0 )   (Mee   0               )  (1 Mee^{-1} Meo)
   *      (Moe Moo)    (Moe Mee^-1    1 )   (0   Moo-Moe Mee^-1 Meo)  (0   1         )
   *                =         L                     D                     U
   *
   * L^-1 = (1              0 )
   *        (-MoeMee^{-1}   1 )   
   * L^{dag} = ( 1       Mee^{-dag} Moe^{dag} )
   *           ( 0       1                    )
   * L^{-d}  = ( 1      -Mee^{-dag} Moe^{dag} )
   *           ( 0       1                    )
   *
   * U^-1 = (1   -Mee^{-1} Meo)
   *        (0    1           )
   * U^{dag} = ( 1                 0)
   *           (Meo^dag Mee^{-dag} 1)
   * U^{-dag} = (  1                 0)
   *            (-Meo^dag Mee^{-dag} 1)
   ***********************
   *     M psi = eta
   ***********************
   *Odd
   * i)   (D_oo)^{\dag} D_oo psi_o = (D_oo)^\dag L^{-1}  eta_o
   *                        eta_o' = D_oo (eta_o - Moe Mee^{-1} eta_e)
   *Even
   * ii)  Mee psi_e + Meo psi_o = src_e
   *
   *   => sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
   *
   */

namespace Grid {

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Take a matrix and form a Red Black solver calling a Herm solver
  // Use of RB info prevents making SchurRedBlackSolve conform to standard interface
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class Field> class SchurRedBlackSolve : public OperatorFunction<Field>{
  private:
    SparseMatrixBase<Field> & _Matrix;
    OperatorFunction<Field> & _HermitianRBSolver;
    int CBfactorise;
  public:

    /////////////////////////////////////////////////////
    // Wrap the usual normal equations Schur trick
    /////////////////////////////////////////////////////
  SchurRedBlackSolve(SparseMatrixBase<Field> &Matrix, OperatorFunction<Field> &HermitianRBSolver) 
    :  _Matrix(Matrix), _HermitianRBSolver(HermitianRBSolver) { 
      CBfactorise=0;
    };

    void operator() (const Field &in, Field &out){

      // FIXME CGdiagonalMee not implemented virtual function
      // FIXME need to make eo grid from full grid.
      // FIXME use CBfactorise to control schur decomp
      const int Even=0;
      const int Odd =1;

      // Make a cartesianRedBlack from full Grid
      GridRedBlackCartesian grid(in._grid);
 
      Field src_e(&grid);
      Field src_o(&grid);
      Field sol_e(&grid);
      Field sol_o(&grid);
      Field   tmp(&grid);
      Field  Mtmp(&grid);
      
      pickCheckerboard(Even,src_e,in);
      pickCheckerboard(Odd ,src_o,in);

      /////////////////////////////////////////////////////
      // src_o = Mdag * (source_o - Moe MeeInv source_e)
      /////////////////////////////////////////////////////
      _Matrix.MooeeInv(src_e,tmp);     //    MooeeInv(source[Even],tmp,DaggerNo,Even);
      _Matrix.Meooe   (tmp,Mtmp);      //    Meo     (tmp,src,Odd,DaggerNo);
      tmp=src_o-Mtmp;                  //    axpy    (tmp,src,source[Odd],-1.0);
      _Matrix.MpcDag(tmp,src_o);       //    Mprec(tmp,src,Mtmp,DaggerYes);  

      //////////////////////////////////////////////////////////////
      // Call the red-black solver
      //////////////////////////////////////////////////////////////
      _HermitianRBSolver(src_o,sol_o); //  CGNE_prec_MdagM(solution[Odd],src);

      ///////////////////////////////////////////////////
      // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
      ///////////////////////////////////////////////////
      _Matrix.Meooe(sol_o,tmp);        // Meo(solution[Odd],tmp,Even,DaggerNo);
      src_e = src_e-tmp;               // axpy(src,tmp,source[Even],-1.0);
      _Matrix.MooeeInv(src_e,sol_e);   // MooeeInv(src,solution[Even],DaggerNo,Even);
     
      setCheckerboard(out,sol_e);
      setCheckerboard(out,sol_o);
 
    }     
  };

}
#endif
