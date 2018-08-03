    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/SchurRedBlack.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
   * i)                 D_oo psi_o =  L^{-1}  eta_o
   *                        eta_o' = (D_oo)^dag (eta_o - Moe Mee^{-1} eta_e)
   *
   * Wilson:
   *      (D_oo)^{\dag} D_oo psi_o = (D_oo)^dag L^{-1}  eta_o
   * Stag:
   *      D_oo psi_o = L^{-1}  eta =    (eta_o - Moe Mee^{-1} eta_e)
   *
   * L^-1 eta_o= (1              0 ) (e
   *             (-MoeMee^{-1}   1 )   
   *
   *Even
   * ii)  Mee psi_e + Meo psi_o = src_e
   *
   *   => sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
   *
   * 
   * TODO: Other options:
   * 
   * a) change checkerboards for Schur e<->o
   *
   * Left precon by Moo^-1
   * b) Doo^{dag} M_oo^-dag Moo^-1 Doo psi_0 =  (D_oo)^dag M_oo^-dag Moo^-1 L^{-1}  eta_o
   *                              eta_o'     = (D_oo)^dag  M_oo^-dag Moo^-1 (eta_o - Moe Mee^{-1} eta_e)
   *
   * Right precon by Moo^-1
   * c) M_oo^-dag Doo^{dag} Doo Moo^-1 phi_0 = M_oo^-dag (D_oo)^dag L^{-1}  eta_o
   *                              eta_o'     = M_oo^-dag (D_oo)^dag (eta_o - Moe Mee^{-1} eta_e)
   *                              psi_o = M_oo^-1 phi_o
   * TODO: Deflation 
   */
namespace Grid {

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Take a matrix and form a Red Black solver calling a Herm solver
  // Use of RB info prevents making SchurRedBlackSolve conform to standard interface
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Now make the norm reflect extra factor of Mee
  template<class Field> class SchurRedBlackStaggeredSolve {
  private:
    OperatorFunction<Field> & _HermitianRBSolver;
    int CBfactorise;
    bool subGuess;
  public:

    /////////////////////////////////////////////////////
    // Wrap the usual normal equations Schur trick
    /////////////////////////////////////////////////////
  SchurRedBlackStaggeredSolve(OperatorFunction<Field> &HermitianRBSolver, const bool initSubGuess = false)  :
     _HermitianRBSolver(HermitianRBSolver) 
    { 
      CBfactorise=0;
      subtractGuess(initSubGuess);
    };
    void subtractGuess(const bool initSubGuess)
    {
      subGuess = initSubGuess;
    }
    bool isSubtractGuess(void)
    {
      return subGuess;
    }

    template<class Matrix>
    void operator() (Matrix & _Matrix,const Field &in, Field &out){
      ZeroGuesser<Field> guess;
      (*this)(_Matrix,in,out,guess);
    }
    template<class Matrix, class Guesser>
    void operator() (Matrix & _Matrix,const Field &in, Field &out, Guesser &guess){

      // FIXME CGdiagonalMee not implemented virtual function
      // FIXME use CBfactorise to control schur decomp
      GridBase *grid = _Matrix.RedBlackGrid();
      GridBase *fgrid= _Matrix.Grid();

      SchurStaggeredOperator<Matrix,Field> _HermOpEO(_Matrix);
 
      Field src_e(grid);
      Field src_o(grid);
      Field sol_e(grid);
      Field sol_o(grid);
      Field   tmp(grid);
      Field  Mtmp(grid);
      Field resid(fgrid);
      
      std::cout << GridLogMessage << " SchurRedBlackStaggeredSolve " <<std::endl;
      pickCheckerboard(Even,src_e,in);
      pickCheckerboard(Odd ,src_o,in);
      pickCheckerboard(Even,sol_e,out);
      pickCheckerboard(Odd ,sol_o,out);
      std::cout << GridLogMessage << " SchurRedBlackStaggeredSolve checkerboards picked" <<std::endl;
    
      /////////////////////////////////////////////////////
      // src_o = (source_o - Moe MeeInv source_e)
      /////////////////////////////////////////////////////
      _Matrix.MooeeInv(src_e,tmp);     assert(  tmp.checkerboard ==Even);
      _Matrix.Meooe   (tmp,Mtmp);      assert( Mtmp.checkerboard ==Odd);     
      tmp=src_o-Mtmp;                  assert(  tmp.checkerboard ==Odd);     

      //src_o = tmp;     assert(src_o.checkerboard ==Odd);
      _Matrix.Mooee(tmp,src_o); // Extra factor of "m" in source from dumb choice of matrix norm.

      //////////////////////////////////////////////////////////////
      // Call the red-black solver
      //////////////////////////////////////////////////////////////
      std::cout<<GridLogMessage << "SchurRedBlackStaggeredSolver calling the Mpc solver" <<std::endl;
      guess(src_o, sol_o);
      Mtmp = sol_o;
      _HermitianRBSolver(_HermOpEO,src_o,sol_o);  assert(sol_o.checkerboard==Odd);
      std::cout<<GridLogMessage << "SchurRedBlackStaggeredSolver called  the Mpc solver" <<std::endl;
      // Fionn A2A boolean behavioural control
      if (subGuess)        sol_o = sol_o-Mtmp;

      ///////////////////////////////////////////////////
      // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
      ///////////////////////////////////////////////////
      _Matrix.Meooe(sol_o,tmp);        assert(  tmp.checkerboard   ==Even);
      src_e = src_e-tmp;               assert(  src_e.checkerboard ==Even);
      _Matrix.MooeeInv(src_e,sol_e);   assert(  sol_e.checkerboard ==Even);
     
      std::cout<<GridLogMessage << "SchurRedBlackStaggeredSolver reconstructed other CB" <<std::endl;
      setCheckerboard(out,sol_e); assert(  sol_e.checkerboard ==Even);
      setCheckerboard(out,sol_o); assert(  sol_o.checkerboard ==Odd );
      std::cout<<GridLogMessage << "SchurRedBlackStaggeredSolver inserted solution" <<std::endl;

      // Verify the unprec residual
      if ( ! subGuess ) {
        _Matrix.M(out,resid); 
        resid = resid-in;
        RealD ns = norm2(in);
        RealD nr = norm2(resid);
        std::cout<<GridLogMessage << "SchurRedBlackStaggered solver true unprec resid "<< std::sqrt(nr/ns) <<" nr "<< nr <<" ns "<<ns << std::endl;
      } else {
        std::cout << GridLogMessage << "Guess subtracted after solve." << std::endl;
      }
    }     
  };
  template<class Field> using SchurRedBlackStagSolve = SchurRedBlackStaggeredSolve<Field>;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Take a matrix and form a Red Black solver calling a Herm solver
  // Use of RB info prevents making SchurRedBlackSolve conform to standard interface
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class Field> class SchurRedBlackDiagMooeeSolve {
  private:
    OperatorFunction<Field> & _HermitianRBSolver;
    int CBfactorise;
    bool subGuess;
  public:

    /////////////////////////////////////////////////////
    // Wrap the usual normal equations Schur trick
    /////////////////////////////////////////////////////
  SchurRedBlackDiagMooeeSolve(OperatorFunction<Field> &HermitianRBSolver,int cb=0, const bool initSubGuess = false)  :  _HermitianRBSolver(HermitianRBSolver) 
  { 
    CBfactorise=cb;
    subtractGuess(initSubGuess);
  };
    void subtractGuess(const bool initSubGuess)
    {
      subGuess = initSubGuess;
    }
    bool isSubtractGuess(void)
    {
      return subGuess;
    }
    template<class Matrix>
    void operator() (Matrix & _Matrix,const Field &in, Field &out){
      ZeroGuesser<Field> guess;
      (*this)(_Matrix,in,out,guess);
    }
    template<class Matrix, class Guesser>
    void operator() (Matrix & _Matrix,const Field &in, Field &out,Guesser &guess){

      // FIXME CGdiagonalMee not implemented virtual function
      // FIXME use CBfactorise to control schur decomp
      GridBase *grid = _Matrix.RedBlackGrid();
      GridBase *fgrid= _Matrix.Grid();

      SchurDiagMooeeOperator<Matrix,Field> _HermOpEO(_Matrix);
 
      Field src_e(grid);
      Field src_o(grid);
      Field sol_e(grid);
      Field sol_o(grid);
      Field   tmp(grid);
      Field  Mtmp(grid);
      Field resid(fgrid);

      pickCheckerboard(Even,src_e,in);
      pickCheckerboard(Odd ,src_o,in);
      pickCheckerboard(Even,sol_e,out);
      pickCheckerboard(Odd ,sol_o,out);
    
      /////////////////////////////////////////////////////
      // src_o = Mdag * (source_o - Moe MeeInv source_e)
      /////////////////////////////////////////////////////
      _Matrix.MooeeInv(src_e,tmp);     assert(  tmp.checkerboard ==Even);
      _Matrix.Meooe   (tmp,Mtmp);      assert( Mtmp.checkerboard ==Odd);     
      tmp=src_o-Mtmp;                  assert(  tmp.checkerboard ==Odd);     

      // get the right MpcDag
      _HermOpEO.MpcDag(tmp,src_o);     assert(src_o.checkerboard ==Odd);       

      //////////////////////////////////////////////////////////////
      // Call the red-black solver
      //////////////////////////////////////////////////////////////
      std::cout<<GridLogMessage << "SchurRedBlack solver calling the MpcDagMp solver" <<std::endl;
      guess(src_o,sol_o);
      Mtmp = sol_o;
      _HermitianRBSolver(_HermOpEO,src_o,sol_o);  assert(sol_o.checkerboard==Odd);
      // Fionn A2A boolean behavioural control
      if (subGuess)        sol_o = sol_o-Mtmp;

      ///////////////////////////////////////////////////
      // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
      ///////////////////////////////////////////////////
      _Matrix.Meooe(sol_o,tmp);        assert(  tmp.checkerboard   ==Even);
      src_e = src_e-tmp;               assert(  src_e.checkerboard ==Even);
      _Matrix.MooeeInv(src_e,sol_e);   assert(  sol_e.checkerboard ==Even);
     
      setCheckerboard(out,sol_e); assert(  sol_e.checkerboard ==Even);
      setCheckerboard(out,sol_o); assert(  sol_o.checkerboard ==Odd );

      // Verify the unprec residual
      if ( ! subGuess ) {
        _Matrix.M(out,resid); 
        resid = resid-in;
        RealD ns = norm2(in);
        RealD nr = norm2(resid);

        std::cout<<GridLogMessage << "SchurRedBlackDiagMooee solver true unprec resid "<< std::sqrt(nr/ns) <<" nr "<< nr <<" ns "<<ns << std::endl;
      } else {
        std::cout << GridLogMessage << "Guess subtracted after solve." << std::endl;
      }
    }     
  };


  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Take a matrix and form a Red Black solver calling a Herm solver
  // Use of RB info prevents making SchurRedBlackSolve conform to standard interface
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class Field> class SchurRedBlackDiagTwoSolve {
  private:
    OperatorFunction<Field> & _HermitianRBSolver;
    int CBfactorise;
    bool subGuess;
  public:

    /////////////////////////////////////////////////////
    // Wrap the usual normal equations Schur trick
    /////////////////////////////////////////////////////
  SchurRedBlackDiagTwoSolve(OperatorFunction<Field> &HermitianRBSolver, const bool initSubGuess = false)  :
     _HermitianRBSolver(HermitianRBSolver) 
    { 
      CBfactorise = 0;
      subtractGuess(initSubGuess);
    };
    void subtractGuess(const bool initSubGuess)
    {
      subGuess = initSubGuess;
    }
    bool isSubtractGuess(void)
    {
      return subGuess;
    }

    template<class Matrix>
    void operator() (Matrix & _Matrix,const Field &in, Field &out){
      ZeroGuesser<Field> guess;
      (*this)(_Matrix,in,out,guess);
    }
    template<class Matrix,class Guesser>
    void operator() (Matrix & _Matrix,const Field &in, Field &out,Guesser &guess){

      // FIXME CGdiagonalMee not implemented virtual function
      // FIXME use CBfactorise to control schur decomp
      GridBase *grid = _Matrix.RedBlackGrid();
      GridBase *fgrid= _Matrix.Grid();

      SchurDiagTwoOperator<Matrix,Field> _HermOpEO(_Matrix);
 
      Field src_e(grid);
      Field src_o(grid);
      Field sol_e(grid);
      Field sol_o(grid);
      Field   tmp(grid);
      Field  Mtmp(grid);
      Field resid(fgrid);

      pickCheckerboard(Even,src_e,in);
      pickCheckerboard(Odd ,src_o,in);
      pickCheckerboard(Even,sol_e,out);
      pickCheckerboard(Odd ,sol_o,out);
    
      /////////////////////////////////////////////////////
      // src_o = Mdag * (source_o - Moe MeeInv source_e)
      /////////////////////////////////////////////////////
      _Matrix.MooeeInv(src_e,tmp);     assert(  tmp.checkerboard ==Even);
      _Matrix.Meooe   (tmp,Mtmp);      assert( Mtmp.checkerboard ==Odd);     
      tmp=src_o-Mtmp;                  assert(  tmp.checkerboard ==Odd);     

      // get the right MpcDag
      _HermOpEO.MpcDag(tmp,src_o);     assert(src_o.checkerboard ==Odd);       

      //////////////////////////////////////////////////////////////
      // Call the red-black solver
      //////////////////////////////////////////////////////////////
      std::cout<<GridLogMessage << "SchurRedBlack solver calling the MpcDagMp solver" <<std::endl;
//      _HermitianRBSolver(_HermOpEO,src_o,sol_o);  assert(sol_o.checkerboard==Odd);
      guess(src_o,tmp);
      Mtmp = tmp;
      _HermitianRBSolver(_HermOpEO,src_o,tmp);  assert(tmp.checkerboard==Odd);
      // Fionn A2A boolean behavioural control
      if (subGuess)      tmp = tmp-Mtmp;
      _Matrix.MooeeInv(tmp,sol_o);       assert(  sol_o.checkerboard   ==Odd);

      ///////////////////////////////////////////////////
      // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
      ///////////////////////////////////////////////////
      _Matrix.Meooe(sol_o,tmp);        assert(  tmp.checkerboard   ==Even);
      src_e = src_e-tmp;               assert(  src_e.checkerboard ==Even);
      _Matrix.MooeeInv(src_e,sol_e);   assert(  sol_e.checkerboard ==Even);
     
      setCheckerboard(out,sol_e); assert(  sol_e.checkerboard ==Even);
      setCheckerboard(out,sol_o); assert(  sol_o.checkerboard ==Odd );

      // Verify the unprec residual
      if ( ! subGuess ) {
        _Matrix.M(out,resid); 
        resid = resid-in;
        RealD ns = norm2(in);
        RealD nr = norm2(resid);

        std::cout<<GridLogMessage << "SchurRedBlackDiagTwo solver true unprec resid "<< std::sqrt(nr/ns) <<" nr "<< nr <<" ns "<<ns << std::endl;
      } else {
        std::cout << GridLogMessage << "Guess subtracted after solve." << std::endl;
      }
    }     
  };
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Take a matrix and form a Red Black solver calling a Herm solver
  // Use of RB info prevents making SchurRedBlackSolve conform to standard interface
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class Field> class SchurRedBlackDiagTwoMixed {
  private:
    LinearFunction<Field> & _HermitianRBSolver;
    int CBfactorise;
    bool subGuess;
  public:

    /////////////////////////////////////////////////////
    // Wrap the usual normal equations Schur trick
    /////////////////////////////////////////////////////
  SchurRedBlackDiagTwoMixed(LinearFunction<Field> &HermitianRBSolver, const bool initSubGuess = false)  :
     _HermitianRBSolver(HermitianRBSolver) 
    { 
      CBfactorise=0;
      subtractGuess(initSubGuess);
    };
    void subtractGuess(const bool initSubGuess)
    {
      subGuess = initSubGuess;
    }
    bool isSubtractGuess(void)
    {
      return subGuess;
    }

    template<class Matrix>
    void operator() (Matrix & _Matrix,const Field &in, Field &out){
      ZeroGuesser<Field> guess;
      (*this)(_Matrix,in,out,guess);
    }
    template<class Matrix, class Guesser>
    void operator() (Matrix & _Matrix,const Field &in, Field &out,Guesser &guess){

      // FIXME CGdiagonalMee not implemented virtual function
      // FIXME use CBfactorise to control schur decomp
      GridBase *grid = _Matrix.RedBlackGrid();
      GridBase *fgrid= _Matrix.Grid();

      SchurDiagTwoOperator<Matrix,Field> _HermOpEO(_Matrix);
 
      Field src_e(grid);
      Field src_o(grid);
      Field sol_e(grid);
      Field sol_o(grid);
      Field   tmp(grid);
      Field  Mtmp(grid);
      Field resid(fgrid);

      pickCheckerboard(Even,src_e,in);
      pickCheckerboard(Odd ,src_o,in);
      pickCheckerboard(Even,sol_e,out);
      pickCheckerboard(Odd ,sol_o,out);
    
      /////////////////////////////////////////////////////
      // src_o = Mdag * (source_o - Moe MeeInv source_e)
      /////////////////////////////////////////////////////
      _Matrix.MooeeInv(src_e,tmp);     assert(  tmp.checkerboard ==Even);
      _Matrix.Meooe   (tmp,Mtmp);      assert( Mtmp.checkerboard ==Odd);     
      tmp=src_o-Mtmp;                  assert(  tmp.checkerboard ==Odd);     

      // get the right MpcDag
      _HermOpEO.MpcDag(tmp,src_o);     assert(src_o.checkerboard ==Odd);       

      //////////////////////////////////////////////////////////////
      // Call the red-black solver
      //////////////////////////////////////////////////////////////
      std::cout<<GridLogMessage << "SchurRedBlack solver calling the MpcDagMp solver" <<std::endl;
//      _HermitianRBSolver(_HermOpEO,src_o,sol_o);  assert(sol_o.checkerboard==Odd);
//      _HermitianRBSolver(_HermOpEO,src_o,tmp);  assert(tmp.checkerboard==Odd);
      guess(src_o,tmp);
      Mtmp = tmp;
      _HermitianRBSolver(_HermOpEO,src_o,tmp);  assert(tmp.checkerboard==Odd);
      // Fionn A2A boolean behavioural control
      if (subGuess)      tmp = tmp-Mtmp;
      _Matrix.MooeeInv(tmp,sol_o);        assert(  sol_o.checkerboard   ==Odd);

      ///////////////////////////////////////////////////
      // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
      ///////////////////////////////////////////////////
      _Matrix.Meooe(sol_o,tmp);        assert(  tmp.checkerboard   ==Even);
      src_e = src_e-tmp;               assert(  src_e.checkerboard ==Even);
      _Matrix.MooeeInv(src_e,sol_e);   assert(  sol_e.checkerboard ==Even);
     
      setCheckerboard(out,sol_e); assert(  sol_e.checkerboard ==Even);
      setCheckerboard(out,sol_o); assert(  sol_o.checkerboard ==Odd );

      // Verify the unprec residual
      if ( ! subGuess ) {
        _Matrix.M(out,resid); 
        resid = resid-in;
        RealD ns = norm2(in);
        RealD nr = norm2(resid);

        std::cout << GridLogMessage << "SchurRedBlackDiagTwo solver true unprec resid " << std::sqrt(nr / ns) << " nr " << nr << " ns " << ns << std::endl;
      } else {
        std::cout << GridLogMessage << "Guess subtracted after solve." << std::endl;
      }
    }     
  };

}
#endif
