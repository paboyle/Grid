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
  // Use base class to share code
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Take a matrix and form a Red Black solver calling a Herm solver
  // Use of RB info prevents making SchurRedBlackSolve conform to standard interface
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class Field> class SchurRedBlackBase {
  protected:
    typedef CheckerBoardedSparseMatrixBase<Field> Matrix;
    OperatorFunction<Field> & _HermitianRBSolver;
    int CBfactorise;
    bool subGuess;
    bool useSolnAsInitGuess; // if true user-supplied solution vector is used as initial guess for solver
  public:

    SchurRedBlackBase(OperatorFunction<Field> &HermitianRBSolver, const bool initSubGuess = false,
        const bool _solnAsInitGuess = false)  :
    _HermitianRBSolver(HermitianRBSolver),
    useSolnAsInitGuess(_solnAsInitGuess)
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

    /////////////////////////////////////////////////////////////
    // Shared code
    /////////////////////////////////////////////////////////////
    void operator() (Matrix & _Matrix,const Field &in, Field &out){
      ZeroGuesser<Field> guess;
      (*this)(_Matrix,in,out,guess);
    }
    void operator()(Matrix &_Matrix, const std::vector<Field> &in, std::vector<Field> &out) 
    {
      ZeroGuesser<Field> guess;
      (*this)(_Matrix,in,out,guess);
    }

    void RedBlackSource(Matrix &_Matrix, const std::vector<Field> &in, std::vector<Field> &src_o) 
    {
      GridBase *grid = _Matrix.RedBlackGrid();
      Field tmp(grid);
      int nblock = in.size();
      for(int b=0;b<nblock;b++){
	RedBlackSource(_Matrix,in[b],tmp,src_o[b]);
      }
    }
    // James can write his own deflated guesser
    // with optimised code for the inner products
    //    RedBlackSolveSplitGrid();
    //    RedBlackSolve(_Matrix,src_o,sol_o); 

    void RedBlackSolution(Matrix &_Matrix, const std::vector<Field> &in, const std::vector<Field> &sol_o, std::vector<Field> &out)
    {
      GridBase *grid = _Matrix.RedBlackGrid();
      Field tmp(grid);
      int nblock = in.size();
      for(int b=0;b<nblock;b++) {
	pickCheckerboard(Even,tmp,in[b]);
	RedBlackSolution(_Matrix,sol_o[b],tmp,out[b]);
      }
    }

    template<class Guesser>
    void operator()(Matrix &_Matrix, const std::vector<Field> &in, std::vector<Field> &out,Guesser &guess) 
    {
      GridBase *grid = _Matrix.RedBlackGrid();
      GridBase *fgrid= _Matrix.Grid();
      int nblock = in.size();

      std::vector<Field> src_o(nblock,grid);
      std::vector<Field> sol_o(nblock,grid);
      
      std::vector<Field> guess_save;

      Field resid(fgrid);
      Field tmp(grid);

      ////////////////////////////////////////////////
      // Prepare RedBlack source
      ////////////////////////////////////////////////
      RedBlackSource(_Matrix,in,src_o);
	//      for(int b=0;b<nblock;b++){
	//	RedBlackSource(_Matrix,in[b],tmp,src_o[b]);
	//      }
      
      ////////////////////////////////////////////////
      // Make the guesses
      ////////////////////////////////////////////////
      if ( subGuess ) guess_save.resize(nblock,grid);

      
      if(useSolnAsInitGuess) {
        for(int b=0;b<nblock;b++){
          pickCheckerboard(Odd, sol_o[b], out[b]);
        }
      } else {
        guess(src_o, sol_o); 
      }

	    if ( subGuess ) { 
        for(int b=0;b<nblock;b++){
          guess_save[b] = sol_o[b];
        }
      }
      //////////////////////////////////////////////////////////////
      // Call the block solver
      //////////////////////////////////////////////////////////////
      std::cout<<GridLogMessage << "SchurRedBlackBase calling the solver for "<<nblock<<" RHS" <<std::endl;
      RedBlackSolve(_Matrix,src_o,sol_o);

      ////////////////////////////////////////////////
      // A2A boolean behavioural control & reconstruct other checkerboard
      ////////////////////////////////////////////////
      for(int b=0;b<nblock;b++) {

	if (subGuess)   sol_o[b] = sol_o[b] - guess_save[b];

	///////// Needs even source //////////////
	pickCheckerboard(Even,tmp,in[b]);
	RedBlackSolution(_Matrix,sol_o[b],tmp,out[b]);

	/////////////////////////////////////////////////
	// Check unprec residual if possible
	/////////////////////////////////////////////////
	if ( ! subGuess ) {
	  _Matrix.M(out[b],resid); 
	  resid = resid-in[b];
	  RealD ns = norm2(in[b]);
	  RealD nr = norm2(resid);
	
	  std::cout<<GridLogMessage<< "SchurRedBlackBase solver true unprec resid["<<b<<"] "<<std::sqrt(nr/ns) << std::endl;
	} else {
	  std::cout<<GridLogMessage<< "SchurRedBlackBase Guess subtracted after solve["<<b<<"] " << std::endl;
	}

      }
    }
    template<class Guesser>
    void operator() (Matrix & _Matrix,const Field &in, Field &out,Guesser &guess){

      // FIXME CGdiagonalMee not implemented virtual function
      // FIXME use CBfactorise to control schur decomp
      GridBase *grid = _Matrix.RedBlackGrid();
      GridBase *fgrid= _Matrix.Grid();

      Field resid(fgrid);
      Field src_o(grid);
      Field src_e(grid);
      Field sol_o(grid);

      ////////////////////////////////////////////////
      // RedBlack source
      ////////////////////////////////////////////////
      RedBlackSource(_Matrix,in,src_e,src_o);

      ////////////////////////////////
      // Construct the guess
      ////////////////////////////////
      if(useSolnAsInitGuess) {
        pickCheckerboard(Odd, sol_o, out);
      } else {
        guess(src_o,sol_o);
      }

      Field  guess_save(grid);
      guess_save = sol_o;

      //////////////////////////////////////////////////////////////
      // Call the red-black solver
      //////////////////////////////////////////////////////////////
      RedBlackSolve(_Matrix,src_o,sol_o);

      ////////////////////////////////////////////////
      // Fionn A2A boolean behavioural control
      ////////////////////////////////////////////////
      if (subGuess)      sol_o= sol_o-guess_save;

      ///////////////////////////////////////////////////
      // RedBlack solution needs the even source
      ///////////////////////////////////////////////////
      RedBlackSolution(_Matrix,sol_o,src_e,out);

      // Verify the unprec residual
      if ( ! subGuess ) {
        _Matrix.M(out,resid); 
        resid = resid-in;
        RealD ns = norm2(in);
        RealD nr = norm2(resid);

        std::cout<<GridLogMessage << "SchurRedBlackBase solver true unprec resid "<< std::sqrt(nr/ns) << std::endl;
      } else {
        std::cout << GridLogMessage << "SchurRedBlackBase Guess subtracted after solve." << std::endl;
      }
    }     
    
    /////////////////////////////////////////////////////////////
    // Override in derived. 
    /////////////////////////////////////////////////////////////
    virtual void RedBlackSource  (Matrix & _Matrix,const Field &src, Field &src_e,Field &src_o)                =0;
    virtual void RedBlackSolution(Matrix & _Matrix,const Field &sol_o, const Field &src_e,Field &sol)          =0;
    virtual void RedBlackSolve   (Matrix & _Matrix,const Field &src_o, Field &sol_o)                           =0;
    virtual void RedBlackSolve   (Matrix & _Matrix,const std::vector<Field> &src_o,  std::vector<Field> &sol_o)=0;

  };

  template<class Field> class SchurRedBlackStaggeredSolve : public SchurRedBlackBase<Field> {
  public:
    typedef CheckerBoardedSparseMatrixBase<Field> Matrix;

    SchurRedBlackStaggeredSolve(OperatorFunction<Field> &HermitianRBSolver, const bool initSubGuess = false,
        const bool _solnAsInitGuess = false) 
      :    SchurRedBlackBase<Field> (HermitianRBSolver,initSubGuess,_solnAsInitGuess) 
    {
    }

    //////////////////////////////////////////////////////
    // Override RedBlack specialisation
    //////////////////////////////////////////////////////
    virtual void RedBlackSource(Matrix & _Matrix,const Field &src, Field &src_e,Field &src_o)
    {
      GridBase *grid = _Matrix.RedBlackGrid();
      GridBase *fgrid= _Matrix.Grid();

      Field   tmp(grid);
      Field  Mtmp(grid);

      pickCheckerboard(Even,src_e,src);
      pickCheckerboard(Odd ,src_o,src);

      /////////////////////////////////////////////////////
      // src_o = (source_o - Moe MeeInv source_e)
      /////////////////////////////////////////////////////
      _Matrix.MooeeInv(src_e,tmp);     assert(  tmp.Checkerboard() ==Even);
      _Matrix.Meooe   (tmp,Mtmp);      assert( Mtmp.Checkerboard() ==Odd);     
      tmp=src_o-Mtmp;                  assert(  tmp.Checkerboard() ==Odd);     

      _Matrix.Mooee(tmp,src_o); // Extra factor of "m" in source from dumb choice of matrix norm.
    }
    virtual void RedBlackSolution(Matrix & _Matrix,const Field &sol_o, const Field &src_e_c,Field &sol)
    {
      GridBase *grid = _Matrix.RedBlackGrid();
      GridBase *fgrid= _Matrix.Grid();

      Field   tmp(grid);
      Field   sol_e(grid);
      Field   src_e(grid);

      src_e = src_e_c; // Const correctness

      ///////////////////////////////////////////////////
      // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
      ///////////////////////////////////////////////////
      _Matrix.Meooe(sol_o,tmp);        assert(  tmp.Checkerboard()   ==Even);
      src_e = src_e-tmp;               assert(  src_e.Checkerboard() ==Even);
      _Matrix.MooeeInv(src_e,sol_e);   assert(  sol_e.Checkerboard() ==Even);
     
      setCheckerboard(sol,sol_e); assert(  sol_e.Checkerboard() ==Even);
      setCheckerboard(sol,sol_o); assert(  sol_o.Checkerboard() ==Odd );
    }
    virtual void RedBlackSolve   (Matrix & _Matrix,const Field &src_o, Field &sol_o)
    {
      SchurStaggeredOperator<Matrix,Field> _HermOpEO(_Matrix);
      this->_HermitianRBSolver(_HermOpEO,src_o,sol_o);  assert(sol_o.Checkerboard()==Odd);
    };
    virtual void RedBlackSolve   (Matrix & _Matrix,const std::vector<Field> &src_o,  std::vector<Field> &sol_o)
    {
      SchurStaggeredOperator<Matrix,Field> _HermOpEO(_Matrix);
      this->_HermitianRBSolver(_HermOpEO,src_o,sol_o); 
    }
  };
  template<class Field> using SchurRedBlackStagSolve = SchurRedBlackStaggeredSolve<Field>;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Site diagonal has Mooee on it.
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class Field> class SchurRedBlackDiagMooeeSolve : public SchurRedBlackBase<Field> {
  public:
    typedef CheckerBoardedSparseMatrixBase<Field> Matrix;

    SchurRedBlackDiagMooeeSolve(OperatorFunction<Field> &HermitianRBSolver, const bool initSubGuess = false,
        const bool _solnAsInitGuess = false)  
      : SchurRedBlackBase<Field> (HermitianRBSolver,initSubGuess,_solnAsInitGuess) {};


    //////////////////////////////////////////////////////
    // Override RedBlack specialisation
    //////////////////////////////////////////////////////
    virtual void RedBlackSource(Matrix & _Matrix,const Field &src, Field &src_e,Field &src_o)
    {
      GridBase *grid = _Matrix.RedBlackGrid();
      GridBase *fgrid= _Matrix.Grid();

      Field   tmp(grid);
      Field  Mtmp(grid);

      pickCheckerboard(Even,src_e,src);
      pickCheckerboard(Odd ,src_o,src);

      /////////////////////////////////////////////////////
      // src_o = Mdag * (source_o - Moe MeeInv source_e)
      /////////////////////////////////////////////////////
      _Matrix.MooeeInv(src_e,tmp);     assert(  tmp.Checkerboard() ==Even);
      _Matrix.Meooe   (tmp,Mtmp);      assert( Mtmp.Checkerboard() ==Odd);     
      tmp=src_o-Mtmp;                  assert(  tmp.Checkerboard() ==Odd);     

      // get the right MpcDag
      SchurDiagMooeeOperator<Matrix,Field> _HermOpEO(_Matrix);
      _HermOpEO.MpcDag(tmp,src_o);     assert(src_o.Checkerboard() ==Odd);       

    }
    virtual void RedBlackSolution(Matrix & _Matrix,const Field &sol_o, const Field &src_e,Field &sol)
    {
      GridBase *grid = _Matrix.RedBlackGrid();
      GridBase *fgrid= _Matrix.Grid();

      Field   tmp(grid);
      Field  sol_e(grid);
      Field  src_e_i(grid);
      ///////////////////////////////////////////////////
      // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
      ///////////////////////////////////////////////////
      _Matrix.Meooe(sol_o,tmp);          assert(  tmp.Checkerboard()   ==Even);
      src_e_i = src_e-tmp;               assert(  src_e_i.Checkerboard() ==Even);
      _Matrix.MooeeInv(src_e_i,sol_e);   assert(  sol_e.Checkerboard() ==Even);
     
      setCheckerboard(sol,sol_e); assert(  sol_e.Checkerboard() ==Even);
      setCheckerboard(sol,sol_o); assert(  sol_o.Checkerboard() ==Odd );
    }
    virtual void RedBlackSolve   (Matrix & _Matrix,const Field &src_o, Field &sol_o)
    {
      SchurDiagMooeeOperator<Matrix,Field> _HermOpEO(_Matrix);
      this->_HermitianRBSolver(_HermOpEO,src_o,sol_o);  assert(sol_o.Checkerboard()==Odd);
    };
    virtual void RedBlackSolve   (Matrix & _Matrix,const std::vector<Field> &src_o,  std::vector<Field> &sol_o)
    {
      SchurDiagMooeeOperator<Matrix,Field> _HermOpEO(_Matrix);
      this->_HermitianRBSolver(_HermOpEO,src_o,sol_o); 
    }
  };

  template<class Field> class NonHermitianSchurRedBlackDiagMooeeSolve : public SchurRedBlackBase<Field> 
  {
    public:
      typedef CheckerBoardedSparseMatrixBase<Field> Matrix;

      NonHermitianSchurRedBlackDiagMooeeSolve(OperatorFunction<Field>& RBSolver, const bool initSubGuess = false,
          const bool _solnAsInitGuess = false)  
      : SchurRedBlackBase<Field>(RBSolver, initSubGuess, _solnAsInitGuess) {};

      //////////////////////////////////////////////////////
      // Override RedBlack specialisation
      //////////////////////////////////////////////////////
      virtual void RedBlackSource(Matrix& _Matrix, const Field& src, Field& src_e, Field& src_o)
      {
        GridBase* grid  = _Matrix.RedBlackGrid();
        GridBase* fgrid = _Matrix.Grid();

        Field  tmp(grid);
        Field Mtmp(grid);

        pickCheckerboard(Even, src_e, src);
        pickCheckerboard(Odd , src_o, src);

        /////////////////////////////////////////////////////
        // src_o = Mdag * (source_o - Moe MeeInv source_e)
        /////////////////////////////////////////////////////
        _Matrix.MooeeInv(src_e, tmp);   assert(   tmp.Checkerboard() == Even );
        _Matrix.Meooe   (tmp, Mtmp);    assert(  Mtmp.Checkerboard() == Odd  );     
        src_o -= Mtmp;                  assert( src_o.Checkerboard() == Odd  );     
      }
      
      virtual void RedBlackSolution(Matrix& _Matrix, const Field& sol_o, const Field& src_e, Field& sol)
      {
        GridBase* grid  = _Matrix.RedBlackGrid();
        GridBase* fgrid = _Matrix.Grid();

        Field     tmp(grid);
        Field   sol_e(grid);
        Field src_e_i(grid);
        
        ///////////////////////////////////////////////////
        // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
        ///////////////////////////////////////////////////
        _Matrix.Meooe(sol_o, tmp);         assert(     tmp.Checkerboard() == Even );
        src_e_i = src_e - tmp;             assert( src_e_i.Checkerboard() == Even );
        _Matrix.MooeeInv(src_e_i, sol_e);  assert(   sol_e.Checkerboard() == Even );
       
        setCheckerboard(sol, sol_e); assert( sol_e.Checkerboard() == Even );
        setCheckerboard(sol, sol_o); assert( sol_o.Checkerboard() == Odd  );
      }

      virtual void RedBlackSolve(Matrix& _Matrix, const Field& src_o, Field& sol_o)
      {
        NonHermitianSchurDiagMooeeOperator<Matrix,Field> _OpEO(_Matrix);
        this->_HermitianRBSolver(_OpEO, src_o, sol_o);  assert(sol_o.Checkerboard() == Odd);
      }

      virtual void RedBlackSolve(Matrix& _Matrix, const std::vector<Field>& src_o, std::vector<Field>& sol_o)
      {
        NonHermitianSchurDiagMooeeOperator<Matrix,Field> _OpEO(_Matrix);
        this->_HermitianRBSolver(_OpEO, src_o, sol_o); 
      }
  };

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Site diagonal is identity, left preconditioned by Mee^inv
  // ( 1 - Mee^inv Meo Moo^inv Moe ) phi = Mee_inv ( Mee - Meo Moo^inv Moe Mee^inv  ) phi =  Mee_inv eta
  //
  // Solve:
  // ( 1 - Mee^inv Meo Moo^inv Moe )^dag ( 1 - Mee^inv Meo Moo^inv Moe ) phi = ( 1 - Mee^inv Meo Moo^inv Moe )^dag  Mee_inv eta
  //
  // Old notation e<->o
  //
  // Left precon by Moo^-1
  //  b) (Doo^{dag} M_oo^-dag) (Moo^-1 Doo) psi_o =  [ (D_oo)^dag M_oo^-dag ] Moo^-1 L^{-1}  eta_o
  //                                   eta_o'     = (D_oo)^dag  M_oo^-dag Moo^-1 (eta_o - Moe Mee^{-1} eta_e)
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class Field> class SchurRedBlackDiagOneSolve : public SchurRedBlackBase<Field> {
  public:
    typedef CheckerBoardedSparseMatrixBase<Field> Matrix;

    /////////////////////////////////////////////////////
    // Wrap the usual normal equations Schur trick
    /////////////////////////////////////////////////////
  SchurRedBlackDiagOneSolve(OperatorFunction<Field> &HermitianRBSolver, const bool initSubGuess = false,
      const bool _solnAsInitGuess = false)  
    : SchurRedBlackBase<Field>(HermitianRBSolver,initSubGuess,_solnAsInitGuess) {};

    virtual void RedBlackSource(Matrix & _Matrix,const Field &src, Field &src_e,Field &src_o)
    {
      GridBase *grid = _Matrix.RedBlackGrid();
      GridBase *fgrid= _Matrix.Grid();

      SchurDiagOneOperator<Matrix,Field> _HermOpEO(_Matrix);
      
      Field   tmp(grid);
      Field  Mtmp(grid);

      pickCheckerboard(Even,src_e,src);
      pickCheckerboard(Odd ,src_o,src);
    
      /////////////////////////////////////////////////////
      // src_o = Mpcdag *MooeeInv * (source_o - Moe MeeInv source_e)
      /////////////////////////////////////////////////////
      _Matrix.MooeeInv(src_e,tmp);     assert(  tmp.Checkerboard() ==Even);
      _Matrix.Meooe   (tmp,Mtmp);      assert( Mtmp.Checkerboard() ==Odd);     
      Mtmp=src_o-Mtmp;                 
      _Matrix.MooeeInv(Mtmp,tmp);      assert( tmp.Checkerboard() ==Odd);     
      
      // get the right MpcDag
      _HermOpEO.MpcDag(tmp,src_o);     assert(src_o.Checkerboard() ==Odd);       
    }

    virtual void RedBlackSolution(Matrix & _Matrix,const Field &sol_o, const Field &src_e,Field &sol)
    {
      GridBase *grid = _Matrix.RedBlackGrid();
      GridBase *fgrid= _Matrix.Grid();

      Field   tmp(grid);
      Field   sol_e(grid);


      ///////////////////////////////////////////////////
      // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
      ///////////////////////////////////////////////////
      _Matrix.Meooe(sol_o,tmp);    assert(  tmp.Checkerboard()   ==Even);
      tmp = src_e-tmp;             assert(  src_e.Checkerboard() ==Even);
      _Matrix.MooeeInv(tmp,sol_e); assert(  sol_e.Checkerboard() ==Even);
     
      setCheckerboard(sol,sol_e);  assert(  sol_e.Checkerboard() ==Even);
      setCheckerboard(sol,sol_o);  assert(  sol_o.Checkerboard() ==Odd );
    };

    virtual void RedBlackSolve   (Matrix & _Matrix,const Field &src_o, Field &sol_o)
    {
      SchurDiagOneOperator<Matrix,Field> _HermOpEO(_Matrix);
      this->_HermitianRBSolver(_HermOpEO,src_o,sol_o);
    };
    virtual void RedBlackSolve   (Matrix & _Matrix,const std::vector<Field> &src_o,  std::vector<Field> &sol_o)
    {
      SchurDiagOneOperator<Matrix,Field> _HermOpEO(_Matrix);
      this->_HermitianRBSolver(_HermOpEO,src_o,sol_o); 
    }
  };

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Site diagonal is identity, right preconditioned by Mee^inv
  // ( 1 - Meo Moo^inv Moe Mee^inv  ) phi =( 1 - Meo Moo^inv Moe Mee^inv  ) Mee psi =  = eta  = eta
  //=> psi = MeeInv phi
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class Field> class SchurRedBlackDiagTwoSolve : public SchurRedBlackBase<Field> {
  public:
    typedef CheckerBoardedSparseMatrixBase<Field> Matrix;

    /////////////////////////////////////////////////////
    // Wrap the usual normal equations Schur trick
    /////////////////////////////////////////////////////
  SchurRedBlackDiagTwoSolve(OperatorFunction<Field> &HermitianRBSolver, const bool initSubGuess = false,
      const bool _solnAsInitGuess = false)  
    : SchurRedBlackBase<Field>(HermitianRBSolver,initSubGuess,_solnAsInitGuess) {};

    virtual void RedBlackSource(Matrix & _Matrix,const Field &src, Field &src_e,Field &src_o)
    {
      GridBase *grid = _Matrix.RedBlackGrid();
      GridBase *fgrid= _Matrix.Grid();

      SchurDiagTwoOperator<Matrix,Field> _HermOpEO(_Matrix);
      
      Field   tmp(grid);
      Field  Mtmp(grid);

      pickCheckerboard(Even,src_e,src);
      pickCheckerboard(Odd ,src_o,src);
    
      /////////////////////////////////////////////////////
      // src_o = Mdag * (source_o - Moe MeeInv source_e)
      /////////////////////////////////////////////////////
      _Matrix.MooeeInv(src_e,tmp);     assert(  tmp.Checkerboard() ==Even);
      _Matrix.Meooe   (tmp,Mtmp);      assert( Mtmp.Checkerboard() ==Odd);     
      tmp=src_o-Mtmp;                  assert(  tmp.Checkerboard() ==Odd);     

      // get the right MpcDag
      _HermOpEO.MpcDag(tmp,src_o);     assert(src_o.Checkerboard() ==Odd);       
    }

    virtual void RedBlackSolution(Matrix & _Matrix,const Field &sol_o, const Field &src_e,Field &sol)
    {
      GridBase *grid = _Matrix.RedBlackGrid();
      GridBase *fgrid= _Matrix.Grid();

      Field   sol_o_i(grid);
      Field   tmp(grid);
      Field   sol_e(grid);

      ////////////////////////////////////////////////
      // MooeeInv due to pecond
      ////////////////////////////////////////////////
      _Matrix.MooeeInv(sol_o,tmp);
      sol_o_i = tmp;

      ///////////////////////////////////////////////////
      // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
      ///////////////////////////////////////////////////
      _Matrix.Meooe(sol_o_i,tmp);    assert(  tmp.Checkerboard()   ==Even);
      tmp = src_e-tmp;               assert(  src_e.Checkerboard() ==Even);
      _Matrix.MooeeInv(tmp,sol_e);   assert(  sol_e.Checkerboard() ==Even);
     
      setCheckerboard(sol,sol_e);    assert(  sol_e.Checkerboard() ==Even);
      setCheckerboard(sol,sol_o_i);  assert(  sol_o_i.Checkerboard() ==Odd );
    };

    virtual void RedBlackSolve   (Matrix & _Matrix,const Field &src_o, Field &sol_o)
    {
      SchurDiagTwoOperator<Matrix,Field> _HermOpEO(_Matrix);
      this->_HermitianRBSolver(_HermOpEO,src_o,sol_o);
    };
    virtual void RedBlackSolve   (Matrix & _Matrix,const std::vector<Field> &src_o,  std::vector<Field> &sol_o)
    {
      SchurDiagTwoOperator<Matrix,Field> _HermOpEO(_Matrix);
      this->_HermitianRBSolver(_HermOpEO,src_o,sol_o); 
    }
  };

  template<class Field> class NonHermitianSchurRedBlackDiagTwoSolve : public SchurRedBlackBase<Field> 
  {
    public:
      typedef CheckerBoardedSparseMatrixBase<Field> Matrix;

      /////////////////////////////////////////////////////
      // Wrap the usual normal equations Schur trick
      /////////////////////////////////////////////////////
      NonHermitianSchurRedBlackDiagTwoSolve(OperatorFunction<Field>& RBSolver, const bool initSubGuess = false,
          const bool _solnAsInitGuess = false)  
      : SchurRedBlackBase<Field>(RBSolver, initSubGuess, _solnAsInitGuess) {};

      virtual void RedBlackSource(Matrix& _Matrix, const Field& src, Field& src_e, Field& src_o)
      {
        GridBase* grid  = _Matrix.RedBlackGrid();
        GridBase* fgrid = _Matrix.Grid();

        Field  tmp(grid);
        Field Mtmp(grid);

        pickCheckerboard(Even, src_e, src);
        pickCheckerboard(Odd , src_o, src);
      
        /////////////////////////////////////////////////////
        // src_o = Mdag * (source_o - Moe MeeInv source_e)
        /////////////////////////////////////////////////////
        _Matrix.MooeeInv(src_e, tmp);   assert(   tmp.Checkerboard() == Even );
        _Matrix.Meooe   (tmp, Mtmp);    assert(  Mtmp.Checkerboard() == Odd  );     
        src_o -= Mtmp;                  assert( src_o.Checkerboard() == Odd  );     
      }

      virtual void RedBlackSolution(Matrix& _Matrix, const Field& sol_o, const Field& src_e, Field& sol)
      {
        GridBase* grid  = _Matrix.RedBlackGrid();
        GridBase* fgrid = _Matrix.Grid();

        Field sol_o_i(grid);
        Field     tmp(grid);
        Field   sol_e(grid);

        ////////////////////////////////////////////////
        // MooeeInv due to pecond
        ////////////////////////////////////////////////
        _Matrix.MooeeInv(sol_o, tmp);
        sol_o_i = tmp;

        ///////////////////////////////////////////////////
        // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
        ///////////////////////////////////////////////////
        _Matrix.Meooe(sol_o_i, tmp);    assert(   tmp.Checkerboard() == Even );
        tmp = src_e - tmp;              assert( src_e.Checkerboard() == Even );
        _Matrix.MooeeInv(tmp, sol_e);   assert( sol_e.Checkerboard() == Even );
       
        setCheckerboard(sol, sol_e);    assert(   sol_e.Checkerboard() == Even );
        setCheckerboard(sol, sol_o_i);  assert( sol_o_i.Checkerboard() == Odd  );
      };

      virtual void RedBlackSolve(Matrix& _Matrix, const Field& src_o, Field& sol_o)
      {
        NonHermitianSchurDiagTwoOperator<Matrix,Field> _OpEO(_Matrix);
        this->_HermitianRBSolver(_OpEO, src_o, sol_o);
      };

      virtual void RedBlackSolve(Matrix& _Matrix, const std::vector<Field>& src_o,  std::vector<Field>& sol_o)
      {
        NonHermitianSchurDiagTwoOperator<Matrix,Field> _OpEO(_Matrix);
        this->_HermitianRBSolver(_OpEO, src_o, sol_o); 
      }
  };
}

#endif
