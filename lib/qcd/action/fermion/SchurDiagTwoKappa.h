#if 1
    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: SchurDiagTwoKappa.h

    Copyright (C) 2017

Author: Christoph Lehner
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
#ifndef  _SCHUR_DIAG_TWO_KAPPA_H
#define  _SCHUR_DIAG_TWO_KAPPA_H

namespace Grid {

  // This is specific to (Z)mobius fermions
  template<class Matrix, class Field>
    class KappaSimilarityTransform {
  public:
    INHERIT_IMPL_TYPES(Matrix);
    std::vector<Coeff_t> kappa, kappaDag, kappaInv, kappaInvDag;

    KappaSimilarityTransform (Matrix &zmob) {
      for (int i=0;i<(int)zmob.bs.size();i++) {
	Coeff_t k = 1.0 / ( 2.0 * (zmob.bs[i] *(4 - zmob.M5) + 1.0) );
	kappa.push_back( k );
	kappaDag.push_back( conj(k) );
	kappaInv.push_back( 1.0 / k );
	kappaInvDag.push_back( 1.0 / conj(k) );
      }
    }

  template<typename vobj>
    void sscale(const Lattice<vobj>& in, Lattice<vobj>& out, Coeff_t* s) {
    GridBase *grid=out._grid;
    out.checkerboard = in.checkerboard;
    assert(grid->_simd_layout[0] == 1); // should be fine for ZMobius for now
    int Ls = grid->_rdimensions[0];
    parallel_for(int ss=0;ss<grid->oSites();ss++){
      vobj tmp = s[ss % Ls]*in._odata[ss];
      vstream(out._odata[ss],tmp);
    }
  }

  RealD sscale_norm(const Field& in, Field& out, Coeff_t* s) {
    sscale(in,out,s);
    return norm2(out);
  }

  virtual RealD M       (const Field& in, Field& out) { return sscale_norm(in,out,&kappa[0]);   }
  virtual RealD MDag    (const Field& in, Field& out) { return sscale_norm(in,out,&kappaDag[0]);}
  virtual RealD MInv    (const Field& in, Field& out) { return sscale_norm(in,out,&kappaInv[0]);}
  virtual RealD MInvDag (const Field& in, Field& out) { return sscale_norm(in,out,&kappaInvDag[0]);}

  };

  template<class Matrix,class Field>
    class SchurDiagTwoKappaOperator :  public SchurOperatorBase<Field> {
  public:
    KappaSimilarityTransform<Matrix, Field> _S;
    SchurDiagTwoOperator<Matrix, Field> _Mat;

    SchurDiagTwoKappaOperator (Matrix &Mat): _S(Mat), _Mat(Mat) {};

    virtual  RealD Mpc      (const Field &in, Field &out) {
      Field tmp(in._grid);

      _S.MInv(in,out);
      _Mat.Mpc(out,tmp);
      return _S.M(tmp,out);

    }
    virtual  RealD MpcDag   (const Field &in, Field &out){
      Field tmp(in._grid);

      _S.MDag(in,out);
      _Mat.MpcDag(out,tmp);
      return _S.MInvDag(tmp,out);
    }
  };

#if 0
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Copied from DiagTwoSolve
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class Field> class SchurRedBlackDiagTwoSolve {
  private:
    OperatorFunction<Field> & _HermitianRBSolver;
    int CBfactorise;
  public:

    /////////////////////////////////////////////////////
    // Wrap the usual normal equations Schur trick
    /////////////////////////////////////////////////////
  SchurRedBlackDiagTwoSolve(OperatorFunction<Field> &HermitianRBSolver)  :
     _HermitianRBSolver(HermitianRBSolver) 
    { 
      CBfactorise=0;
    };

    template<class Matrix>
      void operator() (Matrix & _Matrix,const Field &in, Field &out){

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
      _HermitianRBSolver(_HermOpEO,src_o,tmp);  assert(tmp.checkerboard==Odd);
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
      _Matrix.M(out,resid); 
      resid = resid-in;
      RealD ns = norm2(in);
      RealD nr = norm2(resid);

      std::cout<<GridLogMessage << "SchurRedBlackDiagTwoKappa solver true unprec resid "<< std::sqrt(nr/ns) <<" nr "<< nr <<" ns "<<ns << std::endl;
    }     
  };
#endif
namespace QCD{
    //
    // Determinant is det of middle factor
    // This assumes Mee is indept of U.
    //
    //
    template<class Impl>
    class SchurDifferentiableDiagTwo:  public SchurDiagTwoOperator<FermionOperator<Impl>,typename Impl::FermionField> 
      {
      public:
      INHERIT_IMPL_TYPES(Impl);

 	typedef FermionOperator<Impl> Matrix;

	SchurDifferentiableDiagTwo (Matrix &Mat) : SchurDiagTwoOperator<Matrix,FermionField>(Mat) {};
    };
#if 0
    template<class Impl>
    class SchurDifferentiableDiagTwoKappa :  public SchurDiagTwoKappaOperator<FermionOperator<Impl>,typename Impl::FermionField> 
      {
      public:
      INHERIT_IMPL_TYPES(Impl);

 	typedef FermionOperator<Impl> Matrix;

	SchurDifferentiableDiagTwoKappa (Matrix &Mat) : SchurDiagTwoKappaOperator<Matrix,FermionField>(Mat) {};
    };
#endif
}

}

#endif
#endif
