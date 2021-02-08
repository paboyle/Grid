    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/MADWF.h

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
#pragma once

NAMESPACE_BEGIN(Grid);

template <class Fieldi, class Fieldo,IfNotSame<Fieldi,Fieldo> X=0>
inline void convert(const Fieldi &from,Fieldo &to) 
{
  precisionChange(to,from);
}
template <class Fieldi, class Fieldo,IfSame<Fieldi,Fieldo> X=0>
inline void convert(const Fieldi &from,Fieldo &to) 
{
  to=from;
}

struct MADWFinnerIterCallbackBase{
  virtual void operator()(const RealD current_resid){}
  virtual ~MADWFinnerIterCallbackBase(){}
};

template<class Matrixo,class Matrixi,class PVinverter,class SchurSolver, class Guesser> 
class MADWF 
{
 private:
  typedef typename Matrixo::FermionField FermionFieldo;
  typedef typename Matrixi::FermionField FermionFieldi;

  PVinverter  & PauliVillarsSolvero;// For the outer field
  SchurSolver & SchurSolveri;       // For the inner approx field
  Guesser     & Guesseri;           // To deflate the inner approx solves

  Matrixo & Mato;                   // Action object for outer
  Matrixi & Mati;                   // Action object for inner

  RealD target_resid;
  int   maxiter;

  //operator() is called on "callback" at the end of every inner iteration. This allows for example the adjustment of the inner
  //tolerance to speed up subsequent iteration
  MADWFinnerIterCallbackBase* callback;
  
 public:
  MADWF(Matrixo &_Mato,
	Matrixi &_Mati,
	PVinverter &_PauliVillarsSolvero,
	SchurSolver &_SchurSolveri,
	Guesser & _Guesseri,
	RealD resid,
	int _maxiter,
	MADWFinnerIterCallbackBase* _callback = NULL) :

  Mato(_Mato),Mati(_Mati),
    SchurSolveri(_SchurSolveri),
    PauliVillarsSolvero(_PauliVillarsSolvero),Guesseri(_Guesseri),
    callback(_callback)
    {
      target_resid=resid;
      maxiter     =_maxiter;
    };
   
  void operator() (const FermionFieldo &src,FermionFieldo &sol5)
  {
    std::cout << GridLogMessage<< " ************************************************" << std::endl;
    std::cout << GridLogMessage<< "  MADWF-like algorithm                           " << std::endl;
    std::cout << GridLogMessage<< " ************************************************" << std::endl;

    FermionFieldi    c0i(Mati.GaugeGrid()); // 4d 
    FermionFieldi    y0i(Mati.GaugeGrid()); // 4d
    FermionFieldo    c0 (Mato.GaugeGrid()); // 4d 
    FermionFieldo    y0 (Mato.GaugeGrid()); // 4d

    FermionFieldo    A(Mato.FermionGrid()); // Temporary outer
    FermionFieldo    B(Mato.FermionGrid()); // Temporary outer
    FermionFieldo    b(Mato.FermionGrid()); // 5d source

    FermionFieldo    c(Mato.FermionGrid()); // PVinv source; reused so store
    FermionFieldo    defect(Mato.FermionGrid()); // 5d source

    FermionFieldi   ci(Mati.FermionGrid()); 
    FermionFieldi   yi(Mati.FermionGrid()); 
    FermionFieldi   xi(Mati.FermionGrid()); 
    FermionFieldi srci(Mati.FermionGrid()); 
    FermionFieldi   Ai(Mati.FermionGrid()); 

    RealD m=Mati.Mass();

    ///////////////////////////////////////
    //Import source, include Dminus factors
    ///////////////////////////////////////
    GridBase *src_grid = src.Grid();

    assert( (src_grid == Mato.GaugeGrid()) || (src_grid == Mato.FermionGrid()));

    if ( src_grid == Mato.GaugeGrid() ) {
      Mato.ImportPhysicalFermionSource(src,b);
    } else {
      b=src;
    }
    std::cout << GridLogMessage << " src " <<norm2(src)<<std::endl;
    std::cout << GridLogMessage << " b    " <<norm2(b)<<std::endl;

    defect = b;
    sol5=Zero();
    for (int i=0;i<maxiter;i++) {

      ///////////////////////////////////////
      // Set up c0 from current defect
      ///////////////////////////////////////
      PauliVillarsSolvero(Mato,defect,A);
      Mato.Pdag(A,c);
      ExtractSlice(c0, c, 0 , 0);

      ////////////////////////////////////////////////
      // Solve the inner system with surface term c0
      ////////////////////////////////////////////////
      ci = Zero();  
      convert(c0,c0i); // Possible precison change
      InsertSlice(c0i,ci,0, 0);

      // Dwm P y = Dwm x = D(1) P (c0,0,0,0)^T
      Mati.P(ci,Ai);
      Mati.SetMass(1.0);      Mati.M(Ai,srci);      Mati.SetMass(m);
      SchurSolveri(Mati,srci,xi,Guesseri); 
      Mati.Pdag(xi,yi);
      ExtractSlice(y0i, yi, 0 , 0);
      convert(y0i,y0); // Possible precision change

      //////////////////////////////////////
      // Propagate solution back to outer system
      // Build Pdag PV^-1 Dm P [-sol4,c2,c3... cL]
      //////////////////////////////////////
      c0 = - y0;
      InsertSlice(c0, c, 0   , 0);

      /////////////////////////////
      // Reconstruct the bulk solution Pdag PV^-1 Dm P 
      /////////////////////////////
      Mato.P(c,B);
      Mato.M(B,A);
      PauliVillarsSolvero(Mato,A,B);
      Mato.Pdag(B,A);

      //////////////////////////////
      // Reinsert surface prop
      //////////////////////////////
      InsertSlice(y0,A,0,0);

      //////////////////////////////
      // Convert from y back to x 
      //////////////////////////////
      Mato.P(A,B);

      //         sol5' = sol5 + M^-1 defect
      //               = sol5 + M^-1 src - M^-1 M sol5  ...
      sol5 = sol5 + B;
      std::cout << GridLogMessage << "***************************************" <<std::endl;
      std::cout << GridLogMessage << " Sol5 update "<<std::endl;
      std::cout << GridLogMessage << "***************************************" <<std::endl;
      std::cout << GridLogMessage << " Sol5 now "<<norm2(sol5)<<std::endl;
      std::cout << GridLogMessage << " delta    "<<norm2(B)<<std::endl;

       // New defect  = b - M sol5
       Mato.M(sol5,A);
       defect = b - A;

       std::cout << GridLogMessage << " defect   "<<norm2(defect)<<std::endl;

       double resid = ::sqrt(norm2(defect) / norm2(b));
       std::cout << GridLogMessage << "Residual " << i << ": " << resid  << std::endl;
       std::cout << GridLogMessage << "***************************************" <<std::endl;

       if(callback != NULL) (*callback)(resid);       
       
       if (resid < target_resid) {
	 return;
       }
    }

    std::cout << GridLogMessage << "MADWF : Exceeded maxiter "<<std::endl;
    assert(0);

  }

};

NAMESPACE_END(Grid);
