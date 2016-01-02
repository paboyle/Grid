    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/qdpxx/Test_qdpxx_loops_staples.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>

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
#include <Grid.h>

double calc_grid_p      (Grid::QCD::LatticeGaugeField & lat);
double calc_chroma_p    (Grid::QCD::LatticeGaugeField & lat);
double calc_grid_r      (Grid::QCD::LatticeGaugeField & lat);
double calc_grid_IW     (Grid::QCD::LatticeGaugeField & lat);
double calc_grid_r_dir  (Grid::QCD::LatticeGaugeField & lat);
double calc_chroma_r    (Grid::QCD::LatticeGaugeField & lat);
double calc_chroma_IW   (Grid::QCD::LatticeGaugeField & lat);
void check_grid_r_staple(Grid::QCD::LatticeGaugeField & Umu);
void check_grid_p_staple(Grid::QCD::LatticeGaugeField & Umu);

const double beta=2.6;
const double c1=-0.331;
#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <actions/ferm/invert/syssolver_linop_aggregate.h>



namespace Chroma { 

class ChromaWrapper {
public:

  
  typedef multi1d<LatticeColorMatrix> U;
  
  static void ImportGauge(Grid::QCD::LatticeGaugeField & gr,
			  QDP::multi1d<QDP::LatticeColorMatrix> & ch) 
  {
    Grid::QCD::LorentzColourMatrix LCM;
    Grid::Complex cc;
    QDP::ColorMatrix cm;
    QDP::Complex c;

    std::vector<int> x(4);
    QDP::multi1d<int> cx(4);
    std::vector<int> gd= gr._grid->GlobalDimensions();

    for (x[0]=0;x[0]<gd[0];x[0]++){
    for (x[1]=0;x[1]<gd[1];x[1]++){
    for (x[2]=0;x[2]<gd[2];x[2]++){
    for (x[3]=0;x[3]<gd[3];x[3]++){
      cx[0] = x[0];
      cx[1] = x[1];
      cx[2] = x[2];
      cx[3] = x[3];
      Grid::peekSite(LCM,gr,x);

      for(int mu=0;mu<4;mu++){
	for(int i=0;i<3;i++){
	for(int j=0;j<3;j++){
	  cc = LCM(mu)()(i,j);
	  c = QDP::cmplx(QDP::Real(real(cc)),QDP::Real(imag(cc)));
	  QDP::pokeColor(cm,c,i,j);
	}}
	QDP::pokeSite(ch[mu],cm,cx);
      }

    }}}}
  }

  static Chroma::Handle< Chroma::LinearGaugeAction > GetRBCAction ( U &u )
  {
    Chroma::CreatePeriodicGaugeState<U,U> CPGS;
    Chroma::Handle< Chroma::CreateGaugeState<U,U> > cgs (new Chroma::CreatePeriodicGaugeState<U,U>);
    Chroma::Handle< Chroma::LinearGaugeAction> action = new Chroma::RBCGaugeAct(cgs,beta,c1);
    return action;
  }

  static Chroma::Handle< Chroma::LinearGaugeAction > GetRectangle ( U &u )
  {
    Chroma::CreatePeriodicGaugeState<U,U> CPGS;
    Chroma::Handle< Chroma::CreateGaugeState<U,U> > cgs (new Chroma::CreatePeriodicGaugeState<U,U>);
    Chroma::Handle< Chroma::LinearGaugeAction> action = new Chroma::RectGaugeAct(cgs, Real(c1));
    return action;
  }

  static Chroma::Handle< Chroma::LinearGaugeAction > GetPlaquette ( U &u )
  {
    Chroma::CreatePeriodicGaugeState<U,U> CPGS;
    Chroma::Handle< Chroma::CreateGaugeState<U,U> > cgs (new Chroma::CreatePeriodicGaugeState<U,U>);
    Chroma::Handle< Chroma::LinearGaugeAction> action = new Chroma::WilsonGaugeAct(cgs, Real(beta));
    return action;
  }

};
}

int main (int argc,char **argv )
{
  /********************************************************
   * Setup QDP
   *********************************************************/
  Chroma::initialize(&argc,&argv);
  Chroma::WilsonTypeFermActs4DEnv::registerAll(); 

  /********************************************************
   * Setup Grid
   *********************************************************/
  Grid::Grid_init(&argc,&argv);
  Grid::GridCartesian * UGrid   = Grid::QCD::SpaceTimeGrid::makeFourDimGrid(Grid::GridDefaultLatt(), 
									    Grid::GridDefaultSimd(Grid::QCD::Nd,Grid::vComplex::Nsimd()),
									    Grid::GridDefaultMpi());
  
  std::vector<int> gd = UGrid->GlobalDimensions();
  QDP::multi1d<int> nrow(QDP::Nd);
  for(int mu=0;mu<4;mu++) nrow[mu] = gd[mu];

  QDP::Layout::setLattSize(nrow);
  QDP::Layout::create();

  Grid::QCD::LatticeGaugeField lat(UGrid);

  double s_grid   = calc_grid_p  (lat);

  double s_chroma = calc_chroma_p(lat);
  // Match conventions
  double vol = UGrid->gSites();
  s_chroma+= beta * QDP::Nd * (QDP::Nd-1)*vol *0.5;

  std::cout << " Chroma/Grid plaquette = " <<s_chroma<<" "<<s_grid<<std::endl;
  std::cout << " Chroma avg plaquette = "  <<1.0 - s_chroma / vol/ 6 / beta <<std::endl;

  s_grid   = calc_grid_r  (lat);
  s_chroma = calc_chroma_r(lat); 
  std::cout << std::setprecision(10);
  std::cout<< " bare s_chroma "<<s_chroma<<std::endl;
  s_chroma+= c1 * 12.0 *vol ;
  std::cout<< " adjusted s_chroma "<<s_chroma<<std::endl;
  std::cout<< " adjust "<< c1 * 12.0 *vol <<std::endl;

  std::cout << " Chroma/Grid rectangle = " <<s_chroma<<" "<<s_grid<<std::endl;
  std::cout << " Chroma avg  rectangle = "  <<1.0 - s_chroma / vol/ 12.0 / c1 <<std::endl;

  check_grid_r_staple(lat);

  check_grid_p_staple(lat);

  calc_grid_r_dir(lat);
  
  // Iwasaki case
  std::cout << "Checking Iwasaki action (Grid/Chroma) " << std::endl;
  s_grid   = calc_grid_IW(lat);
  s_chroma = calc_chroma_IW(lat);

  // Adjust for the missing "1" in chroma's action.
  s_chroma+= beta * (1.0-8*c1) * QDP::Nd * (QDP::Nd-1)*vol *0.5;
  s_chroma+= beta * c1 * 12.0 *vol ;

  std::cout << "Iwasaki action (Grid/Chroma) = " << s_grid<< " " << s_chroma <<std::endl;

  Chroma::finalize();
}

double calc_chroma_IW(Grid::QCD::LatticeGaugeField & lat)
{
  typedef QDP::multi1d<QDP::LatticeColorMatrix> U;

  QDP::multi1d<QDP::LatticeColorMatrix> u(4);

  Chroma::ChromaWrapper::ImportGauge(lat,u) ;

  for(int mu=0;mu<4;mu++){
    std::cout <<"Imported Gauge norm ["<<mu<<"] "<< QDP::norm2(u[mu])<<std::endl;
  }

  auto action =Chroma::ChromaWrapper::GetRBCAction(u);

  Chroma::Handle<GaugeState<U,U> > gs(action->getCreateState()(u));

  Real act = action->S(gs);
  
  double s = toDouble(act);

  return s;
}
double calc_chroma_r(Grid::QCD::LatticeGaugeField & lat)
{
  typedef QDP::multi1d<QDP::LatticeColorMatrix> U;

  QDP::multi1d<QDP::LatticeColorMatrix> u(4);

  Chroma::ChromaWrapper::ImportGauge(lat,u) ;

  for(int mu=0;mu<4;mu++){
    std::cout <<"Imported Gauge norm ["<<mu<<"] "<< QDP::norm2(u[mu])<<std::endl;
  }

  auto action =Chroma::ChromaWrapper::GetRectangle(u);

  Chroma::Handle<GaugeState<U,U> > gs(action->getCreateState()(u));

  Real act = action->S(gs);
  
  double s = toDouble(act);

  return s;
}

// Conventions matching:
//
// Chroma:
//
// w_plaq is defined in MesPlq as
// w_plaq =( 2/(V*Nd*(Nd-1)*Nc)) * Sum Re Tr Plaq
//
// S = -(coeff/(Nc) Sum Re Tr Plaq
//   
// S = -coeff * (V*Nd*(Nd-1)/2) w_plaq 
//   = -coeff * (V*Nd*(Nd-1)/2)*(2/(V*Nd*(Nd-1)*Nc))* Sum Re Tr Plaq
//   = -coeff * (1/(Nc)) * Sum Re Tr Plaq
//
// Grid: has 1-plaq as usual to define Fmunu
//
// action = beta*(1.0 -plaq)*(Nd*(Nd-1.0))*vol*0.5;
// action = beta * Nd*Nd-1*vol*0.5 - beta * Nd*Nd-1*vol*0.5*plaq
//
// plaq == sumplaq * 2/(V*Nd*(Nd-1)*Nc)
double calc_chroma_p(Grid::QCD::LatticeGaugeField & lat)
{
  typedef QDP::multi1d<QDP::LatticeColorMatrix> U;

  QDP::multi1d<QDP::LatticeColorMatrix> u(4);

  Chroma::ChromaWrapper::ImportGauge(lat,u) ;

  for(int mu=0;mu<4;mu++){
    std::cout <<"Imported Gauge norm ["<<mu<<"] "<< QDP::norm2(u[mu])<<std::endl;
  }

  auto action =Chroma::ChromaWrapper::GetPlaquette(u);

  Chroma::Handle<GaugeState<U,U> > gs(action->getCreateState()(u));

  Real act = action->S(gs);

  double s = toDouble(act);

  return s;
}



double calc_grid_p(Grid::QCD::LatticeGaugeField & Umu)
{
  std::vector<int> seeds4({1,2,3,4});

  Grid::GridCartesian         * UGrid   = (Grid::GridCartesian *) Umu._grid;
  Grid::GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  Grid::QCD::SU3::HotConfiguration(RNG4,Umu);

  Grid::QCD::LatticeColourMatrix tmp(UGrid); 
  tmp = Grid::zero;

  Grid::QCD::PokeIndex<Grid::QCD::LorentzIndex>(Umu,tmp,2);
  Grid::QCD::PokeIndex<Grid::QCD::LorentzIndex>(Umu,tmp,3);

  Grid::QCD::WilsonGaugeActionR Wilson(beta); // Just take beta = 1.0
  
  return Wilson.S(Umu);
} 
double calc_grid_r(Grid::QCD::LatticeGaugeField & Umu)
{
  Grid::GridCartesian         * UGrid   = (Grid::GridCartesian *) Umu._grid;

  Grid::QCD::PlaqPlusRectangleActionR Wilson(0.0,c1); // Just take beta = 0.0
  
  return Wilson.S(Umu);
} 
double calc_grid_IW(Grid::QCD::LatticeGaugeField & Umu)
{
  Grid::GridCartesian         * UGrid   = (Grid::GridCartesian *) Umu._grid;

  Grid::QCD::IwasakiGaugeActionR Iwasaki(beta);
  
  return Iwasaki.S(Umu);
} 
double calc_grid_r_dir(Grid::QCD::LatticeGaugeField & Umu)
{
  Grid::GridCartesian         * UGrid   = (Grid::GridCartesian *) Umu._grid;

  std::vector<Grid::QCD::LatticeColourMatrix> U(4,UGrid);
  for(int mu=0;mu<Nd;mu++){
    U[mu] = Grid::PeekIndex<Grid::QCD::LorentzIndex>(Umu,mu);
  }

  Grid::QCD::LatticeComplex rect(UGrid);
  Grid::QCD::TComplex trect;
  Grid::QCD::Complex  crect;
  Grid::RealD vol = UGrid->gSites();
  for(int mu=0;mu<Grid::QCD::Nd;mu++){
  for(int nu=0;nu<Grid::QCD::Nd;nu++){
    if ( mu!=nu ) {

      Grid::QCD::WilsonLoops<Grid::QCD::LatticeGaugeField>::traceDirRectangle(rect,U,mu,nu);
      trect = Grid::sum(rect);
      crect = Grid::TensorRemove(trect);
      std::cout<< "mu/nu = "<<mu<<"/"<<nu<<" ; rect = "<<crect/vol/2.0/3.0<<std::endl;

      Grid::GridStopWatch Peter;
      Grid::GridStopWatch Azusa;

      // Staple test
      Peter.Start();
      {                                                  
	Grid::QCD::LatticeColourMatrix Stap(UGrid);
	Grid::QCD::LatticeComplex      SumTrStap(UGrid);
	Grid::QCD::LatticeComplex      TrStap(UGrid);

	/*
	 * Make staple for loops centered at coor of link ; this one is ok.     //     |
	 */
	
	//           __ ___ 
	//          |    __ |
	Stap = 
	  Grid::Cshift(Grid::QCD::CovShiftForward (U[mu],mu,
		       Grid::QCD::CovShiftForward (U[nu],nu,
		       Grid::QCD::CovShiftBackward(U[mu],mu,
                       Grid::QCD::CovShiftBackward(U[mu],mu,
		       Grid::Cshift(adj(U[nu]),nu,-1))))) , mu, 1);

	TrStap = Grid::trace (U[mu]*Stap);
	SumTrStap = TrStap;

	trect = Grid::sum(TrStap);
	crect = Grid::TensorRemove(trect);
	std::cout<< "mu/nu inline trace 2x1 code = "<<mu<<"/"<<nu<<" ; rect = "<<crect/vol/1.0/3.0<<std::endl;


	//              __ 
	//          |__ __ |

	Stap = Grid::Cshift(Grid::QCD::CovShiftForward (U[mu],mu,
		            Grid::QCD::CovShiftBackward(U[nu],nu,
   		            Grid::QCD::CovShiftBackward(U[mu],mu,
                            Grid::QCD::CovShiftBackward(U[mu],mu, U[nu])))) , mu, 1);

	TrStap = Grid::trace (U[mu]*Stap);

	trect = Grid::sum(TrStap);
	crect = Grid::TensorRemove(trect);
	std::cout<< "mu/nu inline trace 2x1 code = "<<mu<<"/"<<nu<<" ; rect = "<<crect/vol/1.0/3.0<<std::endl;

	//           __ 
	//          |__ __ |

	Stap = Grid::Cshift(Grid::QCD::CovShiftBackward(U[nu],nu,
		            Grid::QCD::CovShiftBackward(U[mu],mu,
                            Grid::QCD::CovShiftBackward(U[mu],mu,
   		            Grid::QCD::CovShiftForward(U[nu],nu,U[mu])))) , mu, 1);

	TrStap = Grid::trace (U[mu]*Stap);

	trect = Grid::sum(TrStap);
	crect = Grid::TensorRemove(trect);
	std::cout<< "mu/nu inline trace 2x1 code = "<<mu<<"/"<<nu<<" ; rect = "<<crect/vol/1.0/3.0<<std::endl;


	//           __ ___ 
	//          |__    |

	Stap = Grid::Cshift(Grid::QCD::CovShiftForward (U[nu],nu,
		            Grid::QCD::CovShiftBackward(U[mu],mu,
                            Grid::QCD::CovShiftBackward(U[mu],mu,
                            Grid::QCD::CovShiftBackward(U[nu],nu,U[mu])))) , mu, 1);


	TrStap = Grid::trace (U[mu]*Stap);
	trect = Grid::sum(TrStap);
	crect = Grid::TensorRemove(trect);
	std::cout<< "mu/nu inline trace 2x1 code = "<<mu<<"/"<<nu<<" ; rect = "<<crect/vol/1.0/3.0<<std::endl;

  
	//       --    
	//      |  |              
	//          
	//      |  | 
	//       

	/*
	 * Make staple for loops centered at coor of link ; this one is ok.     //     |
	 */
	//	Stap = 
	//	  Grid::Cshift(Grid::QCD::CovShiftForward(U[nu],nu,U[nu]),mu,1)* // ->||
	//	  Grid::adj(Grid::QCD::CovShiftForward(U[nu],nu,Grid::QCD::CovShiftForward(U[nu],nu,U[mu]))) ;
	Stap = Grid::Cshift(Grid::QCD::CovShiftForward(U[nu],nu,
		            Grid::QCD::CovShiftForward(U[nu],nu,
                            Grid::QCD::CovShiftBackward(U[mu],mu,
                            Grid::QCD::CovShiftBackward(U[nu],nu,  Grid::Cshift(adj(U[nu]),nu,-1))))) , mu, 1);
	  
	TrStap = Grid::trace (U[mu]*Stap);
	SumTrStap += TrStap;
	trect = Grid::sum(TrStap);
	crect = Grid::TensorRemove(trect);
	std::cout<< "mu/nu inline trace 1x2 code = "<<mu<<"/"<<nu<<" ; rect = "<<crect/vol/1.0/3.0<<std::endl;
  


	//          
	//      |  |              
	//          
	//      |  | 
	//       -- 

	Stap = Grid::Cshift(Grid::QCD::CovShiftBackward(U[nu],nu,
		            Grid::QCD::CovShiftBackward(U[nu],nu,
                            Grid::QCD::CovShiftBackward(U[mu],mu,
                            Grid::QCD::CovShiftForward (U[nu],nu,U[nu])))) , mu, 1);

	TrStap = Grid::trace (U[mu]*Stap);
	trect = Grid::sum(TrStap);
	crect = Grid::TensorRemove(trect);
	std::cout<< "mu/nu inline trace 1x2 code = "<<mu<<"/"<<nu<<" ; rect = "<<crect/vol/1.0/3.0<<std::endl;

	trect = Grid::sum(SumTrStap);
	crect = Grid::TensorRemove(trect);
	std::cout<< "mu/nu inline trace 2x1+1x2 code = "<<mu<<"/"<<nu<<" ; rect = "<<crect/vol/2.0/3.0<<std::endl;
      }
      Peter.Stop();
      Azusa.Start();
      {
	Grid::QCD::LatticeComplex RectPlaq_d(UGrid);
	Grid::QCD::LatticeColourMatrix ds_U(UGrid);
	Grid::QCD::LatticeColourMatrix left_2(UGrid);
	Grid::QCD::LatticeColourMatrix upper_l(UGrid);
	Grid::QCD::LatticeColourMatrix upper_staple(UGrid);
	Grid::QCD::LatticeColourMatrix down_l(UGrid);
	Grid::QCD::LatticeColourMatrix down_staple(UGrid);
	Grid::QCD::LatticeColourMatrix tmp(UGrid);
	
	// 2 (mu)x1(nu)
	left_2=  Grid::QCD::CovShiftForward(U[mu],mu,U[mu]);   // Umu(x) Umu(x+mu)
	tmp=Grid::Cshift(U[nu],mu,2);                          // Unu(x+2mu)

	upper_l=  Grid::QCD::CovShiftForward(tmp,nu,Grid::adj(left_2)); //  Unu(x+2mu) Umu^dag(x+mu+nu) Umu^dag(x+nu) 
	//                 __ __ 
	//              =       |
	
	// Unu(x-2mu) Umudag(x-mu+nu) Umudag(x+nu) Unudag(x)
	//                 __ __ 
	// upper_staple = |     |
	//                      v

	upper_staple= upper_l*adj(U[nu]); //  Unu(x+2mu) Umu^dag(x+mu+nu) Umu^dag(x+nu)  Unu^dag(x)

	//                 
	// down_staple  = |__ __|
	//
	tmp = adj(left_2*tmp)*U[nu];
	down_staple= Grid::Cshift(tmp,nu,-1); // Unu^dag((x+2mu-nu)  Umu^dag(x+mu-nu) Umu^dag(x-nu) Unu(x-nu)

	//  __  __
	// |    __|
	//
	tmp=Grid::Cshift(U[mu],mu,1);
	ds_U=tmp*(upper_staple);           // Umu(x+mu) Unu(x+2mu) Umu^dag(x+mu+nu) Umu^dag(x+nu)  Unu^dag(x)

	RectPlaq_d = Grid::trace(U[mu]*ds_U);
	trect = Grid::sum(RectPlaq_d);
	crect = Grid::TensorRemove(trect);
	std::cout<< "mu/nu inline AZUSA trace 2x1 code = "<<mu<<"/"<<nu<<" ; rect = "<<crect/vol/1.0/3.0<<std::endl;

	//  __  __
	// |__    |
	//

	tmp=upper_staple*U[mu];
	ds_U= Grid::Cshift(tmp,mu,-1);

	RectPlaq_d = Grid::trace(U[mu]*ds_U);
	trect = Grid::sum(RectPlaq_d);
	crect = Grid::TensorRemove(trect);
	std::cout<< "mu/nu inline AZUSA trace 2x1 code = "<<mu<<"/"<<nu<<" ; rect = "<<crect/vol/1.0/3.0<<std::endl;

	//      __
	//  |__ __ |
	//

	tmp=Grid::Cshift(U[mu],mu,1);
	ds_U=tmp*(down_staple);           // Umu(x+mu) Unu^dag((x+2mu-nu)  Umu^dag(x+mu-nu) Umu^dag(x-nu) Unu(x-nu) 

	RectPlaq_d = Grid::trace(U[mu]*ds_U);
	trect = Grid::sum(RectPlaq_d);
	crect = Grid::TensorRemove(trect);
	std::cout<< "mu/nu inline AZUSA trace 2x1 code = "<<mu<<"/"<<nu<<" ; rect = "<<crect/vol/1.0/3.0<<std::endl;


	//   __
	//  |__ __ |
	//

	tmp = down_staple*U[mu];
	ds_U=Grid::Cshift(tmp,mu,-1);

	RectPlaq_d = Grid::trace(U[mu]*ds_U);
	trect = Grid::sum(RectPlaq_d);
	crect = Grid::TensorRemove(trect);
	std::cout<< "mu/nu inline AZUSA trace 2x1 code = "<<mu<<"/"<<nu<<" ; rect = "<<crect/vol/1.0/3.0<<std::endl;

	   
	// 1(mu) x 2 (nu)  ** this was ok
	//   _
	//  | |
	//  | |
	Grid::QCD::LatticeColourMatrix up2= Grid::QCD::CovShiftForward(U[nu],nu,U[nu]);

	upper_l= Grid::QCD::CovShiftForward(Grid::Cshift(up2,mu,1),nu,Grid::Cshift(adj(U[mu]),nu,1));
	ds_U= upper_l*Grid::adj(up2);

	RectPlaq_d = Grid::trace(U[mu]*ds_U);
	trect = Grid::sum(RectPlaq_d);
	crect = Grid::TensorRemove(trect);

	std::cout<< "mu/nu inline AZUSA trace 1x2 code = "<<mu<<"/"<<nu<<" ; rect = "<<crect/vol/1.0/3.0<<std::endl;

	// 1(mu) x 2 (nu) ** this was ok
	//   
	//  | |
	//  |_|

/*
                       |
                       V 
   downer_l=           |  
               (x)<----V                 
*/    
	down_l= Grid::adj(Grid::QCD::CovShiftForward(U[mu],mu,up2)); //downer_l
/*
                     ^     |
   down_staple  =    |     V 
                     ^     | 
                     |     V
                  (x)<----
          down_staple= upper*upper_l;
*/    
	tmp= down_l*up2;
	ds_U= Grid::Cshift(tmp,nu,-2);

	RectPlaq_d = Grid::trace(U[mu]*ds_U);
	trect = Grid::sum(RectPlaq_d);
	crect = Grid::TensorRemove(trect);

	std::cout<< "mu/nu inline AZUSA trace 1x2 code = "<<mu<<"/"<<nu<<" ; rect = "<<crect/vol/1.0/3.0<<std::endl;
	
      }
      Azusa.Stop();

      std::chrono::milliseconds pt = Peter.Elapsed();
      std::chrono::milliseconds az = Azusa.Elapsed();

      std::cout<< "Times ";
      std::cout<<pt.count();
      std::cout<<" A ";
      std::cout<<az.count();
      std::cout<<std::endl;

    }
  }}
  Grid::QCD::PlaqPlusRectangleActionR Wilson(0.0,c1); // Just take beta = 0.0
  
  return Wilson.S(Umu);
};

void check_grid_r_staple(Grid::QCD::LatticeGaugeField & Umu)
{

  std::vector<int> seeds4({1,2,3,4});

  Grid::GridCartesian         * UGrid   = (Grid::GridCartesian *) Umu._grid;

  Grid::QCD::PlaqPlusRectangleActionR Wilson(0.0,c1); // Just take beta = 0.0

  Grid::QCD::LatticeColourMatrix staple(UGrid);
  Grid::QCD::LatticeColourMatrix link(UGrid);
  Grid::QCD::LatticeComplex Traced(UGrid);
  Grid::Complex Rplaq(0.0);

  for(int mu=0;mu<Nd;mu++){
    
    Grid::RealD vol = UGrid->gSites();
    
    // for mu, nu!=mu => 12
    // 6 loops contribute to staple for each orientation.
    // Nc=3.
    // Vol as for each site
    Grid::RealD RectScale(1.0/vol/12.0/6.0/3.0); 

    Grid::QCD::WilsonLoops<Grid::QCD::LatticeGaugeField>::RectStaple(staple,Umu,mu);
    
    link = Grid::QCD::PeekIndex<Grid::QCD::LorentzIndex>(Umu,mu);

    Traced = Grid::trace( link*staple) * RectScale;
    Grid::QCD::TComplex Tp = Grid::sum(Traced);
    Grid::Complex p   = Grid::TensorRemove(Tp);

    std::cout<< "Rect from RectStaple "<<p<<std::endl;
    Rplaq = Rplaq+ p;

  }
  std::cout<< "Rect from RectStaple "<<Rplaq<<std::endl;
} 

void check_grid_p_staple(Grid::QCD::LatticeGaugeField & Umu)
{

  std::vector<int> seeds4({1,2,3,4});

  Grid::GridCartesian         * UGrid   = (Grid::GridCartesian *) Umu._grid;

  Grid::QCD::PlaqPlusRectangleActionR Wilson(1.0,0.0); // Just take c1 = 0.0

  Grid::QCD::LatticeColourMatrix staple(UGrid);
  Grid::QCD::LatticeColourMatrix link(UGrid);
  Grid::QCD::LatticeComplex Traced(UGrid);
  Grid::Complex plaq(0.0);

  for(int mu=0;mu<Nd;mu++){
    
    Grid::RealD vol = UGrid->gSites();
    
    // for mu, nu!=mu => 12
    // 2 loops contribute to staple for each orientation.
    // Nc=3.
    // Vol as for each site
    Grid::RealD Scale(1.0/vol/12.0/2.0/3.0); 

    Grid::QCD::WilsonLoops<Grid::QCD::LatticeGaugeField>::Staple(staple,Umu,mu);
    
    link = Grid::QCD::PeekIndex<Grid::QCD::LorentzIndex>(Umu,mu);

    Traced = Grid::trace( link*staple) * Scale;
    Grid::QCD::TComplex Tp = Grid::sum(Traced);
    Grid::Complex p   = Grid::TensorRemove(Tp);

    std::cout<< "Plaq from PlaqStaple "<<p<<std::endl;
    plaq = plaq+ p;

  }
  std::cout<< "Plaq from PlaqStaple "<<plaq<<std::endl;
} 




