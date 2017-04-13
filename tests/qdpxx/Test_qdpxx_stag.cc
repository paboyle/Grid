    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/qdpxx/Test_qdpxx_munprec.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#include <Grid/Grid.h>

double mq=0.1;

typedef Grid::QCD::StaggeredImplR::FermionField FermionField;
typedef Grid::QCD::LatticeGaugeField GaugeField;

void make_gauge     (GaugeField & lat, FermionField &src);
void calc_grid      (GaugeField & lat, GaugeField & uthin,GaugeField & ufat, FermionField &src, FermionField &res,int dag);
void calc_chroma    (GaugeField & lat,GaugeField & uthin,GaugeField & ufat, FermionField &src, FermionField &res,int dag);

#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <actions/ferm/invert/syssolver_linop_aggregate.h>

namespace Chroma { 


class ChromaWrapper {
public:
  
  typedef multi1d<LatticeColorMatrix> U;
  typedef LatticeStaggeredFermion T4;
  
  static void ImportGauge(GaugeField & gr,
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

  static void ExportGauge(GaugeField & gr,
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

      for(int mu=0;mu<4;mu++){
	for(int i=0;i<3;i++){
	for(int j=0;j<3;j++){
	  cm = QDP::peekSite(ch[mu],cx);
	  c  = QDP::peekColor(cm,i,j);
	  cc = Grid::Complex(toDouble(real(c)),toDouble(imag(c)));
	  LCM(mu)()(i,j)= cc;
	}}
      }
      Grid::pokeSite(LCM,gr,x);

    }}}}
  }

  
  static void ImportFermion(FermionField & gr,
			    QDP::LatticeStaggeredFermion & ch  ) 
  {
    Grid::QCD::ColourVector F;
    Grid::Complex c;


    std::vector<int> x(5);
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

      Grid::peekSite(F,gr,x);
      QDP::ColorVector cv;
      for(int j=0;j<3;j++){
	QDP::Complex cc;
	c  = F()()(j) ;
	cc = QDP::cmplx(QDP::Real(real(c)),QDP::Real(imag(c)));
	pokeColor(cv,cc,j);
      }
      QDP::StaggeredFermion cF;
      pokeSpin(cF,cv,0);
      QDP::pokeSite(ch,cF,cx);
    }}}}
  }
  static void ExportFermion(FermionField & gr,
			    QDP::LatticeStaggeredFermion & ch  ) 
  {
    Grid::QCD::ColourVector F;
    Grid::Complex c;

    std::vector<int> x(5);
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

      QDP::StaggeredFermion cF = QDP::peekSite(ch,cx);
      for(int j=0;j<3;j++){
	QDP::ColorVector cS=QDP::peekSpin(cF,0);
	QDP::Complex cc=QDP::peekColor(cS,j);
	c = Grid::Complex(QDP::toDouble(QDP::real(cc)), 
			  QDP::toDouble(QDP::imag(cc)));
	F()()(j) = c;
      }
      Grid::pokeSite(F,gr,x);
    }}}}
  }

  static Handle< Chroma::EvenOddLinearOperator<T4,U,U> >  GetLinOp (U &u,U &u_fat,U &u_triple)
  {
    QDP::Real _mq(mq);
    QDP::multi1d<int> bcs(QDP::Nd);

    bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;

    Chroma::AsqtadFermActParams p; 
    p.Mass = _mq; 
    p.u0 = Real(1.0);


    Chroma::Handle<Chroma::FermBC<T4,U,U> > fbc(new Chroma::SimpleFermBC< T4, U, U >(bcs));
    Chroma::Handle<Chroma::CreateFermState<T4,U,U> > cfs( new Chroma::CreateSimpleFermState<T4,U,U>(fbc));
    Chroma::AsqtadFermAct S_f(cfs,p);
    Chroma::Handle< Chroma::FermState<T4,U,U> >  ffs(  S_f.createState(u) );
    u_fat   =ffs.cast<AsqtadConnectStateBase>()->getFatLinks();
    u_triple=ffs.cast<AsqtadConnectStateBase>()->getTripleLinks();
    return S_f.linOp(ffs);
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

  GaugeField uthin  (UGrid);
  GaugeField ufat   (UGrid);
  GaugeField utriple(UGrid);
  FermionField    src(UGrid);
  FermionField    res_chroma(UGrid);
  FermionField    res_grid  (UGrid);
  

  {

    std::cout << "*****************************"<<std::endl;
    std::cout << "Staggered Action "            <<std::endl;
    std::cout << "*****************************"<<std::endl;

    make_gauge(uthin,src);

    for(int dag=0;dag<2;dag++) {

      std::cout << "Dag =  "<<dag<<std::endl;
      
      calc_chroma(uthin,utriple,ufat,src,res_chroma,dag);

      // Remove the normalisation of Chroma Gauge links ??
      std::cout << "Norm of chroma Asqtad multiply "<<Grid::norm2(res_chroma)<<std::endl;
      calc_grid  (uthin,utriple,ufat,src,res_grid,dag);

      std::cout << "Norm of thin gauge "<< Grid::norm2(uthin) <<std::endl;
      std::cout << "Norm of fat  gauge "<< Grid::norm2(ufat) <<std::endl;

      std::cout << "Norm of Grid Asqtad multiply "<<Grid::norm2(res_grid)<<std::endl;
      
      /*
      std::cout << " site 0 of Uthin  "<<uthin._odata[0] <<std::endl;
      std::cout << " site 0 of Utriple"<<utriple._odata[0] <<std::endl;
      std::cout << " site 0 of Ufat   "<<ufat._odata[0] <<std::endl;

      std::cout << " site 0 of Grid   "<<res_grid._odata[0] <<std::endl;
      std::cout << " site 0 of Chroma "<<res_chroma._odata[0] <<std::endl;
      */

      res_chroma=res_chroma - res_grid;
      std::cout << "Norm of difference "<<Grid::norm2(res_chroma)<<std::endl;
    }
  }

  std::cout << "Finished test "<<std::endl;

  Chroma::finalize();
}

void calc_chroma(GaugeField & lat, GaugeField &uthin, GaugeField &ufat, FermionField &src, FermionField &res,int dag)
{
  typedef QDP::LatticeStaggeredFermion T;
  typedef QDP::multi1d<QDP::LatticeColorMatrix> U;
  
  U u(4);
  U ut(4);
  U uf(4);

  //  Chroma::HotSt(u);
  Chroma::ChromaWrapper::ImportGauge(lat,u) ;

  QDP::LatticeStaggeredFermion  check;
  QDP::LatticeStaggeredFermion  result;
  QDP::LatticeStaggeredFermion  tmp;
  QDP::LatticeStaggeredFermion  psi;

  Chroma::ChromaWrapper::ImportFermion(src,psi);

  auto linop =Chroma::ChromaWrapper::GetLinOp(u,uf,ut);

  Chroma::ChromaWrapper::ExportGauge(uthin,ut) ;
  Chroma::ChromaWrapper::ExportGauge(ufat ,uf) ;

  enum Chroma::PlusMinus isign;
  if ( dag ) {
    isign=Chroma::MINUS;
  } else {
    isign=Chroma::PLUS;
  }

  std::cout << "Calling Chroma Linop "<< std::endl;
  linop->evenEvenLinOp(tmp,psi,isign); check[rb[0]] = tmp;
  linop->oddOddLinOp  (tmp,psi,isign); check[rb[1]] = tmp;
  linop->evenOddLinOp(tmp,psi,isign) ; check[rb[0]]+= tmp;
  linop->oddEvenLinOp(tmp,psi,isign) ; check[rb[1]]+= tmp;

  Chroma::ChromaWrapper::ExportFermion(res,check) ;
}


void make_gauge(GaugeField & Umu,FermionField &src)
{
  using namespace Grid;
  using namespace Grid::QCD;

  std::vector<int> seeds4({1,2,3,4});

  Grid::GridCartesian         * UGrid   = (Grid::GridCartesian *) Umu._grid;
  Grid::GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  Grid::QCD::SU3::HotConfiguration(RNG4,Umu);
  Grid::gaussian(RNG4,src);
}

void calc_grid(GaugeField & Uthin, GaugeField & Utriple, GaugeField & Ufat, FermionField &src, FermionField &res,int dag)
{
  using namespace Grid;
  using namespace Grid::QCD;

  Grid::GridCartesian         * UGrid   = (Grid::GridCartesian *) Uthin._grid;
  Grid::GridRedBlackCartesian * UrbGrid = Grid::QCD::SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  Grid::QCD::ImprovedStaggeredFermionR Dstag(Uthin,Utriple,Ufat,*UGrid,*UrbGrid,mq*2.0);

  std::cout << Grid::GridLogMessage <<" Calling Grid staggered multiply "<<std::endl;

  if ( dag ) 
    Dstag.Mdag(src,res);  
  else 
    Dstag.M(src,res);  

  res = res ; // Convention mismatch to Chroma
  return;
} 





