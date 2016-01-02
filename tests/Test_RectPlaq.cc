    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_RectPlaq.cc

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
#include <Grid.h>

#include <qcd/utils/CovariantCshift.h>
#include <qcd/utils/WilsonLoops.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


/* For Metropolis */
class Metropolis {
public:
  GridSerialRNG & sRNG;
  Metropolis(GridSerialRNG & _sRNG) : sRNG(_sRNG) {};
  bool AcceptReject(const RealD Delta)
  {
    RealD rand;

    if(Delta <=0.0) return true;

    random(sRNG,rand);
    if(rand <= exp(-Delta))
      return true;
    else 
      return false;
  }
};

void RectPlaq(const std::vector<LatticeColourMatrix> &U, LatticeComplex &RectPlaqValue )
{
  RectPlaqValue=zero;
  // 12 * vol loops
  for(int mu=1;mu<Nd;mu++){
    for(int nu=0;nu<mu;nu++){
      RectPlaqValue = RectPlaqValue + trace(
             PeriodicBC::CovShiftForward(U[mu],mu,PeriodicBC::CovShiftForward(U[mu],mu,U[nu]))* // ->->|
         adj(PeriodicBC::CovShiftForward(U[nu],nu,PeriodicBC::CovShiftForward(U[mu],mu,U[mu]))) );
      RectPlaqValue = RectPlaqValue + trace(
             PeriodicBC::CovShiftForward(U[mu],mu,PeriodicBC::CovShiftForward(U[nu],nu,U[nu]))* // ->||
         adj(PeriodicBC::CovShiftForward(U[nu],nu,PeriodicBC::CovShiftForward(U[nu],nu,U[mu]))) );
    }
  }
}

void RectPlaqDeriv(const std::vector<LatticeColourMatrix> &U, LatticeComplex &RectPlaqValue )
{

}
int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  std::vector<int> latt_size  ({16,16,16,32});
  std::vector<int> clatt_size  ({4,4,4,8});
  int orthodir=3;
  int orthosz =latt_size[orthodir];
    
  GridCartesian     Fine(latt_size,simd_layout,mpi_layout);
  GridCartesian     Coarse(clatt_size,simd_layout,mpi_layout);

  LatticeGaugeField Umu(&Fine);

  std::vector<LatticeColourMatrix> U(4,&Fine);
  
  NerscField header;
  
  std::string file("./ckpoint_lat.4000");
  NerscIO::readConfiguration(Umu,header,file);

  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }

  // Painful ; fix syntactical niceness : to check reader
  LatticeComplex LinkTrace(&Fine);
  LinkTrace=zero;
  for(int mu=0;mu<Nd;mu++){
    LinkTrace = LinkTrace + trace(U[mu]);
  }


  LatticeComplex Plaq(&Fine);
  LatticeComplex cPlaq(&Coarse);
  Plaq = zero;
  for(int mu=1;mu<Nd;mu++){
    for(int nu=0;nu<mu;nu++){
      Plaq = Plaq + trace(PeriodicBC::CovShiftForward(U[mu],mu,U[nu])*adj(PeriodicBC::CovShiftForward(U[nu],nu,U[mu])));
    }
  }

  LatticeComplex RectPlaqValue(&Fine);
  double vol = Fine.gSites();
  Complex PlaqScale(1.0/vol/6.0/3.0);
  Complex RectScale(1.0/vol/12.0/3.0);

  std::cout<<GridLogMessage <<"PlaqScale" << PlaqScale<<std::endl;
  std::cout<<GridLogMessage <<"RectScale" << RectScale<<std::endl;

  RectPlaq(U,RectPlaqValue);

  TComplex TRp = sum(RectPlaqValue);
  Complex rp  = TensorRemove(TRp);
  std::cout<<GridLogMessage << "calculated Rect plaquettes A " <<rp*RectScale<<std::endl;


 // Rect Plaq Calc Deriv

  LatticeComplex RectPlaq_d(&Fine);
  RectPlaq_d = zero;
  LatticeColourMatrix ds_U(&Fine);
  LatticeColourMatrix left_2(&Fine);
  LatticeColourMatrix upper_l(&Fine);
  LatticeColourMatrix upper_staple(&Fine);
  LatticeColourMatrix down_staple(&Fine);
  LatticeColourMatrix tmp(&Fine);

  // 2x1 // Each link has 2*(Nd-1) + 4*(Nd-1) = 6(Nd-1) , 1x2 and 2x1 loops attached.
  //     //
  //     // For producing the rectangle term normalised to number of loops
  //     // there are Vol x Nd.(Nd-1) x 2 / 2 distinct loops total. (mu<nu, mu>nu)
  //     // 
  //     // Expect scale factor to be 
  //     // 
  for(int mu=0;mu<Nd;mu++){

    ds_U=zero; // dS / dUmu

    for(int nu=0;nu<Nd;nu++){

    if ( nu != mu ) {
/*     
           (x)  ---> --->  : U(x,mu)*U(x+mu, mu)
*/
           left_2=  PeriodicBC::CovShiftForward(U[mu],mu,U[mu]);
/*
  upper_l =   <---- <---
                         ^
                         | =>tmp
                         (x+2mu)
	   Unu(x+2mu) Umudag(x+mu+nu) Umudag(x+nu)
*/

          tmp=Cshift(U[nu],mu,2);
          upper_l=  PeriodicBC::CovShiftForward(tmp,nu,adj(left_2)); // i.e. upper_l
/*
 upper_staple=               <---- <--- ^
                             |          |
                             V (x)       (x + 2mu)
*/
	  // Unu(x+2mu) Umudag(x+mu+nu) Umudag(x+nu) Unudag(x)
          upper_staple= upper_l*adj(U[nu]);
/*
 down_staple=             ^            
                          |            |
                      (x) <----- <---- V x + 2mu
*/
          down_staple= adj(left_2*tmp)*U[nu];
/*
     ds_U+=                 <---- <--- ^ 
                            |          |
                     (x-mu) V----->    (x + mu)
*/
          tmp=upper_staple*U[mu];
          ds_U+= Cshift(tmp,mu,-1);
/*
     ds_U+=          (x-mu) ^---->      (x + mu)
                            |          |
                            <-----<--- V 
*/          
          tmp=PeriodicBC::CovShiftBackward(U[mu],nu,down_staple);
          ds_U+=Cshift(tmp,mu,-1); 
/*
     ds_U+=                 <----<---- ^ 
                            |          |
                        (x) V     -----> (x + 2mu)
*/
          tmp=Cshift(U[mu],mu,1);
/*
     ds_U+=             (x) ^      ----> (x + 2mu)
                            |          |
                            <---- <----V 
*/
     ds_U+=tmp*(upper_staple+down_staple);

/*****Part 2********/
/*
                     ^
                     | 
   upper=            ^ 
                     | 
                   (x)
*/    
     LatticeColourMatrix up2= PeriodicBC::CovShiftForward(U[nu],nu,U[nu]);
/*
                  <----^
                       | 
   upper_l=            ^ 
                       | 
                 (x)
*/
    // Unu(x+mu)Unu(x+mu+nu) UmuDag(x+nu+nu) lives at X
    upper_l= PeriodicBC::CovShiftForward(Cshift(up2,mu,1),nu,Cshift(adj(U[mu]),nu,1));
/*
                     |<----^
   upper_staple =    V     | 
                     |     ^ 
                 (x) V     | 
*/    
  
    ds_U+= upper_l*adj(up2);

/*
                       |
                       V 
   downer_l=           |  
               (x)<----V                 
*/    
    upper_l= adj(PeriodicBC::CovShiftForward(U[mu],mu,up2)); //downer_l
/*
                     ^     |
   down_staple  =    |     V 
                     ^     | 
                     |     V
                  (x)<----
          down_staple= upper*upper_l;
*/    
    tmp= upper_l*up2;
    ds_U+= Cshift(tmp,nu,-2);


    //TRp = sum(RectPlaq_d);
    //rp  = TensorRemove(TRp);
    //std::cout << GridLogMessage<< " Rect[" << " " << "] = "<<  TensorRemove(TRp) <<std::endl;
    }}

    RectPlaq_d += trace( U[mu]*ds_U) * 0.25;
  }

  TRp = sum(RectPlaq_d);
  rp  = TensorRemove(TRp);
  std::cout<<GridLogMessage << "calculated Rect plaquettes_d " <<rp*RectScale<<std::endl;


  std::vector<TComplex> Plaq_T(orthosz);
  sliceSum(Plaq,Plaq_T,Nd-1);
  int Nt = Plaq_T.size();


  TComplex Plaq_T_sum; 
  Plaq_T_sum=zero;
  for(int t=0;t<Nt;t++){
    Plaq_T_sum = Plaq_T_sum+Plaq_T[t];
    Complex Pt=TensorRemove(Plaq_T[t]);
    std::cout<<GridLogMessage << "sliced ["<<t<<"]" <<Pt*PlaqScale*Real(Nt) << std::endl;
  }

  {
    Complex Pt = TensorRemove(Plaq_T_sum);
    std::cout<<GridLogMessage << "total " <<Pt*PlaqScale<<std::endl;
  }  

  TComplex Tp = sum(Plaq);
  Complex p  = TensorRemove(Tp);
  std::cout<<GridLogMessage << "calculated plaquettes " <<p*PlaqScale<<std::endl;

  RealD avg_plaq = ColourWilsonLoops::avgPlaquette(Umu);
  std::cout<<GridLogMessage << "NEW : calculated real plaquettes " <<avg_plaq<<std::endl;
  // Staple Plaq  
  RealD   StapScale(1.0/vol/6.0/3.0);

  RealD stap_plaq=0.0;
  LatticeColourMatrix stap(&Fine);
  LatticeComplex stap_tr(&Fine);
  for(int mu=0;mu<Nd;mu++){
    ColourWilsonLoops::Staple(stap,Umu,mu);
    stap_tr = trace(U[mu]*stap);
    TComplex Ts = sum(stap_tr);
    Complex s  = TensorRemove(Ts);
    stap_plaq+=real(s);
  }
  std::cout<<GridLogMessage << "NEW : plaquette via staples"<< stap_plaq*StapScale*0.25<< std::endl;
  Complex LinkTraceScale(1.0/vol/4.0/3.0);
  TComplex Tl = sum(LinkTrace);
  Complex l  = TensorRemove(Tl);
  std::cout<<GridLogMessage << "calculated link trace " <<l*LinkTraceScale<<std::endl;

  blockSum(cPlaq,Plaq);
  TComplex TcP = sum(cPlaq);
  Complex ll= TensorRemove(TcP);
  std::cout<<GridLogMessage << "coarsened plaquettes sum to " <<ll*PlaqScale<<std::endl;


  Grid_finalize();
  }


