/*
  2f Full det MdagM 10^6 force ~ 1.3e7
rid : Message : 1767.283471 s : +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Grid : Message : 1767.283476 s : S1 : 1.52885e+09
Grid : Message : 1767.283480 s : S2 : 1.52886e+09
Grid : Message : 1767.283482 s : dS : 8877.34
Grid : Message : 1767.283483 s : dSpred : 8877.7
Grid : Message : 1767.283484 s : diff : -0.360484
Grid : Message : 1767.283485 s : *********************************************************

  2f Full det MpcdagMpc 10^6 force ~ 1.8e6
Grid : Message : 2399.576962 s : +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Grid : Message : 2399.576968 s : S1 : 1.52885e+09
Grid : Message : 2399.576972 s : S2 : 1.52886e+09
Grid : Message : 2399.576974 s : dS : 9728.49
Grid : Message : 2399.576975 s : dSpred : 9726.58
Grid : Message : 2399.576976 s : diff : 1.90683
Grid : Message : 2399.576977 s : *********************************************************

  2f bdy MdagM 1500 force Force ~ 2800
Grid : Message : 4622.385061 s : +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Grid : Message : 4622.385067 s : S1 : 1.52885e+09
Grid : Message : 4622.385071 s : S2 : 1.52885e+09
Grid : Message : 4622.385072 s : dS : 25.4944
Grid : Message : 4622.385073 s : dSpred : 25.4672
Grid : Message : 4622.385074 s : diff : 0.0271414
Grid : Message : 4622.385075 s : *********************************************************

  2f bdy MpcdagMpc 10^6 force   Force ~ 2200
Grid : Message : 4622.385061 s : +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Grid : Message : 4622.385067 s : S1 : 1.52885e+09
Grid : Message : 4622.385071 s : S2 : 1.52885e+09
Grid : Message : 4622.385072 s : dS : 25.4944
Grid : Message : 4622.385073 s : dSpred : 25.4672
Grid : Message : 4622.385074 s : diff : 0.0271414
Grid : Message : 4622.385075 s : *********************************************************
  
  1f Bdy Det
//
// These all had tol set by OFRp,  not through MDpoles
// So assumptions it was Remez might be wrong.
//
Optimisation log:  looser rational AND MD tolerances sloppy
MobiusForce.221179 -- same as HMC. dS is mispredicted Forece  ~2.8
Grid : Message : 6582.258991 s : dS : 0.024478
Grid : Message : 6582.258992 s : dSpred : 0.00791876
Grid : Message : 6582.258994 s : diff : 0.0165592

MobiusForce.221193 -- tight rational AND MD tolerances to 1e-8 ~ 2.8 same
Grid : Message : 1964.939209 s : S1 : 7.64404e+08
Grid : Message : 1964.939213 s : S2 : 7.64404e+08
Grid : Message : 1964.939215 s : dS : -0.00775838 <--- too loose even on action
Grid : Message : 1964.939216 s : dSpred : -0.00416793 
Grid : Message : 1964.939217 s : diff : -0.00359045

MobiusForce.221394 -- tight rational, MD tol sloppy Force ~ 2.8
Grid : Message : 2376.921950 s : S1 : 764404436.44069
Grid : Message : 2376.921954 s : S2 : 764404436.43299
Grid : Message : 2376.921956 s : dS : -0.0076971054077148
Grid : Message : 2376.921958 s : dSpred : -0.0041610472282526
Grid : Message : 2376.921959 s : diff : -0.0035360581794623


MobiusForce.221587 -- slightly sloppier action, coming from tol array
                   -- much sloppier force
		   -- degree 18
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-6,3.0e-7,1.0e-7,1.0e-7,  // Orig sloppy
Grid : Message : 2438.875507 s : S1 : 764404436.42251
Grid : Message : 2438.875512 s : S2 : 764404436.4148
Grid : Message : 2438.875514 s : dS : -0.0077102184295654
Grid : Message : 2438.875516 s : dSpred : -0.0075684496959103
Grid : Message : 2438.875517 s : diff : -0.00014176873365508

MobiusForce.221639       3.0e-6,1.0e-6,1.0e-7,1.0e-7, // soften convergence more

Grid : Message : 2373.927550 s : +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Grid : Message : 2373.927600 s : S1 : 764404436.42251
Grid : Message : 2373.927640 s : S2 : 764404436.4148
Grid : Message : 2373.927660 s : dS : -0.0077102184295654
Grid : Message : 2373.927680 s : dSpred : -0.0075993463919849
Grid : Message : 2373.927690 s : diff : -0.00011087203758051
Grid : Message : 2373.927700 s : *********************************************************


Grid : Message : 69.269319 s :  ApproxPowerMD shift[0]  pole    9.5166866092503e-06 residue -2.0047722631555e-08 tol     3e-06
Grid : Message : 69.269321 s :  ApproxPowerMD shift[1]  pole    4.7123486192778e-05 residue -1.316766030683e-07 tol     1e-06
Grid : Message : 69.269323 s :  ApproxPowerMD shift[2]  pole    0.00014860967743736 residue -6.109883117444e-07 tol     1e-07
Grid : Message : 69.269325 s :  ApproxPowerMD shift[3]  pole    0.00041055696132763 residue -2.6088717433891e-06 tol     1e-07
Grid : Message : 69.269327 s :  ApproxPowerMD shift[4]  pole    0.0010822555692906 residue -1.0853799412802e-05 tol     1e-08
Grid : Message : 69.269329 s :  ApproxPowerMD shift[5]  pole    0.0028029613512087 residue -4.4741734470158e-05 tol     1e-08
Grid : Message : 69.269331 s :  ApproxPowerMD shift[6]  pole    0.0072103567378527 residue -0.00018380499193253 tol     1e-08

rusher 96I]$ more MobiusForce.221887
      1.0e-5,3.0e-6,3.0e-7,1.0e-7, // soften convergence more more
// <-- this is the dirichlet solve, why poorer conditioned???
Grid : Message : 1627.226206 s : ConjugateGradientMultiShift k=3643 Shift 3 has converged
Grid : Message : 1667.373045 s : ConjugateGradientMultiShift k=5381 Shift 2 has converged
Grid : Message : 1705.236992 s : ConjugateGradientMultiShift k=7063 Shift 1 has converged
Grid : Message : 1752.493182 s : ConjugateGradientMultiShift k=9220 Shift 0 has converged
// 
//Grid : Message : 1414.837250 s : OneFlavourEvenOddRatioRationalPseudoFermionAction deriv: doing (M^dag M)^{-1/2} ( (V^dag V)^{1/4} Phi)
Grid : Message : 1523.416680 s : ConjugateGradientMultiShift k=3846 Shift 2 has converged
Grid : Message : 1530.798503 s : ConjugateGradientMultiShift k=4143 Shift 1 has converged
Grid : Message : 1536.153421 s : ConjugateGradientMultiShift k=4353 Shift 0 has converged <-- this is the non-dirichlet solve

Grid : Message : 2339.927565 s : +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Grid : Message : 2339.927571 s : S1 : 764404436.42251
Grid : Message : 2339.927575 s : S2 : 764404436.4148
Grid : Message : 2339.927577 s : dS : -0.0077102184295654
Grid : Message : 2339.927579 s : dSpred : -0.0068752425267964
Grid : Message : 2339.927580 s : diff : -0.00083497590276901
Grid : Message : 2339.927581 s : *********************************************************
Grid : Message : 2339.927582 s : Done
Grid : Message : 2339.927582 s : *********************************************************

Force 76 S {S {S {(9.0175185326468,-3.5764415623768e-36)}}}
Force 77 S {S {S {(4.1289977678493,-4.3364721285803e-37)}}}
Force 78 S {S {S {(3.2299269465841,6.0391022273495e-37)}}}
Force 79 S {S {S {(3.0051199649288,-9.6243599973575e-37)}}}
Force 80 S {S {S {(2.8924316727872,-1.3371248240604e-37)}}}
Force 81 S {S {S {(2.8270868791781,1.792628885004e-37)}}}
Force 82 S {S {S {(2.8676819960087,-1.3518185034456e-36)}}}
Force 83 S {S {S {(2.7724152154523,1.4950818774521e-37)}}}
Force 84 S {S {S {(3.0204624534964,-9.6475025423893e-36)}}}
Force 85 S {S {S {(2.8631304063459,2.2426228161781e-37)}}}
Force 86 S {S {S {(2.9025673908905,-1.3942465026706e-36)}}}
Force 87 S {S {S {(2.8553405232646,-2.0938493124022e-38)}}}
Force 88 S {S {S {(3.2820184381375,-1.422348164495e-36)}}}
Force 89 S {S {S {(3.8974980085791,1.1682209795266e-35)}}}
Force 90 S {S {S {(4.660053618223,-1.4399805797573e-37)}}}
Force 91 S {S {S {(6.7993872372366,1.4524702072348e-36)}}}
Full
Grid : Message : 1523.416680 s : ConjugateGradientMultiShift k=3846 Shift 2 has converged
Grid : Message : 1530.798503 s : ConjugateGradientMultiShift k=4143 Shift 1 has converged
Grid : Message : 1536.153421 s : ConjugateGradientMultiShift k=4353 Shift 0 has converged
PV solve depth 3
Grid : Message : 1667.373045 s : ConjugateGradientMultiShift k=5381 Shift 2 has converged
Grid : Message : 1705.236992 s : ConjugateGradientMultiShift k=7063 Shift 1 has converged
Grid : Message : 1752.493182 s : ConjugateGradientMultiShift k=9220 Shift 0 has converged

MobiusForce.222490 depth 1
Grid : Message : 2155.595070 s : +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Grid : Message : 2155.595076 s : S1 : 764404436.37475
Grid : Message : 2155.595080 s : S2 : 764404436.21131
Grid : Message : 2155.595082 s : dS : -0.16344606876373
Grid : Message : 2155.595084 s : dSpred : -0.16235663327375
Grid : Message : 2155.595085 s : diff : -0.0010894354899788

Force 4 S {S {S {(24.512489110423,-7.4203080895657e-36)}}}
Force 5 S {S {S {(14.442663101577,7.3909207307951e-37)}}}
Force 6 S {S {S {(12.298567945213,2.1989091200069e-36)}}}
Force 7 S {S {S {(11.582362859271,-2.2540104177017e-36)}}}
Force 8 S {S {S {(11.465725500906,-2.9512255045332e-36)}}}
Force 9 S {S {S {(10.869067954412,-2.8388188572358e-36)}}}
Force 10 S {S {S {(10.937111429576,-3.3530976357206e-36)}}}
Force 11 S {S {S {(11.23500117508,-1.4487967873885e-36)}}}
Force 12 S {S {S {(10.900736551834,5.1427877848475e-36)}}} Force is bigger
Force 13 S {S {S {(10.951921323651,-1.2098775605838e-35)}}}
Force 14 S {S {S {(10.676529230575,-2.50527233519e-36)}}}
Force 15 S {S {S {(10.98568474467,3.2193851533145e-36)}}}
Force 16 S {S {S {(11.931707726568,-8.5223340434616e-37)}}}
Force 17 S {S {S {(13.751904678482,7.6337337826369e-36)}}}
Force 18 S {S {S {(17.518955473833,1.8073225643893e-36)}}}
Force 19 S {S {S {(20.36519304598,-2.5184966466368e-36)}}}
Full solve
Grid : Message : 1441.297575 s : ConjugateGradientMultiShift k=3846 Shift 2 has converged
Grid : Message : 1449.206520 s : ConjugateGradientMultiShift k=4143 Shift 1 has converged
Grid : Message : 1454.352909 s : ConjugateGradientMultiShift k=4353 Shift 0 has converged

Dirichlet solve -- why so expensive??
Spectral radius worse?
Grid : Message : 1571.887003 s : ConjugateGradientMultiShift k=5195 Shift 2 has converged
Grid : Message : 1599.543760 s : ConjugateGradientMultiShift k=6508 Shift 1 has converged
Grid : Message : 1625.368198 s : ConjugateGradientMultiShift k=7819 Shift 0 has converged


dS is much bigger.


MobiusForce.223606
Grid : Message : 1123.276405 s : ConjugateGradientMultiShift k=3273 Shift 0 has converged
Grid : Message : 1125.945359 s : ConjugateGradientMultiShift k=3407 Shift 1 has converged
Grid : Message : 1127.896580 s : ConjugateGradientMultiShift k=3508 Shift 2 has converged <-- 2 takes longer
first (bdy) hasenbusch mass raised to 0.005 -- reduces Dirchlet solve cost
Force looks ok still
Grid : Message : 1510.884960 s : OneFlavourEvenOddRatioRationalPseudoFermionAction compute action: complete
Grid : Message : 1510.969380 s : +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Grid : Message : 1510.969440 s : S1 : 764404436.37475
Grid : Message : 1510.969480 s : S2 : 764404436.17379
Grid : Message : 1510.969500 s : dS : -0.20095825195312
Grid : Message : 1510.969520 s : dSpred : -0.20025674631954
Grid : Message : 1510.969530 s : diff : -0.00070150563358654
Force 76 S {S {S {(24.161229317675,2.0147973173094e-35)}}}
Force 77 S {S {S {(15.841085162729,3.983456481349e-36)}}}
Force 78 S {S {S {(11.031761776856,9.0394046210295e-35)}}}
Force 79 S {S {S {(12.177830066719,1.583978637733e-36)}}}
Force 80 S {S {S {(9.8372072482222,6.4284847310594e-37)}}}
Force 81 S {S {S {(9.6588863493149,1.0501572656659e-35)}}}
Force 82 S {S {S {(10.623076227724,-4.4161853392455e-35)}}}
Force 83 S {S {S {(8.9477003784221,-7.067659784319e-37)}}}
Force 84 S {S {S {(9.7663166497594,-2.1014900256825e-35)}}}
Force 85 S {S {S {(8.9992648919057,-4.7107936109203e-36)}}}
Force 86 S {S {S {(9.0399987268337,6.4652189295226e-37)}}}
Force 87 S {S {S {(9.1319052497073,7.9566273871284e-37)}}}
Force 88 S {S {S {(10.094569606113,-1.263656427134e-37)}}}
Force 89 S {S {S {(11.563679905523,-1.2777623593438e-35)}}}
Force 90 S {S {S {(13.653150474463,2.9093485182852e-37)}}}
Force 91 S {S {S {(16.303719912019,2.9857556510886e-36)}}}

MobiusForce.223749
first (bdy) hasenbusch mass raised to 0.01 -- reduces Dirchlet solve cost
Grid : Message : 1374.472462 s : OneFlavourEvenOddRatioRationalPseudoFermionAction compute action: complete
Grid : Message : 1374.479206 s : +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Grid : Message : 1374.479211 s : S1 : 764404436.37428
Grid : Message : 1374.479215 s : S2 : 764404436.20009
Grid : Message : 1374.479217 s : dS : -0.17418932914734
Grid : Message : 1374.479219 s : dSpred : -0.17358090105485
Grid : Message : 1374.479220 s : diff : -0.00060842809248995
Force 76 S {S {S {(27.006858541753,4.2141472476979e-36)}}}
Force 77 S {S {S {(19.388701462694,-5.1620365048422e-35)}}}
Force 78 S {S {S {(13.502424539662,-2.4038859474316e-35)}}}
Force 79 S {S {S {(15.555776987064,6.0567346426118e-36)}}}
Force 80 S {S {S {(12.752116522904,-2.3720006631655e-35)}}}
Force 81 S {S {S {(12.656857824233,1.6912424972456e-35)}}}
Force 82 S {S {S {(15.159284452724,5.0898905390605e-36)}}}
Force 83 S {S {S {(12.222695136014,-2.2061824913027e-35)}}}
Force 84 S {S {S {(12.92077598466,9.6287681011731e-36)}}}
Force 85 S {S {S {(11.884630495484,2.822655809912e-36)}}}
Force 86 S {S {S {(11.896353116174,1.0926219990893e-35)}}}
Force 87 S {S {S {(11.557019282287,2.1532117771187e-35)}}}
Force 88 S {S {S {(11.945108384613,-3.0210204816133e-36)}}}
Force 89 S {S {S {(13.295373801078,7.3115748621146e-36)}}}
Force 90 S {S {S {(15.373728471417,-7.4923071185536e-36)}}}
Force 91 S {S {S {(17.348173714234,1.0344350287236e-36)}}}

MobiusForce.223829
      1.0e-5,5.0e-6,1.0e-6,1.0e-7, // soften convergence more more
Grid : Message : 1000.951387 s : ConjugateGradientMultiShift k=1881 Shift 0 has converged
Grid : Message : 1002.619542 s : ConjugateGradientMultiShift k=1960 Shift 1 has converged
Grid : Message : 1003.726982 s : ConjugateGradientMultiShift k=2014 Shift 4 has converged
Grid : Message : 1005.698741 s : ConjugateGradientMultiShift k=2113 Shift 2 has converged
Grid : Message : 1007.320875 s : ConjugateGradientMultiShift k=2197 Shift 3 has converged
Grid : Message : 1351.171259 s : S1 : 764404436.37428
Grid : Message : 1351.171263 s : S2 : 764404436.20009
Grid : Message : 1351.171265 s : dS : -0.17418932914734
Grid : Message : 1351.171266 s : dSpred : -0.1743248065338
Grid : Message : 1351.171267 s : diff : 0.00013547738646566
Force 76 S {S {S {(27.004288088317,6.035575744297e-35)}}}
Force 77 S {S {S {(19.388023720604,-6.9736202362532e-36)}}}
Force 78 S {S {S {(13.502663916173,6.4067380855692e-35)}}}
Force 79 S {S {S {(15.55135748152,1.7219522871608e-35)}}}
Force 80 S {S {S {(12.75135802213,-1.1303847551095e-35)}}}
Force 81 S {S {S {(12.655732786276,1.689773129307e-36)}}}
Force 82 S {S {S {(15.158469055699,-6.7205950772387e-35)}}}
Force 83 S {S {S {(12.222907191126,-1.6775773754173e-35)}}}
Force 84 S {S {S {(12.916025368247,-1.9641041234302e-35)}}}
Force 85 S {S {S {(11.881879452577,-2.3054382955502e-36)}}}
Force 86 S {S {S {(11.897253557199,-3.3617669065579e-35)}}}
Force 87 S {S {S {(11.55717723524,-1.8690360178074e-36)}}}
Force 88 S {S {S {(11.945590605851,-6.7208889508264e-36)}}}
Force 89 S {S {S {(13.298173932749,-1.0322309768158e-35)}}}
Force 90 S {S {S {(15.373845416836,7.4158999857501e-36)}}}
Force 91 S {S {S {(17.348058307158,-1.8514036025451e-36)}}}
-- could make the stopping condition mandatory if shift 0 is converged.
-- Save 20% of iterations and single tunable
*/

//
/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_double_ratio.cc

    Copyright (C) 2022

Author: Peter Boyle <pboyle@bnl.gov>

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

using namespace std;
using namespace Grid;

typedef MobiusFermionD FermionAction;
typedef WilsonImplD FimplD;
typedef WilsonImplD FermionImplPolicy;

template<class Gimpl>
void ForceTest(Action<LatticeGaugeField> &action,LatticeGaugeField & U,MomentumFilterBase<LatticeGaugeField> &Filter)
{
  GridBase *UGrid = U.Grid();

  std::vector<int> seeds({1,2,3,5});
  GridSerialRNG            sRNG;         sRNG.SeedFixedIntegers(seeds);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds);

  LatticeColourMatrix Pmu(UGrid); 
  LatticeGaugeField P(UGrid); 
  LatticeGaugeField UdSdU(UGrid); 

  std::cout << GridLogMessage << "*********************************************************"<<std::endl;
  std::cout << GridLogMessage << " Force test for "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "*********************************************************"<<std::endl;
  
  RealD eps=0.005;

  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout << GridLogMessage << " Refresh "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  
  Gimpl::generate_momenta(P,sRNG,RNG4);
  Filter.applyFilter(P);

  FieldMetaData header;
  std::string file("./ckpoint_lat.2000");
  NerscIO::readConfiguration(U,header,file);
 
  action.refresh(U,sRNG,RNG4);

  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout << GridLogMessage << " Action "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;

  RealD S1 = action.S(U);

  Gimpl::update_field(P,U,eps);

  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout << GridLogMessage << " Derivative "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  action.deriv(U,UdSdU);
  UdSdU = Ta(UdSdU);
  Filter.applyFilter(UdSdU);

  DumpSliceNorm("Force",UdSdU,Nd-1);
  
  Gimpl::update_field(P,U,eps);
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout << GridLogMessage << " Action "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  
  RealD S2 = action.S(U);

  // Use the derivative
  LatticeComplex dS(UGrid); dS = Zero();
  for(int mu=0;mu<Nd;mu++){
    auto UdSdUmu = PeekIndex<LorentzIndex>(UdSdU,mu);
    Pmu= PeekIndex<LorentzIndex>(P,mu);
    dS = dS - trace(Pmu*UdSdUmu)*eps*2.0*2.0;
  }
  ComplexD dSpred    = sum(dS);
  RealD diff =  S2-S1-dSpred.real();

  std::cout<< GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout<< GridLogMessage << "S1 : "<< S1    <<std::endl;
  std::cout<< GridLogMessage << "S2 : "<< S2    <<std::endl;
  std::cout<< GridLogMessage << "dS : "<< S2-S1 <<std::endl;
  std::cout<< GridLogMessage << "dSpred : "<< dSpred.real() <<std::endl;
  std::cout<< GridLogMessage << "diff : "<< diff<<std::endl;
  std::cout<< GridLogMessage << "*********************************************************"<<std::endl;
  //  assert(diff<1.0);
  std::cout<< GridLogMessage << "Done" <<std::endl;
  std::cout << GridLogMessage << "*********************************************************"<<std::endl;
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::cout << std::setprecision(14);
  Coordinate latt_size   = GridDefaultLatt();
  Coordinate mpi_layout  = GridDefaultMpi();
  Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  Coordinate shm;
  GlobalSharedMemory::GetShmDims(mpi_layout,shm);

  const int Ls=12;
  const int Nt = latt_size[3];
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);
 
  ////////////////////////////////////////////////////////////////
  // Domain decomposed operator
  ////////////////////////////////////////////////////////////////
  Coordinate CommDim(Nd);
  for(int d=0;d<Nd;d++) CommDim[d]= (mpi_layout[d]/shm[d])>1 ? 1 : 0;

  Coordinate NonDirichlet(Nd+1,0);
  Coordinate Dirichlet(Nd+1,0);
  Dirichlet[1] = CommDim[0]*latt_size[0]/mpi_layout[0] * shm[0];
  Dirichlet[2] = CommDim[1]*latt_size[1]/mpi_layout[1] * shm[1];
  Dirichlet[3] = CommDim[2]*latt_size[2]/mpi_layout[2] * shm[2];
  Dirichlet[4] = CommDim[3]*latt_size[3]/mpi_layout[3] * shm[3];

  Coordinate Block4(Nd);
  Block4[0] = Dirichlet[1];
  Block4[1] = Dirichlet[2];
  Block4[2] = Dirichlet[3];
  Block4[3] = Dirichlet[4];

  std::vector<Complex> boundary = {1,1,1,-1};
  FermionAction::ImplParams Params(boundary);
  FermionAction::ImplParams ParamsDir(boundary);
  Params.dirichlet=NonDirichlet;
  ParamsDir.dirichlet=Dirichlet;
  ParamsDir.partialDirichlet=1;

  ///////////////////// Gauge Field and Gauge Forces ////////////////////////////
  LatticeGaugeField U(UGrid);

  RealD beta=6.0;
  WilsonGaugeActionR PlaqAction(beta);
  IwasakiGaugeActionR RectAction(beta);

  MomentumFilterNone<LatticeGaugeField> FilterNone;
  ForceTest<GimplTypesR>(PlaqAction,U,FilterNone);
  ForceTest<GimplTypesR>(RectAction,U,FilterNone);

  ////////////////////////////////////
  // Action
  ////////////////////////////////////
  RealD mass=0.00078; 
  RealD dmass=0.01;
  RealD pvmass=1.0; 
  RealD M5=1.8; 
  RealD b=1.5;
  RealD c=0.5;
  
  // Double versions
  FermionAction DdwfPeriodic(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c,Params);
  FermionAction PVPeriodic  (U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,pvmass,M5,b,c,Params);
  FermionAction DdwfDirichlet(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,dmass,M5,b,c,ParamsDir);

  double StoppingCondition = 1.0e-8;
  double MaxCGIterations = 50000;
  ConjugateGradient<LatticeFermion>  CG(StoppingCondition,MaxCGIterations);
  
  //////////////////// Two Flavour Determinant Ratio ///////////////////////////////
  TwoFlavourRatioPseudoFermionAction<FimplD> Nf2(PVPeriodic, DdwfPeriodic,CG,CG);
  //  ForceTest<GimplTypesR>(Nf2,U,FilterNone);

  //////////////////// Two Flavour Determinant force test Even Odd ///////////////////////////////
  TwoFlavourEvenOddRatioPseudoFermionAction<FimplD> Nf2eo(PVPeriodic, DdwfPeriodic,CG,CG);
  //  ForceTest<GimplTypesR>(Nf2eo,U,FilterNone);

  //////////////////// Domain forces ////////////////////
  int Width=4;
  DDHMCFilter<WilsonImplD::Field> DDHMCFilter(Block4,Width);
  
  //////////////////// Two flavour boundary det  ////////////////////
  TwoFlavourRatioPseudoFermionAction<FimplD> BdyNf2(DdwfDirichlet, DdwfPeriodic,CG,CG);
  //  ForceTest<GimplTypesR>(BdyNf2,U,DDHMCFilter);

  //////////////////// Two flavour eo boundary det  ////////////////////
  TwoFlavourEvenOddRatioPseudoFermionAction<FimplD> BdyNf2eo(DdwfDirichlet, DdwfPeriodic,CG,CG);
  //  ForceTest<GimplTypesR>(BdyNf2eo,U,DDHMCFilter);

  //////////////////// One flavour boundary det  ////////////////////
  RationalActionParams OFRp; // Up/down
  OFRp.lo       = 6.0e-5;
  OFRp.hi       = 90.0;
  OFRp.inv_pow  = 2;
  OFRp.MaxIter  = SP_iters; // get most shifts by 2000, stop sharing space
  OFRp.action_tolerance= 1.0e-8;
  OFRp.action_degree   = 18;
  OFRp.md_tolerance= 1.0e-5;
  OFRp.md_degree   = 14;
  //  OFRp.degree   = 20; converges
  //  OFRp.degree   = 16;
  OFRp.precision= 80;
  OFRp.BoundsCheckFreq=0;
  /*
  OneFlavourRationalParams OFRp; // Up/down
  OFRp.lo       = 4.0e-5;
  OFRp.hi       = 90.0;
  OFRp.MaxIter  = 60000;
  OFRp.tolerance= 1.0e-9;
  OFRp.mdtolerance= 1.0e-8;
  OFRp.degree   = 18;
  OFRp.precision= 80;
  OFRp.BoundsCheckFreq=0;
  */
  std::vector<RealD> ActionTolByPole({
      1.0e-7,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8
    });
  std::vector<RealD> MDTolByPole({
      1.6e-5,5.0e-6,1.0e-6,3.0e-7, // soften convergence more more
      //      1.0e-6,3.0e-7,1.0e-7,1.0e-7,
      //      3.0e-6,1.0e-6,1.0e-7,1.0e-7, // soften convergence
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8
    });
  /*
  std::vector<RealD> ActionTolByPole({
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8
    });
  std::vector<RealD> MDTolByPole({
      1.0e-5,5.0e-6,1.0e-6,1.0e-7, // soften convergence more more
      //      3.0e-6,1.0e-6,1.0e-7,1.0e-7, // soften convergence more
      //      1.0e-6,3.0e-7,1.0e-7,1.0e-7,  // Orig sloppy
      //      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8
    });
  */
  OneFlavourEvenOddRatioRationalPseudoFermionAction<FermionImplPolicy> BdySqrt(DdwfDirichlet,DdwfPeriodic,OFRp);
  BdySqrt.SetTolerances(ActionTolByPole,MDTolByPole);
  ForceTest<GimplTypesR>(BdySqrt,U,DDHMCFilter);

  Grid_finalize();
}
