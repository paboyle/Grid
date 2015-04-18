#ifndef GRID_LATTICE_TRANSFER_H
#define GRID_LATTICE_TRANSFER_H

namespace Grid {

  ////////////////////////////////////////////////////////////////////////////////////////////
  // remove and insert a half checkerboard
  ////////////////////////////////////////////////////////////////////////////////////////////
  template<class vobj> inline void pickCheckerboard(int cb,Lattice<vobj> &half,const Lattice<vobj> &full){
    half.checkerboard = cb;
    int ssh=0;
#pragma omp parallel for
    for(int ss=0;ss<full._grid->oSites();ss++){
      std::vector<int> coor;
      int cbos;
      
      full._grid->oCoorFromOindex(coor,ss);
      cbos=half._grid->CheckerBoard(coor);
      
      if (cbos==cb) {
	
	half._odata[ssh] = full._odata[ss];
	ssh++;
      }
    }
  }
  template<class vobj> inline void setCheckerboard(Lattice<vobj> &full,const Lattice<vobj> &half){
    int cb = half.checkerboard;
    int ssh=0;
#pragma omp parallel for
    for(int ss=0;ss<full._grid->oSites();ss++){
      std::vector<int> coor;
      int cbos;
      
      full._grid->oCoorFromOindex(coor,ss);
      cbos=half._grid->CheckerBoard(coor);
      
      if (cbos==cb) {
	full._odata[ss]=half._odata[ssh];
	ssh++;
      }
    }
  }
  
}
#endif
