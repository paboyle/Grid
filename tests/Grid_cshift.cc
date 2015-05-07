#include <Grid.h>
#include <parallelIO/GridNerscIO.h>

using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> simd_layout({1,1,2,2});
  std::vector<int> mpi_layout ({2,2,1,4});
  std::vector<int> latt_size  ({8,8,8,16});
    
  GridCartesian     Fine(latt_size,simd_layout,mpi_layout);
  GridParallelRNG           FineRNG(&Fine);
  FineRNG.SeedRandomDevice();

  LatticeComplex U(&Fine);
  LatticeComplex ShiftU(&Fine);

  LatticeComplex lex(&Fine);
  lex=zero;
  Integer stride =1;
  {
    double nrm;
    LatticeComplex coor(&Fine);

    for(int d=0;d<4;d++){
      LatticeCoordinate(coor,d);
      lex = lex + coor*stride;
      stride=stride*latt_size[d];
    }
    U=lex;
  }

  TComplex cm;
  
  for(int dir=0;dir<4;dir++){
    for(int shift=0;shift<latt_size[dir];shift++){
      if ( Fine.IsBoss() ) 
	std::cout<<"Shifting by "<<shift<<" in direction"<<dir<<std::endl;

	ShiftU  = Cshift(U,dir,shift);    // Shift everything

	std::vector<int> coor(4);

	for(coor[3]=0;coor[3]<latt_size[3];coor[3]++){
	for(coor[2]=0;coor[2]<latt_size[2];coor[2]++){
	for(coor[1]=0;coor[1]<latt_size[1];coor[1]++){
	for(coor[0]=0;coor[0]<latt_size[0];coor[0]++){
	  
	  peekSite(cm,ShiftU,coor);

	  double nrm=norm2(U);

	  std::vector<int> scoor(coor);
	  scoor[dir] = (scoor[dir]+shift)%latt_size[dir];
	  
	  Integer slex = scoor[0]
	    + latt_size[0]*scoor[1]
	    + latt_size[0]*latt_size[1]*scoor[2]
	    + latt_size[0]*latt_size[1]*latt_size[2]*scoor[3];

	  Complex scm(slex);
	  
	  nrm = abs(scm-cm()()());
	  std::vector<int> peer(4);
	  int index=real(cm);
	  Fine.CoorFromIndex(peer,index,latt_size);

	  if (nrm > 0){
	    std::cerr<<"FAIL shift "<< shift<<" in dir "<< dir<<" ["<<coor[0]<<","<<coor[1]<<","<<coor[2]<<","<<coor[3]<<"] = "<< cm()()()<<" expect "<<scm<<"  "<<nrm<<std::endl;
	    std::cerr<<"Got    "<<index<<" " << peer[0]<<","<<peer[1]<<","<<peer[2]<<","<<peer[3]<<std::endl;
	    index=real(scm);
	    Fine.CoorFromIndex(peer,index,latt_size);
	    std::cerr<<"Expect "<<index<<" " << peer[0]<<","<<peer[1]<<","<<peer[2]<<","<<peer[3]<<std::endl;
	  }
	}}}}
    }
  }

  Grid_finalize();
}
