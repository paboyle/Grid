#include <Grid.h>

using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  int Nd = latt_size.size();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();

  std::vector<int> mask(Nd,1);
  mask[0]=0;

  GridCartesian         Fine  (latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian RBFine(latt_size,simd_layout,mpi_layout,mask,1);

  GridParallelRNG      FineRNG(&Fine);  FineRNG.SeedRandomDevice();

  LatticeComplex U(&Fine);
  LatticeComplex ShiftU(&Fine);
  LatticeComplex rbShiftU(&Fine);
  LatticeComplex Ue(&RBFine); 
  LatticeComplex Uo(&RBFine);
  LatticeComplex ShiftUe(&RBFine);
  LatticeComplex ShiftUo(&RBFine);
  LatticeComplex lex(&Fine);
  lex=zero;
  Integer stride =1;
  {
    double nrm;
    LatticeComplex coor(&Fine);

    for(int d=0;d<Nd;d++){
      //      Integer i=10000;
      Integer i=0;
      LatticeCoordinate(coor,d);
      lex = lex + coor*stride+i;
      stride=stride*latt_size[d];
    }
    U=lex;
  }

  pickCheckerboard(Even,Ue,U);
  pickCheckerboard(Odd,Uo,U);

  //  std::cout << U<<std::endl;
  std::cout << "Ue " <<norm2(Ue)<<std::endl;
  std::cout << "Uo " <<norm2(Uo)<<std::endl;


  TComplex cm;
  for(int dir=0;dir<Nd;dir++){
    if ( dir!=1 ) continue;
    for(int shift=0;shift<latt_size[dir];shift++){

	std::cout<<"Shifting by "<<shift<<" in direction"<<dir<<std::endl;

	//	std::cout<<"Even grid"<<std::endl;
	ShiftUe = Cshift(Ue,dir,shift);    // Shift everything cb by cb
	//	std::cout << "\tShiftUe " <<norm2(ShiftUe)<<std::endl;

	//	std::cout<<"Odd grid"<<std::endl;
	ShiftUo = Cshift(Uo,dir,shift);    
	//	std::cout << "\tShiftUo " <<norm2(ShiftUo)<<std::endl;

	//	std::cout<<"Recombined Even/Odd grids"<<std::endl;
	setCheckerboard(rbShiftU,ShiftUe);
	setCheckerboard(rbShiftU,ShiftUo);
	//	std::cout << "\trbShiftU " <<norm2(rbShiftU)<<std::endl;

	//	std::cout<<"Full grid shift"<<std::endl;
	ShiftU  = Cshift(U,dir,shift);    // Shift everything
	//	std::cout << "\tShiftU " <<norm2(rbShiftU)<<std::endl;

	std::vector<int> coor(4);

	std::cout << "Checking the non-checkerboard shift"<<std::endl;
	for(coor[3]=0;coor[3]<latt_size[3];coor[3]++){
	for(coor[2]=0;coor[2]<latt_size[2];coor[2]++){
	for(coor[1]=0;coor[1]<latt_size[1];coor[1]++){
	for(coor[0]=0;coor[0]<latt_size[0];coor[0]++){
	  
	  peekSite(cm,ShiftU,coor);

	  /////////	  double nrm=norm2(U);

	  std::vector<int> scoor(coor);
	  scoor[dir] = (scoor[dir]+shift)%latt_size[dir];
	  
	  Integer slex = scoor[0]
	    + latt_size[0]*scoor[1]
	    + latt_size[0]*latt_size[1]*scoor[2]
	    + latt_size[0]*latt_size[1]*latt_size[2]*scoor[3];

	  Complex scm(slex);
	  
	  double nrm = abs(scm-cm()()());
	  std::vector<int> peer(4);
	  Complex ctmp = cm;
	  Integer index=real(ctmp);
	  Fine.CoorFromIndex(peer,index,latt_size);

	  if (nrm > 0){
	    std::cerr<<"FAIL shift "<< shift<<" in dir "<< dir
		     <<" ["<<coor[0]<<","<<coor[1]<<","<<coor[2]<<","<<coor[3]<<"] = "
		     << cm()()()<<" expect "<<scm<<"  "<<nrm<<std::endl;
	    std::cerr<<"Got    "<<index<<" " << peer[0]<<","<<peer[1]<<","<<peer[2]<<","<<peer[3]<<std::endl;
	    index=real(scm);
	    Fine.CoorFromIndex(peer,index,latt_size);
	    std::cerr<<"Expect "<<index<<" " << peer[0]<<","<<peer[1]<<","<<peer[2]<<","<<peer[3]<<std::endl;
	    exit(-1);
	  }
	}}}}


	std::cout << "Checking the checkerboard shift"<<std::endl;
	for(coor[3]=0;coor[3]<latt_size[3];coor[3]++){
	for(coor[2]=0;coor[2]<latt_size[2];coor[2]++){
	for(coor[1]=0;coor[1]<latt_size[1];coor[1]++){
	for(coor[0]=0;coor[0]<latt_size[0];coor[0]++){
	  
	  peekSite(cm,rbShiftU,coor);

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
	  Complex ctmp=cm;
	  Integer index=real(ctmp);
	  Fine.CoorFromIndex(peer,index,latt_size);

	  if (nrm > 0){
	    std::cerr<<"FAIL shift "<< shift<<" in dir "<< dir
		     <<" ["<<coor[0]<<","<<coor[1]<<","<<coor[2]<<","<<coor[3]<<"] = "
		     << cm()()()<<" expect "<<scm<<"  "<<nrm<<std::endl;
	    std::cerr<<"Got    "<<index<<" " << peer[0]<<","<<peer[1]<<","<<peer[2]<<","<<peer[3]<<std::endl;
	    index=real(scm);
	    Fine.CoorFromIndex(peer,index,latt_size);
	    std::cerr<<"Expect "<<index<<" " << peer[0]<<","<<peer[1]<<","<<peer[2]<<","<<peer[3]<<std::endl;
	    exit(-1);
	  } else if (0) { 
	    std::cout<<"PASS shift "<< shift<<" in dir "<< dir
		     <<" ["<<coor[0]<<","<<coor[1]<<","<<coor[2]<<","<<coor[3]<<"] = "
		     << cm()()()<<" expect "<<scm<<"  "<<nrm<<std::endl;
	  }
	}}}}

    }
  }

  Grid_finalize();
}
