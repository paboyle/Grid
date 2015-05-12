#include "Grid.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);


  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplexF::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();

  double volume = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
    
  GridCartesian Fine(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian rbFine(latt_size,simd_layout,mpi_layout);
  GridParallelRNG       fRNG(&Fine);
  //  fRNG.SeedRandomDevice();
  std::vector<int> seeds({1,2,3,4});
  fRNG.SeedFixedIntegers(seeds);
  
  LatticeColourMatrix Foo(&Fine);
  LatticeColourMatrix Bar(&Fine);
  LatticeColourMatrix Check(&Fine);
  LatticeColourMatrix Diff(&Fine);
  
  random(fRNG,Foo);
  gaussian(fRNG,Bar);


    for(int dir=0;dir<4;dir++){
      for(int disp=0;disp<Fine._fdimensions[dir];disp++){

	std::cout << "Using stencil to shift dim "<<dir<< " by "<<disp<<std::endl;
	// start to test the Cartesian npoint stencil infrastructure
	int npoint=1;
	std::vector<int> directions(npoint,dir);
	std::vector<int> displacements(npoint,disp);

	CartesianStencil myStencil(&Fine,npoint,0,directions,displacements);

	std::vector<int> ocoor(4);
	for(int o=0;o<Fine.oSites();o++){
	  Fine.oCoorFromOindex(ocoor,o);
	  ocoor[dir]=(ocoor[dir]+disp)%Fine._rdimensions[dir];
	}
	
	std::vector<vColourMatrix,alignedAllocator<vColourMatrix> >  comm_buf(myStencil._unified_buffer_size);
	SimpleCompressor<vColourMatrix> compress;
	myStencil.HaloExchange(Foo,comm_buf,compress);

	Bar = Cshift(Foo,dir,disp);

	// Implement a stencil code that should agree with cshift!
	for(int i=0;i<Check._grid->oSites();i++){

	  int offset = myStencil._offsets [0][i];
	  int  local = myStencil._is_local[0][i];
	  int permute_type = myStencil._permute_type[0];
	  int perm =myStencil._permute[0][i];
	  if ( local && perm )
	    permute(Check._odata[i],Foo._odata[offset],permute_type);
	  else if (local)
	    Check._odata[i] = Foo._odata[offset];
	  else 
	    Check._odata[i] = comm_buf[offset];
	}

	Real nrmC = norm2(Check);
	Real nrmB = norm2(Bar);
	Diff = Check-Bar;
	Real nrm  = norm2(Diff);
	std::cout<<"N2diff ="<<nrm<<" "<<nrmC<<" " <<nrmB<<std::endl;

	Real snrmC =0;
	Real snrmB =0;
	Real snrm  =0;

	std::vector<int> coor(4);
	for(coor[3]=0;coor[3]<latt_size[3]/mpi_layout[3];coor[3]++){
	for(coor[2]=0;coor[2]<latt_size[2]/mpi_layout[2];coor[2]++){
	for(coor[1]=0;coor[1]<latt_size[1]/mpi_layout[1];coor[1]++){
	for(coor[0]=0;coor[0]<latt_size[0]/mpi_layout[0];coor[0]++){

	  Complex diff;
	  ColourMatrix check,bar;
	  peekSite(check,Check,coor);
	  peekSite(bar,Bar,coor);

	  for(int r=0;r<3;r++){
	  for(int c=0;c<3;c++){
            diff =check()()(r,c)-bar()()(r,c);
            double nn=real(conj(diff)*diff);
            if ( nn > 0){
	      printf("Coor (%d %d %d %d) \t rc %d%d \t %le (%le,%le) %le\n",
		     coor[0],coor[1],coor[2],coor[3],r,c,
		     nn,
		     real(check()()(r,c)),
		     imag(check()()(r,c)),
		     real(bar()()(r,c))
		     );
	    }
	    snrmC=snrmC+real(conj(check()()(r,c))*check()()(r,c));
	    snrmB=snrmB+real(conj(bar()()(r,c))*bar()()(r,c));
	    snrm=snrm+nn;
	  }}
	 
	}}}}

	std::cout<<"scalar N2diff = "<<snrm<<" " <<snrmC<<" "<<snrmB<<std::endl;


      }
    }

 Grid_finalize();
}
