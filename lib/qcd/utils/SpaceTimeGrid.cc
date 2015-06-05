#include <Grid.h>

namespace Grid { 
  namespace QCD {

/////////////////////////////////////////////////////////////////
// Public interface
/////////////////////////////////////////////////////////////////
GridCartesian *SpaceTimeGrid::makeFourDimGrid(const std::vector<int> & latt,const std::vector<int> &simd,const std::vector<int> &mpi)
{
  return new GridCartesian(latt,simd,mpi); 
}
GridRedBlackCartesian *SpaceTimeGrid::makeFourDimRedBlackGrid(const GridCartesian *FourDimGrid)
{
  return new GridRedBlackCartesian(FourDimGrid); 
}

GridCartesian         *SpaceTimeGrid::makeFiveDimGrid(int Ls,const GridCartesian *FourDimGrid)
{
  int N4=FourDimGrid->_ndimension;

  std::vector<int> latt5(1,Ls);
  std::vector<int> simd5(1,1);
  std::vector<int>  mpi5(1,1);
  
  for(int d=0;d<N4;d++){
    latt5.push_back(FourDimGrid->_fdimensions[d]);
    simd5.push_back(FourDimGrid->_simd_layout[d]);
     mpi5.push_back(FourDimGrid->_processors[d]);
  }
  return new GridCartesian(latt5,simd5,mpi5); 
}

GridRedBlackCartesian *SpaceTimeGrid::makeFiveDimRedBlackGrid(int Ls,const GridCartesian *FourDimGrid)
{
  int N4=FourDimGrid->_ndimension;
  int cbd=1;
  std::vector<int> latt5(1,Ls);
  std::vector<int> simd5(1,1);
  std::vector<int>  mpi5(1,1);
  std::vector<int>   cb5(1,0);
    
  for(int d=0;d<N4;d++){
    latt5.push_back(FourDimGrid->_fdimensions[d]);
    simd5.push_back(FourDimGrid->_simd_layout[d]);
     mpi5.push_back(FourDimGrid->_processors[d]);
      cb5.push_back(  1);
    }
  return new GridRedBlackCartesian(latt5,simd5,mpi5,cb5,cbd); 
}

}}
