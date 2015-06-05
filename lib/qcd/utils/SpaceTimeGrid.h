#ifndef GRID_QCD_SPACE_TIME_GRID_H
#define GRID_QCD_SPACE_TIME_GRID_H
namespace Grid {
namespace QCD {

class SpaceTimeGrid {
 public:

  static GridCartesian         *makeFourDimGrid(const std::vector<int> & latt,const std::vector<int> &simd,const std::vector<int> &mpi);
  static GridRedBlackCartesian *makeFourDimRedBlackGrid       (const GridCartesian *FourDimGrid);
  static GridCartesian         *makeFiveDimGrid        (int Ls,const GridCartesian *FourDimGrid);
  static GridRedBlackCartesian *makeFiveDimRedBlackGrid(int Ls,const GridCartesian *FourDimGrid);

};

}}

#endif
