#ifndef GRID_LEBESGUE_H
#define GRID_LEBESGUE_H

#include<vector>

// Lebesgue, Morton, Z-graph ordering assistance
namespace Grid {
  
  class LebesgueOrder { 
  public:

    typedef int32_t IndexInteger;
    
    static int UseLebesgueOrder;
    GridBase *grid;

  public:
    LebesgueOrder(GridBase *_grid);

    inline IndexInteger Reorder(IndexInteger ss) { 
      return _LebesgueReorder[ss] ;
    };

    ////////////////////////////
    // Space filling fractal for cache oblivious
    ////////////////////////////
    void ZGraph(void);
    IndexInteger alignup(IndexInteger n);

    /////////////////////////////////
    // Cartesian stencil blocking strategy
    /////////////////////////////////
    static std::vector<int> Block;
    void CartesianBlocking(void);
    void IterateO(int ND,int dim,
		  std::vector<IndexInteger> & xo,
		  std::vector<IndexInteger> & xi,
		  std::vector<IndexInteger> &dims);
    void IterateI(int ND,int dim,
		  std::vector<IndexInteger> & xo,
		  std::vector<IndexInteger> & xi,
		  std::vector<IndexInteger> &dims);

  private:
    std::vector<IndexInteger> _LebesgueReorder;

  };    
}
#endif
