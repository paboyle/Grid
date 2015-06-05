#ifndef GRID_LEBESGUE_H
#define GRID_LEBESGUE_H

#include<vector>

// Lebesgue, Morton, Z-graph ordering assistance
namespace Grid {
  
  class LebesgueOrder { 
  public:

    static int UseLebesgueOrder;

    typedef uint32_t IndexInteger;

    inline IndexInteger Reorder(IndexInteger ss) { 
      return UseLebesgueOrder ? _LebesgueReorder[ss] : ss; 
    };

    IndexInteger alignup(IndexInteger n);

    LebesgueOrder(GridBase *grid);

  private:
    std::vector<IndexInteger> _LebesgueReorder;

  };    
}
#endif
