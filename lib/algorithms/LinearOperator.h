#ifndef  GRID_ALGORITHM_LINEAR_OP_H
#define  GRID_ALGORITHM_LINEAR_OP_H

#include <Grid.h>

namespace Grid {

  // Red black cases?
    template<class vtype> class LinearOperatorBase {
    public:
      void M(const Lattice<vtype> &in, Lattice<vtype> &out){ assert(0);}
      void Mdag(const Lattice<vtype> &in, Lattice<vtype> &out){ assert(0);}
      void MdagM(const Lattice<vtype> &in, Lattice<vtype> &out){ assert(0);}
    };

}

#endif
