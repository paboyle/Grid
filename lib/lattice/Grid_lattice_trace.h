#ifndef GRID_LATTICE_TRACE_H
#define GRID_LATTICE_TRACE_H

///////////////////////////////////////////////
// Tracing, transposing, peeking, poking
///////////////////////////////////////////////

namespace Grid {

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Trace
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    template<class vobj>
    inline auto trace(const Lattice<vobj> &lhs)
      -> Lattice<decltype(trace(lhs._odata[0]))>
    {
      Lattice<decltype(trace(lhs._odata[0]))> ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = trace(lhs._odata[ss]);
        }
        return ret;
    };
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Trace Index level dependent operation
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    template<int Index,class vobj>
    inline auto traceIndex(const Lattice<vobj> &lhs)
      -> Lattice<decltype(traceIndex<Index>(lhs._odata[0]))>
    {
      Lattice<decltype(traceIndex<Index>(lhs._odata[0]))> ret(lhs._grid);
#pragma omp parallel for
      for(int ss=0;ss<lhs._grid->oSites();ss++){
	ret._odata[ss] = traceIndex<Index>(lhs._odata[ss]);
      }
      return ret;
    };


}
#endif

