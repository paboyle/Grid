#ifndef GRID_LATTICE_COORDINATE_H
#define GRID_LATTICE_COORDINATE_H

namespace Grid {

  /*
  depbase=`echo Grid_main.o | sed 's|[^/]*$|.deps/&|;s|\.o$||'`;\
        icpc -DHAVE_CONFIG_H -I. -I../lib    -I../lib -mmic -O3 -std=c++11 -fopenmp -MT Grid_main.o -MD -MP -MF $depbase.Tpo -c -o Grid_main.o Grid_main.cc &&\
        mv -f $depbase.Tpo $depbase.Po
	  ../lib/lattice/Grid_lattice_coordinate.h(25): error: no suitable user-defined conversion from "vector_type" to "const Grid::iScalar<Grid::iScalar<Grid::iScalar<Grid::vInteger>>>" exists
    l._odata[o]=vI;
                    ^
          detected during instantiation of "void Grid::LatticeCoordinate(Grid::Lattice<iobj> &, int) [with iobj=Grid::QCD::vTInteger]" at line 283 of "Grid_main.cc"

	    compilation aborted for Grid_main.cc (code 2)
*/
    template<class iobj> inline void LatticeCoordinate(Lattice<iobj> &l,int mu)
    {
      typedef typename iobj::scalar_object scalar_object;
      typedef typename iobj::scalar_type scalar_type;
      typedef typename iobj::vector_type vector_type;

      GridBase *grid = l._grid;
      int Nsimd = grid->iSites();

      std::vector<int> gcoor;
      std::vector<scalar_type> mergebuf(Nsimd);

      vector_type vI;
      for(int o=0;o<grid->oSites();o++){
	for(int i=0;i<grid->iSites();i++){
	  grid->RankIndexToGlobalCoor(grid->ThisRank(),o,i,gcoor);
	  mergebuf[i]=(Integer)gcoor[mu];
	}
	merge<vector_type,scalar_type>(vI,mergebuf);
	l._odata[o]._internal._internal._internal=vI;
      }
    };

    // LatticeCoordinate();
    // FIXME for debug; deprecate this; made obscelete by 
    template<class vobj> void lex_sites(Lattice<vobj> &l){
      Real *v_ptr = (Real *)&l._odata[0];
      size_t o_len = l._grid->oSites();
      size_t v_len = sizeof(vobj)/sizeof(vRealF);
      size_t vec_len = vRealF::Nsimd();

      for(int i=0;i<o_len;i++){
	for(int j=0;j<v_len;j++){
          for(int vv=0;vv<vec_len;vv+=2){
	    v_ptr[i*v_len*vec_len+j*vec_len+vv  ]= i+vv*500;
	    v_ptr[i*v_len*vec_len+j*vec_len+vv+1]= i+vv*500;
	  }
	}}
    }


}
#endif
