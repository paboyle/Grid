#ifndef GRID_LATTICE_COORDINATE_H
#define GRID_LATTICE_COORDINATE_H

namespace Grid {

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
	AmergeA<vector_type,scalar_type>(vI,mergebuf);
	l._odata[o]=vI;
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
