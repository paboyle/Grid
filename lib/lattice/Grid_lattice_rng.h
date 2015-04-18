#ifndef GRID_LATTICE_RNG_H
#define GRID_LATTICE_RNG_H

namespace Grid {

    // FIXME Randomise; deprecate this
    template <class vobj> inline void random(Lattice<vobj> &l){
        Real *v_ptr = (Real *)&l._odata[0];
        size_t v_len = l._grid->oSites()*sizeof(vobj);
        size_t d_len = v_len/sizeof(Real);
	
        for(int i=0;i<d_len;i++){

            v_ptr[i]=drand48();
        }
    };
    
    // FIXME Implement a consistent seed management strategy
    template <class vobj> inline void gaussian(Lattice<vobj> &l){
        // Zero mean, unit variance.
        std::normal_distribution<double> distribution(0.0,1.0);
        Real *v_ptr = (Real *)&l._odata[0];
        size_t v_len = l._grid->oSites()*sizeof(vobj);
        size_t d_len = v_len/sizeof(Real);

        for(int i=0;i<d_len;i++){
	  v_ptr[i]= drand48();
        }
    };

}
#endif
