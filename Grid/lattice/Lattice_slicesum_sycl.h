#pragma once

NAMESPACE_BEGIN(Grid);

template <class vobj>
inline void sliceSum_sycl(const Lattice<vobj> &Data, std::vector<typename vobj::scalar_object> &result, int orthogdim)
{
    typedef typename vobj::scalar_object sobj;
    typedef typename vobj::scalar_object::scalar_type scalar_type;

    GridBase *grid = Data.Grid();
    assert(grid!=NULL);

    const int Nd = grid->_ndimension;
    const size_t Nsimd = grid->Nsimd();

    assert(orthogdim >= 0);
    assert(orthogdim < Nd);

    int fd=grid->_fdimensions[orthogdim];
    int ld=grid->_ldimensions[orthogdim];
    int rd=grid->_rdimensions[orthogdim];

    int e1=    grid->_slice_nblock[orthogdim];
    int e2=    grid->_slice_block [orthogdim];
    int stride=grid->_slice_stride[orthogdim];
    int ostride=grid->_ostride[orthogdim];
    size_t subvol_size = e1*e2;

    vobj *mysum = (vobj *) malloc_shared(sizeof(vobj),*theGridAccelerator);
    
    result.resize(fd);

    Vector<vobj> lvSum(rd); 
    Vector<sobj> lsSum(ld,Zero());                    
    commVector<vobj> reduction_buffer(rd*subvol_size);
    ExtractBuffer<sobj> extracted(Nsimd);      
    vobj vobj_zero;
    zeroit(vobj_zero);

    for(int r=0;r<rd;r++){
        lvSum[r]=Zero();
    }

    auto rb_p = &reduction_buffer[0];

    autoView(Data_v, Data, AcceleratorRead);

    //prepare reduction buffer 
    accelerator_for2d( s,subvol_size, r,rd, Nsimd,{ 
    
        int n = s / e2;
        int b = s % e2;
        int so=r*ostride; // base offset for start of plane 
        int ss= so+n*stride+b;

        coalescedWrite(rb_p[r*subvol_size+s], coalescedRead(Data_v[ss]));

    });
  
    for (int r = 0; r < rd; r++) {
        mysum[0] = vobj_zero; //dirty hack: cannot pass vobj_zero as identity to sycl::reduction as its not device_copyable
        theGridAccelerator->submit([&](cl::sycl::handler &cgh) {
            auto Reduction = cl::sycl::reduction(mysum,std::plus<>());
            cgh.parallel_for(cl::sycl::range<1>{subvol_size},
            Reduction,
            [=](cl::sycl::id<1> item, auto &sum) {
                auto s = item[0];
                sum += rb_p[r*subvol_size+s];
            });
        });
        theGridAccelerator->wait();
        lvSum[r] = mysum[0];
    }

    Coordinate icoor(Nd);

    for(int rt=0;rt<rd;rt++){

        extract(lvSum[rt],extracted);

        for(int idx=0;idx<Nsimd;idx++){

        grid->iCoorFromIindex(icoor,idx);

        int ldx =rt+icoor[orthogdim]*rd;
        
        lsSum[ldx]=lsSum[ldx]+extracted[idx];

        }
    }

    // sum over nodes.
    for(int t=0;t<fd;t++){
        int pt = t/ld; // processor plane
        int lt = t%ld;
        if ( pt == grid->_processor_coor[orthogdim] ) {
        result[t]=lsSum[lt];
        } else {
        result[t]=Zero();
        }

    }
    scalar_type * ptr = (scalar_type *) &result[0];
    int words = fd*sizeof(sobj)/sizeof(scalar_type);
    grid->GlobalSumVector(ptr, words);

}

NAMESPACE_END(Grid);