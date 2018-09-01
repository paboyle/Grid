/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/A2AMesonFieldKernels.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#ifndef Hadrons_MContraction_A2AMesonFieldKernels_hpp_
#define Hadrons_MContraction_A2AMesonFieldKernels_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Grid/Eigen/unsupported/CXX11/Tensor>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MContraction)

////////////////////////////////////////////////////////////////////////////////
// Cache blocked arithmetic routine
// Could move to Grid ???
////////////////////////////////////////////////////////////////////////////////
template <typename Field, typename MesonField>
void makeMesonFieldBlock(MesonField &mat, 
                         const Field *lhs_wi,
                         const Field *rhs_vj,
                         std::vector<Gamma::Algebra> gamma,
                         const std::vector<LatticeComplex> &mom,
                         int orthogdim,
                         ModuleBase *caller = nullptr) 
{
    typedef typename Field::vector_object vobj;
    typedef typename vobj::scalar_object  sobj;
    typedef typename vobj::scalar_type    scalar_type;
    typedef typename vobj::vector_type    vector_type;

    typedef iSpinMatrix<vector_type> SpinMatrix_v;
    typedef iSpinMatrix<scalar_type> SpinMatrix_s;
    
    int Lblock = mat.dimension(3); 
    int Rblock = mat.dimension(4);

    GridBase *grid = lhs_wi[0]._grid;
    
    const int    Nd = grid->_ndimension;
    const int Nsimd = grid->Nsimd();

    int Nt     = grid->GlobalDimensions()[orthogdim];
    int Ngamma = gamma.size();
    int Nmom   = mom.size();

    int fd=grid->_fdimensions[orthogdim];
    int ld=grid->_ldimensions[orthogdim];
    int rd=grid->_rdimensions[orthogdim];

    // will locally sum vectors first
    // sum across these down to scalars
    // splitting the SIMD
    int MFrvol = rd*Lblock*Rblock*Nmom;
    int MFlvol = ld*Lblock*Rblock*Nmom;

    Vector<SpinMatrix_v > lvSum(MFrvol);
    parallel_for (int r = 0; r < MFrvol; r++)
    {
        lvSum[r] = zero;
    }

    Vector<SpinMatrix_s > lsSum(MFlvol);             
    parallel_for (int r = 0; r < MFlvol; r++)
    {
        lsSum[r]=scalar_type(0.0);
    }

    int e1=    grid->_slice_nblock[orthogdim];
    int e2=    grid->_slice_block [orthogdim];
    int stride=grid->_slice_stride[orthogdim];

    if (caller) caller->startTimer("contraction: colour trace & mom.");
    // Nested parallelism would be ok
    // Wasting cores here. Test case r
    parallel_for(int r=0;r<rd;r++)
    {
        int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

        for(int n=0;n<e1;n++)
        for(int b=0;b<e2;b++)
        {
            int ss= so+n*stride+b;

            for(int i=0;i<Lblock;i++)
            {
                auto left = conjugate(lhs_wi[i]._odata[ss]);

                for(int j=0;j<Rblock;j++)
                {
                    SpinMatrix_v vv;
                    auto right = rhs_vj[j]._odata[ss];

                    for(int s1=0;s1<Ns;s1++)
                    for(int s2=0;s2<Ns;s2++)
                    {
                        vv()(s1,s2)() = left()(s2)(0) * right()(s1)(0)
                                        + left()(s2)(1) * right()(s1)(1)
                                        + left()(s2)(2) * right()(s1)(2);
                    }
                    
                    // After getting the sitewise product do the mom phase loop
                    int base = Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*r;

                    for ( int m=0;m<Nmom;m++)
                    {
                        int idx = m+base;
                        auto phase = mom[m]._odata[ss];
                        mac(&lvSum[idx],&vv,&phase);
                    }
                }
            }
        }
    }
    if (caller) caller->stopTimer("contraction: colour trace & mom.");

    // Sum across simd lanes in the plane, breaking out orthog dir.
    if (caller) caller->startTimer("contraction: local space sum");
    parallel_for(int rt=0;rt<rd;rt++)
    {
        std::vector<int> icoor(Nd);
        std::vector<SpinMatrix_s> extracted(Nsimd);               

        for(int i=0;i<Lblock;i++)
        for(int j=0;j<Rblock;j++)
        for(int m=0;m<Nmom;m++)
        {

            int ij_rdx = m+Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*rt;

            extract(lvSum[ij_rdx],extracted);
            for(int idx=0;idx<Nsimd;idx++)
            {
                grid->iCoorFromIindex(icoor,idx);

                int ldx    = rt+icoor[orthogdim]*rd;
                int ij_ldx = m+Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*ldx;

                lsSum[ij_ldx]=lsSum[ij_ldx]+extracted[idx];
            }
        }
    }
    if (caller) caller->stopTimer("contraction: local space sum");

    // ld loop and local only??
    if (caller) caller->startTimer("contraction: spin trace");
    int pd = grid->_processors[orthogdim];
    int pc = grid->_processor_coor[orthogdim];
    parallel_for_nest2(int lt=0;lt<ld;lt++)
    {
        for(int pt=0;pt<pd;pt++)
        {
            int t = lt + pt*ld;
            if (pt == pc)
            {
                for(int i=0;i<Lblock;i++)
                for(int j=0;j<Rblock;j++)
                for(int m=0;m<Nmom;m++)
                {
                    int ij_dx = m+Nmom*i + Nmom*Lblock * j + Nmom*Lblock * Rblock * lt;

                    for(int mu=0;mu<Ngamma;mu++)
                    {
                        // this is a bit slow
                        mat(m,mu,t,i,j) = trace(lsSum[ij_dx]*Gamma(gamma[mu]));
                    }
                }
            } 
            else 
            { 
                const scalar_type zz(0.0);

                for(int i=0;i<Lblock;i++)
                for(int j=0;j<Rblock;j++)
                for(int mu=0;mu<Ngamma;mu++)
                for(int m=0;m<Nmom;m++)
                {
                    mat(m,mu,t,i,j) =zz;
                }
            }
        }
    }
    if (caller) caller->stopTimer("contraction: spin trace");
    ////////////////////////////////////////////////////////////////////
    // This global sum is taking as much as 50% of time on 16 nodes
    // Vector size is 7 x 16 x 32 x 16 x 16 x sizeof(complex) = 2MB - 60MB depending on volume
    // Healthy size that should suffice
    ////////////////////////////////////////////////////////////////////
    if (caller) caller->startTimer("contraction: global sum");
    grid->GlobalSumVector(&mat(0,0,0,0,0),Nmom*Ngamma*Nt*Lblock*Rblock);
    if (caller) caller->stopTimer("contraction: global sum");
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif //Hadrons_MContraction_A2AMesonField_hpp_
