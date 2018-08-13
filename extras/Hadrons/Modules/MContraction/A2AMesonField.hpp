/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MContraction/A2AMesonField.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#ifndef Hadrons_MContraction_A2AMesonField_hpp_
#define Hadrons_MContraction_A2AMesonField_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/A2AVectors.hpp>
#include <Grid/Eigen/unsupported/CXX11/Tensor>
#include <Grid/Hadrons/Modules/MSolver/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     All-to-all meson field creation                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaPair;


class A2AMesonFieldPar: Serializable
{
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMesonFieldPar,
                                    int, cacheBlock,
                                    int, block,
                                    std::string, v,
                                    std::string, w,
                                    std::string, output,
                                    std::string, gammas,
                                    std::vector<std::string>, mom);
};

template <typename FImpl>
class TA2AMesonField : public Module<A2AMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef Eigen::TensorMap<Eigen::Tensor<Complex, 5, Eigen::RowMajor>> MesonField;
public:
    // constructor
    TA2AMesonField(const std::string name);
    // destructor
    virtual ~TA2AMesonField(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    // Arithmetic kernel. Move to Grid??
    void makeBlock(MesonField &mat, 
                   const FermionField *lhs,
                   const FermionField *rhs,
                   std::vector<Gamma::Algebra> gamma,
                   const std::vector<LatticeComplex> &mom,
                   int orthogdim,
                   double &t0,
                   double &t1,
                   double &t2,
                   double &t3);
    // IO
    std::string ioname(unsigned int m, unsigned int g) const;
    std::string filename(unsigned int m, unsigned int g) const;
    void initFile(unsigned int m, unsigned int g);
    void saveBlock(const MesonField &mf,
                   unsigned int m, unsigned int g, 
                   unsigned int i, unsigned int j);
private:
    bool                             hasPhase_{false};
    std::string                      momphName_;
    std::vector<Gamma::Algebra>      gamma_;
    std::vector<std::vector<double>> mom_;
};

MODULE_REGISTER(A2AMesonField, ARG(TA2AMesonField<FIMPL>), MContraction);
MODULE_REGISTER(ZA2AMesonField, ARG(TA2AMesonField<ZFIMPL>), MContraction);

/******************************************************************************
*                  TA2AMesonField implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AMesonField<FImpl>::TA2AMesonField(const std::string name)
: Module<A2AMesonFieldPar>(name)
, momphName_(name + "_momph")
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AMesonField<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().v, par().w};

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2AMesonField<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AMesonField<FImpl>::setup(void)
{
    gamma_.clear();
    mom_.clear();
    if (par().gammas == "all")
    {
    gamma_ = {
        Gamma::Algebra::Gamma5,
        Gamma::Algebra::Identity,    
        Gamma::Algebra::GammaX,
        Gamma::Algebra::GammaY,
        Gamma::Algebra::GammaZ,
        Gamma::Algebra::GammaT,
        Gamma::Algebra::GammaXGamma5,
        Gamma::Algebra::GammaYGamma5,
        Gamma::Algebra::GammaZGamma5,
        Gamma::Algebra::GammaTGamma5,
        Gamma::Algebra::SigmaXY,
        Gamma::Algebra::SigmaXZ,
        Gamma::Algebra::SigmaXT,
        Gamma::Algebra::SigmaYZ,
        Gamma::Algebra::SigmaYT,
        Gamma::Algebra::SigmaZT
    };
    }
    else
    {
        gamma_ = strToVec<Gamma::Algebra>(par().gammas);
    }
    for (auto &pstr: par().mom)
    {
        auto p = strToVec<Real>(pstr);

        if (p.size() != env().getNd() - 1)
        {
            HADRONS_ERROR(Size, "Momentum has " + std::to_string(p.size())
                                + " components instead of " 
                                + std::to_string(env().getNd() - 1));
        }
        mom_.push_back(p);
    }
    
    envCache(std::vector<LatticeComplex>, momphName_, 1, 
             par().mom.size(), env().getGrid());
    envTmpLat(LatticeComplex, "coor");
    // preallocate memory for meson field block
    auto tgp = env().getDim().back()*gamma_.size()*mom_.size();

    envTmp(Vector<Complex>, "mfBuf", 1, tgp*par().block*par().block);
    envTmp(Vector<Complex>, "mfCache", 1, tgp*par().cacheBlock*par().cacheBlock);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AMesonField<FImpl>::execute(void)
{
    auto &v = envGet(std::vector<FermionField>, par().v);
    auto &w = envGet(std::vector<FermionField>, par().w);

    int nt         = env().getDim().back();
    int N_i        = w.size();
    int N_j        = v.size();
    int ngamma     = gamma_.size();
    int nmom       = mom_.size();
    int block      = par().block;
    int cacheBlock = par().cacheBlock;

    LOG(Message) << "Computing all-to-all meson fields" << std::endl;
    LOG(Message) << "W: '" << par().w << "' V: '" << par().v << "'" << std::endl;
    LOG(Message) << "Meson field size: " << nt << "*" << N_i << "*" << N_j 
                 << " (" << sizeString(nt*N_i*N_j*sizeof(Complex)) << ")" << std::endl;
    LOG(Message) << "Momenta:" << std::endl;
    for (auto &p: mom_)
    {
        LOG(Message) << "  " << p << std::endl;
    }
    LOG(Message) << "Spin structures:" << std::endl;
    for (auto &g: gamma_)
    {
        LOG(Message) << "  " << g << std::endl;
    }

    ///////////////////////////////////////////////
    // Momentum setup
    ///////////////////////////////////////////////
    auto &ph = envGet(std::vector<LatticeComplex>, momphName_);

    if (!hasPhase_)
    {
        MODULE_TIMER("Momentum phases");
        for (unsigned int j = 0; j < nmom; ++j)
        {
            Complex           i(0.0,1.0);
            std::vector<Real> p;

            envGetTmp(LatticeComplex, coor);
            ph[j] = zero;
            for(unsigned int mu = 0; mu < mom_[j].size(); mu++)
            {
                LatticeCoordinate(coor, mu);
                ph[j] = ph[j] + (mom_[j][mu]/env().getDim(mu))*coor;
            }
            ph[j] = exp((Real)(2*M_PI)*i*ph[j]);
        }
        hasPhase_ = true;
    }
    

    //////////////////////////////////////////////////////////////////////////
    // i,j   is first  loop over SchurBlock factors reusing 5D matrices
    // ii,jj is second loop over cacheBlock factors for high perf contractoin
    // iii,jjj are loops within cacheBlock
    // Total index is sum of these  i+ii+iii etc...
    //////////////////////////////////////////////////////////////////////////
    
    double flops = 0.0;
    double bytes = 0.0;
    double vol   = env().getVolume();
    double t_schur=0;
    double t_contr=0;
    double t_int_0=0;
    double t_int_1=0;
    double t_int_2=0;
    double t_int_3=0;

    envGetTmp(Vector<Complex>, mfBuf);
    envGetTmp(Vector<Complex>, mfCache);
    
    double t0    = usecond();
    int NBlock_i = N_i/block + (((N_i % block) != 0) ? 1 : 0);
    int NBlock_j = N_j/block + (((N_j % block) != 0) ? 1 : 0);

    for(int i=0;i<N_i;i+=block)
    for(int j=0;j<N_j;j+=block)
    {
        // Get the W and V vectors for this block^2 set of terms
        int N_ii = MIN(N_i-i,block);
        int N_jj = MIN(N_j-j,block);

        t_schur-=usecond();
        t_schur+=usecond();

        LOG(Message) << "Meson field block " 
                    << j/block + NBlock_j*i/block + 1 
                    << "/" << NBlock_i*NBlock_j << " [" << i <<" .. " 
                    << i+N_ii-1 << ", " << j <<" .. " << j+N_jj-1 << "]" 
                    << std::endl;

        MesonField mfBlock(mfBuf.data(),nmom,ngamma,nt,N_ii,N_jj);

        // Series of cache blocked chunks of the contractions within this block
        for(int ii=0;ii<N_ii;ii+=cacheBlock)
        for(int jj=0;jj<N_jj;jj+=cacheBlock)
        {
            int N_iii = MIN(N_ii-ii,cacheBlock);
            int N_jjj = MIN(N_jj-jj,cacheBlock);
            MesonField mfCacheBlock(mfCache.data(),nmom,ngamma,nt,N_iii,N_jjj);    

            t_contr-=usecond();
            makeBlock(mfCacheBlock, &w[i+ii], &v[j+jj], gamma_, ph, 
                      env().getNd() - 1, t_int_0, t_int_1, t_int_2, t_int_3);
            t_contr+=usecond();
            
            // flops for general N_c & N_s
            flops += vol * ( 2 * 8.0 + 6.0 + 8.0*nmom) * N_iii*N_jjj*ngamma;
            bytes += vol * (12.0 * sizeof(Complex) ) * N_iii*N_jjj
                +  vol * ( 2.0 * sizeof(Complex) *nmom ) * N_iii*N_jjj* ngamma;

            MODULE_TIMER("Cache copy");
            for(int iii=0;iii< N_iii;iii++)
            for(int jjj=0;jjj< N_jjj;jjj++)
            for(int m =0;m< nmom;m++)
            for(int g =0;g< ngamma;g++)
            for(int t =0;t< nt;t++)
            {
                mfBlock(m,g,t,ii+iii,jj+jjj) = mfCacheBlock(m,g,t,iii,jjj);
            }
        }

        // IO
        double blockSize, ioTime;
        
        MODULE_TIMER("IO");
        for(int m = 0; m < nmom; m++)
        for(int g = 0; g < ngamma; g++)
        {
            if ((i == 0) and (j == 0))
            {
                initFile(m, g);
            }
            saveBlock(mfBlock, m, g, i, j);
        }
        blockSize = static_cast<double>(nmom*ngamma*nt*N_ii*N_jj*sizeof(Complex));
        ioTime    = static_cast<double>(this->getTimer("IO").count());
        LOG(Message) << "HDF5 IO " << blockSize/ioTime*1.0e6/1024/1024
                     << " MB/s" << std::endl;
    }

    double nodes    = env().getGrid()->NodeCount();
    double t_kernel = t_int_0 + t_int_1;

    LOG(Message) << "Perf " << flops/(t_kernel)/1.0e3/nodes << " Gflop/s/node "  << std::endl;
    LOG(Message) << "Perf " << bytes/(t_kernel)/1.0e3/nodes << " GB/s/node "  << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////
// Cache blocked arithmetic routine
// Could move to Grid ???
//////////////////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AMesonField<FImpl>::makeBlock(MesonField &mat, 
					                  const FermionField *lhs_wi,
					                  const FermionField *rhs_vj,
					                  std::vector<Gamma::Algebra> gamma,
					                  const std::vector<LatticeComplex> &mom,
					                  int orthogdim,
					                  double &t0,
					                  double &t1,
					                  double &t2,
					                  double &t3) 
{
    typedef typename FImpl::SiteSpinor vobj;

    typedef typename vobj::scalar_object sobj;
    typedef typename vobj::scalar_type scalar_type;
    typedef typename vobj::vector_type vector_type;

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

    t0-=usecond();
    MODULE_TIMER("Colour trace * mom.");
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
    t0+=usecond();

    // Sum across simd lanes in the plane, breaking out orthog dir.
    MODULE_TIMER("Local space sum");
    t1-=usecond();
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
    t1+=usecond();
    assert(mat.dimension(0) == Nmom);
    assert(mat.dimension(1) == Ngamma);
    assert(mat.dimension(2) == Nt);
    t2-=usecond();

    // ld loop and local only??
    MODULE_TIMER("Spin trace");
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
    t2+=usecond();
    ////////////////////////////////////////////////////////////////////
    // This global sum is taking as much as 50% of time on 16 nodes
    // Vector size is 7 x 16 x 32 x 16 x 16 x sizeof(complex) = 2MB - 60MB depending on volume
    // Healthy size that should suffice
    ////////////////////////////////////////////////////////////////////
    t3-=usecond();
    MODULE_TIMER("Global sum");
    grid->GlobalSumVector(&mat(0,0,0,0,0),Nmom*Ngamma*Nt*Lblock*Rblock);
    t3+=usecond();
}

// IO
template <typename FImpl>
std::string TA2AMesonField<FImpl>::ioname(unsigned int m, unsigned int g) const
{
    std::stringstream ss;

    ss << gamma_[g] << "_";
    for (unsigned int mu = 0; mu < mom_[m].size(); ++mu)
    {
        ss << mom_[m][mu] << ((mu == mom_[m].size() - 1) ? "" : "_");
    }

    return ss.str();
}

template <typename FImpl>
std::string TA2AMesonField<FImpl>::filename(unsigned int m, unsigned int g) const
{
    return par().output + "." + std::to_string(vm().getTrajectory()) 
           + "/" + ioname(m, g) + ".h5";
}

template <typename FImpl>
void TA2AMesonField<FImpl>::initFile(unsigned int m, unsigned int g)
{
#ifdef HAVE_HDF5
    std::string  f     = filename(m, g);
    GridBase     *grid = env().getGrid();
    auto         &v    = envGet(std::vector<FermionField>, par().v);
    auto         &w    = envGet(std::vector<FermionField>, par().w);
    int          nt    = env().getDim().back();
    int          N_i   = w.size();
    int          N_j   = v.size();

    makeFileDir(f, grid);
    if (grid->IsBoss())
    {
        Hdf5Writer           writer(f);
        std::vector<hsize_t> dim = {static_cast<hsize_t>(nt), 
                                    static_cast<hsize_t>(N_i), 
                                    static_cast<hsize_t>(N_j)};
        H5NS::DataSpace      dataspace(dim.size(), dim.data());
        H5NS::DataSet        dataset;
        
        push(writer, ioname(m, g));
        write(writer, "momentum", mom_[m]);
        write(writer, "gamma", gamma_[g]);
        auto &group = writer.getGroup();
        dataset = group.createDataSet("mesonField", Hdf5Type<Complex>::type(),
                                      dataspace);
    }
#else
    HADRONS_ERROR(Implementation, "meson field I/O needs HDF5 library");
#endif
}

template <typename FImpl>
void TA2AMesonField<FImpl>::saveBlock(const MesonField &mf,
                                      unsigned int m, unsigned int g, 
                                      unsigned int i, unsigned int j)
{
#ifdef HAVE_HDF5
    std::string  f     = filename(m, g);
    GridBase     *grid = env().getGrid();

    if (grid->IsBoss())
    {
        Hdf5Reader           reader(f);
        hsize_t              nt = mf.dimension(2),
                             Ni = mf.dimension(3),
                             Nj = mf.dimension(4);
        std::vector<hsize_t> count = {nt, Ni, Nj},
                             offset = {0, static_cast<hsize_t>(i),
                                       static_cast<hsize_t>(j)},
                             stride = {1, 1, 1},
                             block  = {1, 1, 1}; 
        H5NS::DataSpace      memspace(count.size(), count.data()), dataspace;
        H5NS::DataSet        dataset;
        size_t               shift;

        push(reader, ioname(m, g));
        auto &group = reader.getGroup();
        dataset   = group.openDataSet("mesonField");
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                                  stride.data(), block.data());
        shift = (m*mf.dimension(1) + g)*nt*Ni*Nj;
        dataset.write(mf.data() + shift, Hdf5Type<Complex>::type(), memspace, 
                      dataspace);
    }
#else
    HADRONS_ERROR(Implementation, "meson field I/O needs HDF5 library");
#endif
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AMesonField_hpp_
