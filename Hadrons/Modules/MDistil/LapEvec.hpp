/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Hadrons/Modules/MDistil/LapEvec.hpp
 
 Copyright (C) 2019
 
 Author: Felix Erben <ferben@ed.ac.uk>
 Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 
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

#ifndef Hadrons_MDistil_LapEvec_hpp_
#define Hadrons_MDistil_LapEvec_hpp_

#include <Hadrons/Modules/MDistil/Distil.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MDistil)

/******************************************************************************
 
 Laplacian eigenvectors - parameters

 Computes the eigenvectors of the 3D-Laplacian, built from stout-smeared 
 gauge links with the specified number of steps and smearing parameter rho. 
 The smearing is only applied to the spatial components of the gauge field,
 i.e. rho_{4i} = rho_{i4} = rho_{44} = 0. 

 Chebyshev-preconditioning is needed for convergence of the nvec lowest 
 eigenvectors.
 
 ******************************************************************************/

struct StoutParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(StoutParameters,
                                    int, steps,
                                    double, rho)
    StoutParameters() = default;
    template <class ReaderClass> StoutParameters(Reader<ReaderClass>& Reader){read(Reader,"StoutSmearing",*this);}
};

struct ChebyshevParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(ChebyshevParameters,
                                    int, PolyOrder,
                                    double, alpha,
                                    double, beta)
    ChebyshevParameters() = default;
    template <class ReaderClass> ChebyshevParameters(Reader<ReaderClass>& Reader){read(Reader,"Chebyshev",*this);}
};

struct LanczosParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(LanczosParameters,
                                    int, Nvec,
                                    int, Nk,
                                    int, Np,
                                    int, MaxIt,
                                    double, resid,
                                    int, IRLLog)
    LanczosParameters() = default;
    template <class ReaderClass> LanczosParameters(Reader<ReaderClass>& Reader){read(Reader,"Lanczos",*this);}
};

// These are the actual parameters passed to the module during construction

struct LapEvecPar: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(LapEvecPar
                                    ,std::string,         gauge
                                    ,StoutParameters,     Stout
                                    ,ChebyshevParameters, Cheby
                                    ,LanczosParameters,   Lanczos
                                    ,std::string,         FileName)
};

/******************************************************************************
 
 Laplacian eigenvectors - Module (class) definition
 
 ******************************************************************************/

template <typename GImpl>
class TLapEvec: public Module<LapEvecPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
    // constructor
    TLapEvec(const std::string name);
    // destructor
    virtual ~TLapEvec(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
protected:
    std::unique_ptr<GridCartesian> gridLD; // Owned by me, so I must delete it
};

MODULE_REGISTER_TMP(LapEvec, TLapEvec<GIMPL>, MDistil);

/******************************************************************************
 TLapEvec implementation
 ******************************************************************************/

// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TLapEvec<GImpl>::TLapEvec(const std::string name) : Module<LapEvecPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TLapEvec<GImpl>::getInput(void)
{
    return std::vector<std::string>{par().gauge};
}

template <typename GImpl>
std::vector<std::string> TLapEvec<GImpl>::getOutput(void)
{
    return {getName()}; // This is the higher dimensional eigenpack
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TLapEvec<GImpl>::setup(void)
{
    GridCartesian * gridHD = env().getGrid();
    MakeLowerDimGrid(gridLD,gridHD);
    const int Ntlocal{gridHD->LocalDimensions()[Tdir]};
    // Temporaries
    envTmpLat(GaugeField, "Umu_stout");
    envTmpLat(GaugeField, "Umu_smear");
    envTmp(LatticeGaugeField, "UmuNoTime",1,LatticeGaugeField(gridLD.get()));
    envTmp(LatticeColourVector, "src",1,LatticeColourVector(gridLD.get()));
    envTmp(std::vector<LapEvecs>, "eig",1,std::vector<LapEvecs>(Ntlocal));
    // Output objects
    envCreate(LapEvecs, getName(), 1, par().Lanczos.Nvec, gridHD);
}

/*************************************************************************************
 
 -Grad^2 (Peardon, 2009, pg 2, equation 3, https://arxiv.org/abs/0905.2160)
 Field      Type of field the operator will be applied to
 GaugeField Gauge field the operator will smear using
 
 *************************************************************************************/

template<typename Field, typename GaugeField=LatticeGaugeField>
class Laplacian3D : public LinearOperatorBase<Field>, public LinearFunction<Field> {
    typedef typename GaugeField::vector_type vCoeff_t;
public:
    int          nd; // number of spatial dimensions
    std::vector<Lattice<iColourMatrix<vCoeff_t> > > U;
    // Construct this operator given a gauge field and the number of dimensions it should act on
    Laplacian3D( GaugeField& gf, int dimSpatial = Tdir ) : nd{dimSpatial}
    {
        if (dimSpatial<1)
        {
            HADRONS_ERROR(Range,"Must be at least one spatial dimension");
        }
        for (int mu = 0 ; mu < nd ; mu++)
            U.push_back(PeekIndex<LorentzIndex>(gf,mu));
    }
    
    // Apply this operator to "in", return result in "out"
    void operator()(const Field& in, Field& out) {
        if (nd > in.Grid()->Nd())
        {
            HADRONS_ERROR(Range,"nd too large");
        }
        conformable( in, out );
        out = ( ( Real ) ( 2 * nd ) ) * in;
        Field tmp_(in.Grid());
        typedef typename GaugeField::vector_type vCoeff_t;
        for (int mu = 0 ; mu < nd ; mu++)
        {
            out -= U[mu] * Cshift( in, mu, 1);
            tmp_ = adj( U[mu] ) * in;
            out -= Cshift(tmp_,mu,-1);
        }
    }
    
    void OpDiag (const Field &in, Field &out) { HADRONS_ERROR(Definition, "OpDiag() undefined"); };
    void OpDir  (const Field &in, Field &out,int dir,int disp) { HADRONS_ERROR(Definition, "OpDir() undefined"); };
    void Op     (const Field &in, Field &out) { HADRONS_ERROR(Definition, "Op() undefined"); };
    void AdjOp  (const Field &in, Field &out) { HADRONS_ERROR(Definition, "AdjOp() undefined"); };
    void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2) { HADRONS_ERROR(Definition, "HermOpAndNorm() undefined"); };
    void HermOp(const Field &in, Field &out) { operator()(in,out); };
};

template<typename Field>
class Laplacian3DHerm : public LinearFunction<Field> {
public:
    OperatorFunction<Field>   & poly_;
    LinearOperatorBase<Field> &Linop_;
    Laplacian3DHerm(OperatorFunction<Field> & poly,LinearOperatorBase<Field>& linop)
    : poly_{poly}, Linop_{linop} {}
    void operator()(const Field& in, Field& out)
    {
        poly_(Linop_,in,out);
    }
};

/******************************************************************************
 Calculate low-mode eigenvalues of the Laplacian
 ******************************************************************************/

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TLapEvec<GImpl>::execute(void)
{
    const ChebyshevParameters &ChebPar{par().Cheby};
    const LanczosParameters   &LPar{par().Lanczos};
    
    // Disable IRL logging if requested
    LOG(Message) << "IRLLog=" << LPar.IRLLog << std::endl;
    const int PreviousIRLLogState{GridLogIRL.isActive()};
    GridLogIRL.Active( LPar.IRLLog == 0 ? 0 : 1 );
    
    // Stout smearing
    envGetTmp(GaugeField, Umu_smear);
    Umu_smear = envGet(GaugeField, par().gauge); // The smeared field starts off as the Gauge field
    LOG(Message) << "Initial plaquette: " << WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu_smear) << std::endl;
    const StoutParameters &Stout{par().Stout};
    if( Stout.steps )
    {
        envGetTmp(GaugeField, Umu_stout);
        Smear_Stout<PeriodicGimplR> LS(Stout.rho, Tdir); // spatial smearing only
        for (int i = 0; i < Stout.steps; i++) {
            LS.smear(Umu_stout, Umu_smear);
            Umu_smear = Umu_stout;
        }
        LOG(Message) << "Smeared plaquette: " << WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu_smear) << std::endl;
    }
    
    ////////////////////////////////////////////////////////////////////////
    // Invert nabla operator separately on each time-slice
    ////////////////////////////////////////////////////////////////////////
    
    auto & eig4d = envGet(LapEvecs, getName() );
    envGetTmp(std::vector<LapEvecs>, eig);   // Eigenpack for each timeslice
    envGetTmp(LatticeGaugeField, UmuNoTime); // Gauge field without time dimension
    envGetTmp(LatticeColourVector, src);
    GridCartesian * gridHD = env().getGrid();
    const int Ntlocal{gridHD->LocalDimensions()[Tdir]};
    const int Ntfirst{gridHD->LocalStarts()[Tdir]};
    uint32_t ConvergenceErrors{0};
    const int NtFull{env().getDim(Tdir)};
    TimesliceEvals Evals{ NtFull, LPar.Nvec };
    for (int t = 0; t < NtFull; t++)
        for (int v = 0; v < LPar.Nvec; v++)
            Evals.tensor( t, v ) = 0;
    for (int t = 0; t < Ntlocal; t++ )
    {
        LOG(Message) << "------------------------------------------------------------" << std::endl;
        LOG(Message) << " Compute eigenpack, local timeslice = " << t << " / " << Ntlocal << std::endl;
        LOG(Message) << "------------------------------------------------------------" << std::endl;
        eig[t].resize(LPar.Nk+LPar.Np,gridLD.get());
        
        // Construct smearing operator
        ExtractSliceLocal(UmuNoTime,Umu_smear,0,t,Tdir); // switch to 3d/4d objects
        Laplacian3D<LatticeColourVector> Nabla(UmuNoTime);
        LOG(Message) << "Chebyshev preconditioning to order " << ChebPar.PolyOrder
                     << " with parameters (alpha,beta) = (" << ChebPar.alpha << "," << ChebPar.beta << ")" << std::endl;
        Chebyshev<LatticeColourVector> Cheb(ChebPar.alpha,ChebPar.beta,ChebPar.PolyOrder);
        
        // Construct source vector according to Test_dwf_compressed_lanczos.cc
        src = 11.0; // NB: This is a dummy parameter and just needs to be non-zero
        RealD nn = norm2(src);
        nn = Grid::sqrt(nn);
        src = src * (1.0/nn);
        
        Laplacian3DHerm<LatticeColourVector> NablaCheby(Cheb,Nabla);
        ImplicitlyRestartedLanczos<LatticeColourVector>
        IRL(NablaCheby,Nabla,LPar.Nvec,LPar.Nk,LPar.Nk+LPar.Np,LPar.resid,LPar.MaxIt);
        int Nconv = 0;
        IRL.calc(eig[t].eval,eig[t].evec,src,Nconv);
        if (Nconv < LPar.Nvec)
        {
            // NB: Can't assert here since we are processing local slices - i.e. not all nodes would assert
            ConvergenceErrors = 1;
            LOG(Error) << "MDistil::LapEvec : Not enough eigenvectors converged. If this occurs in practice, we should modify the eigensolver to iterate once more to ensure the second convergence test does not take us below the requested number of eigenvectors" << std::endl;
        }
        if( Nconv != LPar.Nvec )
            eig[t].resize(LPar.Nvec, gridLD.get());
        RotateEigen( eig[t].evec ); // Rotate the eigenvectors into our phase convention
        
        for (int i=0;i<LPar.Nvec;i++){
            InsertSliceLocal(eig[t].evec[i],eig4d.evec[i],0,t,Tdir);
            if(t==0 && Ntfirst==0)
                eig4d.eval[i] = eig[t].eval[i]; // TODO: Discuss: is this needed? Is there a better way?
            if(gridLD->IsBoss()) // Only do this on one node per timeslice, so a global sum will work
                Evals.tensor(t + Ntfirst,i) = eig[t].eval[i];
        }
    }
    GridLogIRL.Active( PreviousIRLLogState );
    gridHD->GlobalSum(ConvergenceErrors);
    if(ConvergenceErrors!=0)
    {
        HADRONS_ERROR(Program,"The eingensolver failed to find enough eigenvectors on at least one node");
    }
    // Now write out the 4d eigenvectors
    std::string sEigenPackName(par().FileName);
    if( !sEigenPackName.empty() )
    {
        eig4d.record.solverXml = parString();
        ModuleBase * b{vm().getModule(par().gauge)};
        std::string sOperatorXml{ "<module><id><type>" };
        sOperatorXml.append( b->getRegisteredName() );
        sOperatorXml.append( "</type></id><options>" );
        sOperatorXml.append( b->parString() );
        sOperatorXml.append( "</options></module>" );
        eig4d.record.operatorXml = sOperatorXml;
        sEigenPackName.append(1, '.');
        std::size_t NameLen{ sEigenPackName.length() };
        const std::string sTrajNum{std::to_string(vm().getTrajectory())};
        sEigenPackName.append(sTrajNum);
        eig4d.write(sEigenPackName,false);
        // Communicate eig[t].evec to boss-node, save into new object evecs
        gridHD->GlobalSumVector(EigenIO::getFirstScalar(Evals.tensor),
                                static_cast<int>(EigenIO::getScalarCount(Evals.tensor)));
        if(gridHD->IsBoss())
        {
            sEigenPackName.resize(NameLen);
            sEigenPackName.append("evals.");
            sEigenPackName.append(sTrajNum);
            Evals.write( sEigenPackName );
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_LapEvec_hpp_
