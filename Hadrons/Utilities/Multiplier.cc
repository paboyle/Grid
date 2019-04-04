#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DiskVector.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/Utilities/Multiplier.hpp>


#ifdef GRID_COMMS_MPI3
#define GET_RANK(rank, nMpi) \
MPI_Comm_size(MPI_COMM_WORLD, &(nMpi));\
MPI_Comm_rank(MPI_COMM_WORLD, &(rank))
#define BARRIER() MPI_Barrier(MPI_COMM_WORLD)
#define GLOBAL_DSUM(x) MPI_Allreduce(MPI_IN_PLACE, &x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD)
#define GLOBAL_DMAX(x) MPI_Allreduce(MPI_IN_PLACE, &x, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD)
#define INIT() MPI_Init(NULL, NULL)
#define FINALIZE() MPI_Finalize()
#else
#define GET_RANK(rank, nMpi) (nMpi) = 1; (rank) = 0
#define BARRIER()
#define GLOBAL_DSUM(x)
#define GLOBAL_DMAX(x)
#define INIT()
#define FINALIZE()
#endif

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

#ifndef HADRONS_A2AM_IO_TYPE
#define HADRONS_A2AM_IO_TYPE ComplexF
#endif

#define TIME_MOD(t) (((t) + par.global.nt) % par.global.nt)

struct MultiplierPar
{
    Multiplier::GlobalPar                  global;
    std::vector<Multiplier::A2AMatrixPar>  a2aMatrix;
    std::vector<Multiplier::ProductPar>    product;
};

struct IoHelper
{
    A2AMatrixIo<HADRONS_A2AM_IO_TYPE> io;
    Multiplier::MultiplierMetadata md;
    unsigned int e, s, i, j;
};

void makeTimeSeq(std::vector<std::vector<unsigned int>> &timeSeq,
                 const std::vector<std::set<unsigned int>> &times,
                 std::vector<unsigned int> &current,
                 const unsigned int depth)
{
    if (depth > 0)
    {
        for (auto t : times[times.size() - depth])
        {
            current[times.size() - depth] = t;
            makeTimeSeq(timeSeq, times, current, depth - 1);
        }
    }
    else
    {
        timeSeq.push_back(current);
    }
}

void makeTimeSeq(std::vector<std::vector<unsigned int>> &timeSeq,
                 const std::vector<std::set<unsigned int>> &times)
{
    std::vector<unsigned int> current(times.size());

    makeTimeSeq(timeSeq, times, current, times.size());
}

void makeOutFilename(std::string &outFilename, std::stringstream &ioName,
                     std::string momentum, std::string &outStem, std::vector<std::string> term,
                     const std::vector<Gamma::Algebra> &gammas, const int traj,
                     const std::vector<unsigned int> &dt)
{
    std::string termString;

    if (gammas.size() != dt.size())
            {
                HADRONS_ERROR(Size, "number of gammas (" + std::to_string(gammas.size()) 
                            + ") different from number of times (" 
                            + std::to_string(dt.size()) + ")");
            }
    for (unsigned int g = 0; g < gammas.size(); g++)
    {
        if (g < gammas.size() - 1)
        {
            ioName << std::to_string(dt[g]) << "_" << gammas[g] << "_";
        }
        else
        {
            ioName << std::to_string(dt[g]) << "_" << gammas[g];
        }
        
        
    }
    ioName << momentum;

    for(unsigned int i = 0; i < term.size(); i++)
    {
        if ( i < term.size() - 1)
        {
            termString += term[i] + "_";
        }
        else
        {
            termString += term[i];
        }
        
    }
    
    outFilename = outStem + "." + std::to_string(traj) + "/" + termString + "/" + ioName.str() + ".h5";
}

void makeIoHelper(IoHelper &helper, std::string momentum, std::string &outStem, std::vector<std::string> term,
                  const std::vector<Gamma::Algebra> &gammas,
                  const std::vector<unsigned int> dt, const int nt, const int N_i, const int N_j,
                  const int blockSize, const int traj)
{
    std::string filenameOut;
    std::stringstream ioName;
    makeOutFilename(filenameOut, ioName, momentum, outStem, term, gammas, traj, dt);
    helper.io = A2AMatrixIo<HADRONS_A2AM_IO_TYPE>(filenameOut,
                                             ioName.str(), nt, N_i, N_j);

    std::replace(momentum.begin(), momentum.end(), '_', ' ');
    std::vector<std::vector<Real>> mom;
    auto p = strToVec<Real>(momentum);
    mom.push_back(p);

    helper.md.momentum = mom[0];
    helper.md.gammas = gammas;

    std::cout << "Saving matrix product to '" << filenameOut << "'" << std::endl;
    makeFileDir(filenameOut);
    helper.io.initFile(helper.md, blockSize);
}

void parseDataset(std::string &dataset, std::string &matrix, std::string &momentum)
{
    std::regex matRex("(^[^_]+)");
    std::regex momRex("_.*$");

    std::smatch smMatrix;
    std::smatch smMom;
    std::string s;

    std::regex_search(dataset, smMatrix, matRex);
    matrix = smMatrix[0].str();
    std::regex_search(dataset, smMom, momRex);
    momentum = smMom[0].str();
}

void printPerf(const double bytes, const double usec)
{
    double maxt;

    maxt = usec;
    GLOBAL_DMAX(maxt);
    std::cout << maxt/1.0e6 << " sec " << bytes/maxt*1.0e6/1024/1024/1024 << " GB/s";
}

void printPerf(const double bytes, const double busec, 
               const double flops, const double fusec)
{
    double maxt;

    printPerf(bytes, busec);
    std::cout << " ";
    maxt = fusec;
    GLOBAL_DMAX(maxt);
    std::cout << flops/fusec/1.0e3 << " GFlop/s";
}

std::set<unsigned int> parseTimeRange(const std::string str, const unsigned int nt)
{
    std::regex rex("([0-9]+)|(([0-9]+)\\.\\.([0-9]+))");
    std::smatch sm;
    std::vector<std::string> rstr = strToVec<std::string>(str);
    std::set<unsigned int> tSet;

    for (auto &s : rstr)
    {
        std::regex_match(s, sm, rex);
        if (sm[1].matched)
        {
            unsigned int t;

            t = std::stoi(sm[1].str());
            if (t >= nt)
            {
                HADRONS_ERROR(Range, "time out of range (from expression '" + str + "')");
            }
            tSet.insert(t);
        }
        else if (sm[2].matched)
        {
            unsigned int ta, tb;

            ta = std::stoi(sm[3].str());
            tb = std::stoi(sm[4].str());
            if ((ta >= nt) or (tb >= nt))
            {
                HADRONS_ERROR(Range, "time out of range (from expression '" + str + "')");
            }
            for (unsigned int ti = ta; ti <= tb; ++ti)
            {
                tSet.insert(ti);
            }
        }
    }

    return tSet;
}
struct Sec
{
    Sec(const double usec)
    {
        seconds = usec/1.0e6;
    }
    
    double seconds;
};

inline std::ostream & operator<< (std::ostream& s, const Sec &&sec)
{
    s << std::setw(10) << sec.seconds << " sec";

    return s;
}

struct Flops
{
    Flops(const double flops, const double fusec)
    {
        gFlopsPerSec = flops/fusec/1.0e3;
    }
    
    double gFlopsPerSec;
};

inline std::ostream & operator<< (std::ostream& s, const Flops &&f)
{
    s << std::setw(10) << f.gFlopsPerSec << " GFlop/s";

    return s;
}

struct Bytes
{
    Bytes(const double bytes, const double busec)
    {
        gBytesPerSec = bytes/busec*1.0e6/1024/1024/1024;
    }
    
    double gBytesPerSec;
};

inline std::ostream & operator<< (std::ostream& s, const Bytes &&b)
{
    s << std::setw(10) << b.gBytesPerSec << " GB/s";

    return s;
}

int main(int argc, char* argv[])
{
    // MPI init
    int nMpi, rank;

    INIT();
    GET_RANK(rank, nMpi);
    if (rank != 0)
    {
        std::cout.setstate(std::ios::badbit);
    }

    // parse command line
    std::string   parFilename;

    if (argc != 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file>";
        std::cerr << std::endl;
        
        return EXIT_FAILURE;
    }
    parFilename = argv[1];

    // parse parameter file
    MultiplierPar par;
    unsigned int  nMat, nCont;
    XmlReader     reader(parFilename);

    read(reader, "global",    par.global);
    read(reader, "a2aMatrix", par.a2aMatrix);
    read(reader, "product",   par.product);
    nMat  = par.a2aMatrix.size();
    nCont = par.product.size();

    // create diskvectors
    std::map<std::string, EigenDiskVector<ComplexD>> a2aMat;
    std::map<std::string, std::string>               a2aMatDataset;
    unsigned int                                     cacheSize;
    bool                                             clean = true;
    std::string                                      prodFileName = "product";

    for (auto &p: par.a2aMatrix)
    {
        if (p.name == prodFileName)
        {
            HADRONS_ERROR(Implementation, "Parameter name must differ from product file name: '"
                                          + prodFileName + "'");
        }
        std::string dirName = par.global.diskVectorDir + "/" + p.name + "." + std::to_string(rank);

        a2aMat.emplace(p.name, EigenDiskVector<ComplexD>(dirName, par.global.nt, par.global.cacheSize, clean));
        a2aMatDataset.emplace(p.name, p.dataset);
    }

    std::string dirName = par.global.diskVectorDir + "/" + prodFileName + "." + std::to_string(rank);
    EigenDiskVector<ComplexD> prod(dirName, par.global.nt, par.global.cacheSize, clean);

    // trajectory loop
    std::vector<unsigned int> tList = par.global.trajCounter.getTrajectoryList();
    unsigned int              indi, inde, indPerRank;

    indPerRank = tList.size()/nMpi;
    indi       = rank*indPerRank;
    BARRIER();
    std::cout << ":::::::: nMpi " << nMpi << std::endl;
    for (unsigned int tInd = indi; tInd < indi + indPerRank; tInd++)
    {
        unsigned int traj;

        if (tInd < tList.size())
        {
            traj = tList[tInd];
        }
        else
        {
            traj = tList.back();
        }
        if (nMpi > 1)
        {
            if (rank == 0)
            {
                std::cout << ":::::::: Trajectories ";
                for (unsigned int r = 0; r < nMpi - 1; ++r)
                {
                    std::cout << tList[tInd + r*indPerRank] << " ";
                }
                if (tInd + (nMpi - 1)*indPerRank < tList.size())
                {
                    std::cout << tList[tInd + (nMpi - 1)*indPerRank];
                }
                std::cout << std::endl;
            }
        }
        else
        {
            std::cout << ":::::::: Trajectory " << traj << std::endl;
        }
    
        int nt = par.global.nt;
        int blockSize = par.global.blockSize;

        Vector<HADRONS_A2AM_IO_TYPE> mfBuf;
        int mfBufSize = nt * blockSize * blockSize;
        mfBuf.resize(mfBufSize);

        // load data
        for (auto &p: par.a2aMatrix)
        {
            std::string filename = p.file;
            double      t, size;

            tokenReplace(filename, "traj", traj);
            std::cout << "======== Loading '" << filename << "'" << std::endl;

            A2AMatrixIo<HADRONS_A2AM_IO_TYPE> a2aIo(filename, p.dataset, nt);

            a2aIo.load(a2aMat.at(p.name), &t);
            std::cout << "Read " << a2aIo.getSize() << " bytes in " << t/1.0e6 
                    << " sec, " << a2aIo.getSize()/t*1.0e6/1024/1024 << " MB/s" << std::endl;
        }

        // Multiply
        for (auto &p: par.product)
        {
            std::vector<std::string>               term = strToVec<std::string>(p.terms);
            std::vector<std::set<unsigned int>>    times;
            std::vector<std::vector<unsigned int>> timeSeq;
            A2AMatrix<ComplexD>                    prodT, tmp;
            TimerArray                             tAr;
            double                                 fusec, busec, flops, bytes, tusec;

            tAr.startTimer("Total");
            if (term.size() != p.times.size())
            {
                HADRONS_ERROR(Size, "number of terms (" + std::to_string(term.size()) 
                            + ") different from number of times (" 
                            + std::to_string(p.times.size()) + ")");
            }
            for (auto &s : p.times)
            {
                times.push_back(parseTimeRange(s, nt));
            }

            makeTimeSeq(timeSeq, times);
            std::vector<unsigned int> dt = timeSeq[0];

            // A2A Eigen Tensor
            prodT = a2aMat.at(term[0])[0];
            int N_i = prodT.rows();
            prodT = a2aMat.at(term[term.size() - 1])[0];
            int N_j = prodT.cols();

            std::vector<Gamma::Algebra> gammas, gamma;
            std::string matString, momString;
            std::vector<std::string> momStrings;

            for (unsigned int j = 0; j < term.size(); ++j)
            {
                parseDataset(a2aMatDataset.at(term[j]), matString, momString);

                gamma = strToVec<Gamma::Algebra>(matString);
                gammas.push_back(gamma[0]);

                momStrings.push_back(momString);
                
                if (momStrings[0] != momStrings[j])
                {
                    HADRONS_ERROR(Implementation, "Momentum mismatch between term 0: (" + momStrings[0] +
                     ") and term " + std::to_string(j) + ": (" + momStrings[j] + ").");
                }
            }

            IoHelper h;
            makeIoHelper(h, momStrings[0], par.global.output, term, gammas, dt, nt, N_i, N_j, blockSize, traj);

            flops  = 0.;
            bytes  = 0.;
            fusec  = tAr.getDTimer("A*B algebra");
            busec  = tAr.getDTimer("A*B total");
            tAr.startTimer("Linear algebra");
            for (unsigned int t = 0; t < nt; t++)
            {
                tAr.startTimer("Disk vector overhead");
                prodT = a2aMat.at(term[0])[TIME_MOD(t + dt[0])];
                tAr.stopTimer("Disk vector overhead");
                for (unsigned int j = 1; j < term.size(); ++j)
                {
                    tAr.startTimer("Disk vector overhead");
                    const A2AMatrix<ComplexD> &ref = a2aMat.at(term[j])[TIME_MOD(t + dt[j])];
                    tAr.stopTimer("Disk vector overhead");

                    tAr.startTimer("A*B total");
                    tAr.startTimer("A*B algebra");
                    A2AContraction::mul(tmp, prodT, ref);
                    tAr.stopTimer("A*B algebra");
                    flops += A2AContraction::mulFlops(prodT, ref);
                    prodT = tmp;
                    tAr.stopTimer("A*B total");
                    bytes += 3. * tmp.rows() * tmp.cols() * sizeof(ComplexD);
                }
                if (term.size() > 1)
                {
                    printPerf(bytes*nMpi, tAr.getDTimer("A*B total") - busec,
                              flops*nMpi, tAr.getDTimer("A*B algebra") - fusec);
                    std::cout << std::endl;
                }
                prod[t] = prodT;
            }
            tAr.stopTimer("Linear algebra");

            tAr.startTimer("Matrix I/o");
            for (unsigned int i = 0; i < N_i; i += blockSize){
                std::cout << "i = " << i << " of " << N_i << std::endl;
                for (unsigned int j = 0; j < N_j; j += blockSize)
                {
                    int N_ii = MIN(N_i - i, blockSize);
                    int N_jj = MIN(N_j - j, blockSize);
                    A2AMatrixSet<HADRONS_A2AM_IO_TYPE> mfBlock(mfBuf.data(), 1, 1, nt, N_ii, N_jj);
                    tAr.startTimer("Matrix copy");
                    for (unsigned int t = 0; t < nt; t++)
                        for (unsigned int ii = 0; ii < N_ii; ii++)
                            for (unsigned int jj = 0; jj < N_jj; jj++)
                            {
                                mfBlock(0, 0, t, ii, jj) = prod(t,i + ii, j + jj);
                            }
                    tAr.stopTimer("Matrix copy");

                    h.i = i;
                    h.j = j;
                    h.e = 0;
                    h.s = 0;
                    h.io.saveBlock(mfBlock, h.e, h.s, h.i, h.j);
                }
            }
            tAr.stopTimer("Matrix I/o");
            tAr.stopTimer("Total");
            printTimeProfile(tAr.getTimings(), tAr.getTimer("Total"));
        }
    }
    FINALIZE();
    
    return EXIT_SUCCESS;
}
