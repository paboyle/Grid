#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#ifdef USE_MKL
#include "mkl.h"
#include "mkl_cblas.h"
#endif

using namespace Grid;
using namespace Hadrons;

#ifdef GRID_COMMS_MPI3
#define GET_RANK(rank, nMpi) \
MPI_Comm_size(MPI_COMM_WORLD, &(nMpi));\
MPI_Comm_rank(MPI_COMM_WORLD, &(rank))
#define BARRIER() MPI_Barrier(MPI_COMM_WORLD)
#define INIT() MPI_Init(NULL, NULL)
#define FINALIZE() MPI_Finalize()
#else
#define GET_RANK(rank, nMpi) (nMpi) = 1; (rank) = 0
#define BARRIER()
#define INIT()
#define FINALIZE()
#endif

template <typename Function, typename MatLeft, typename MatRight>
inline void trBenchmark(const std::string name, const MatLeft &left,
                        const MatRight &right, const ComplexD ref, Function fn)
{
    double       t, flops, bytes, n = left[0].rows()*left[0].cols();
    unsigned int nMat = left.size();
    int          nMpi, rank;
    ComplexD     buf;

    t = 0.;
    GET_RANK(rank, nMpi);
    t = -usecond();
    BARRIER();
    for (unsigned int i = rank*nMat/nMpi; i < (rank+1)*nMat/nMpi; ++i)
    {
        fn(buf, left[i], right[i]);      
    }
    BARRIER();
    t += usecond();
    flops = nMat*(6.*n + 2.*(n - 1.));
    bytes = nMat*(2.*n*sizeof(ComplexD));

    if (rank == 0)
    {
        std::cout << std::setw(34) << name << ": diff= "
                  << std::setw(12) << std::norm(buf-ref)
                  << std::setw(10) << t/1.0e6 << " sec "
                  << std::setw(10) << flops/t/1.0e3 << " GFlop/s " 
                  << std::setw(10) << bytes/t*1.0e6/1024/1024/1024 << " GB/s "
                  << std::endl;
    }
    ::sleep(1);
}

template <typename Function, typename MatV, typename Mat>
inline void mulBenchmark(const std::string name, const MatV &left,
                         const MatV &right, const Mat &ref, Function fn)
{
    double       t, flops, bytes;
    double       nr = left[0].rows(), nc = left[0].cols(), n = nr*nc;
    unsigned int nMat = left.size();
    int          nMpi, rank;
    Mat          buf(left[0].rows(), left[0].rows());

    t = 0.;
    GET_RANK(rank, nMpi);
    t = -usecond();
    BARRIER();
    for (unsigned int i = rank*nMat/nMpi; i < (rank+1)*nMat/nMpi; ++i)
    {
        fn(buf, left[i], right[i]);
    }
    BARRIER();
    t += usecond();
    flops = nMat*(nr*nr*(6.*nc + 2.*(nc - 1.)));
    bytes = nMat*(2*nc*nr*sizeof(ComplexD));

    if (rank == 0)
    {
        std::cout << std::setw(34) << name << ": diff= "
                  << std::setw(12) << (buf-ref).squaredNorm()
                  << std::setw(10) << t/1.0e6 << " sec "
                  << std::setw(10) << flops/t/1.0e3 << " GFlop/s " 
                  << std::setw(10) << bytes/t*1.0e6/1024/1024/1024 << " GB/s "
                  << std::endl;
    }
    ::sleep(1);
}

#ifdef USE_MKL
template <typename MatLeft, typename MatRight>
static inline void zdotuRow(ComplexD &res, const unsigned int aRow,
                            const MatLeft &a, const MatRight &b)
{
    const ComplexD *aPt, *bPt;
    unsigned int   aInc, bInc;

    if (MatLeft::Options == Eigen::RowMajor)
    {
        aPt  = a.data() + aRow*a.cols();
        aInc = 1;
    }
    else if (MatLeft::Options == Eigen::ColMajor)
    {
        aPt  = a.data() + aRow;
        aInc = a.rows();
    }
    if (MatRight::Options == Eigen::RowMajor)
    {
        bPt  = b.data() + aRow;
        bInc = b.cols();
    }
    else if (MatRight::Options == Eigen::ColMajor)
    {
        bPt  = b.data() + aRow*b.rows();
        bInc = 1;
    }
    cblas_zdotu_sub(a.cols(), aPt, aInc, bPt, bInc, &res);
}

template <typename MatLeft, typename MatRight>
static inline void zdotuCol(ComplexD &res, const unsigned int aCol,
                            const MatLeft &a, const MatRight &b)
{
    const ComplexD *aPt, *bPt;
    unsigned int   aInc, bInc;

    if (MatLeft::Options == Eigen::RowMajor)
    {
        aPt  = a.data() + aCol;
        aInc = a.cols();
    }
    else if (MatLeft::Options == Eigen::ColMajor)
    {
        aPt  = a.data() + aCol*a.rows();
        aInc = 1;
    }
    if (MatRight::Options == Eigen::RowMajor)
    {
        bPt  = b.data() + aCol*b.cols();
        bInc = 1;
    }
    else if (MatRight::Options == Eigen::ColMajor)
    {
        bPt  = b.data() + aCol;
        bInc = b.rows();
    }
    cblas_zdotu_sub(a.rows(), aPt, aInc, bPt, bInc, &res);
}
#endif

template <typename MatLeft, typename MatRight>
void fullTrBenchmark(const unsigned int ni, const unsigned int nj, const unsigned int nMat)
{
    std::vector<MatLeft>  left;
    std::vector<MatRight> right;
    MatRight              buf;
    ComplexD              ref;
    int                   rank, nMpi;

    left.resize(nMat, MatLeft::Random(ni, nj));
    right.resize(nMat, MatRight::Random(nj, ni));
    GET_RANK(rank, nMpi);
    if (rank == 0)
    {
        std::cout << "==== tr(A*B) benchmarks" << std::endl;
        std::cout << "A matrices use ";
        if (MatLeft::Options == Eigen::RowMajor)
        {
            std::cout << "row-major ordering" << std::endl;
        }
        else if (MatLeft::Options == Eigen::ColMajor)
        {
            std::cout << "col-major ordering" << std::endl;
        }
        std::cout << "B matrices use ";
        if (MatRight::Options == Eigen::RowMajor)
        {
            std::cout << "row-major ordering" << std::endl;
        }
        else if (MatRight::Options == Eigen::ColMajor)
        {
            std::cout << "col-major ordering" << std::endl;
        }
        std::cout << std::endl;
    }
    BARRIER();
    ref = (left.back()*right.back()).trace();
    trBenchmark("Hadrons A2AContraction::accTrMul", left, right, ref,
    [](ComplexD &res, const MatLeft &a, const MatRight &b)
    { 
        res = 0.;
        A2AContraction::accTrMul(res, a, b);
    });
    trBenchmark("Naive loop rows first", left, right, ref,
    [](ComplexD &res, const MatLeft &a, const MatRight &b)
    { 
        auto nr = a.rows(), nc = a.cols();
        
        res = 0.;
        parallel_for (unsigned int i = 0; i < nr; ++i)
        {
            ComplexD tmp = 0.;

            for (unsigned int j = 0; j < nc; ++j)
            {
                tmp += a(i, j)*b(j, i);
            }
            parallel_critical
            {
                res += tmp;
            }
        }
    });
    trBenchmark("Naive loop cols first", left, right, ref,
    [](ComplexD &res, const MatLeft &a, const MatRight &b)
    {
        auto nr = a.rows(), nc = a.cols();
        
        res = 0.;
        parallel_for (unsigned int j = 0; j < nc; ++j)
        {
            ComplexD tmp = 0.;

            for (unsigned int i = 0; i < nr; ++i)
            {
                tmp += a(i, j)*b(j, i);
            }        
            parallel_critical
            {
                res += tmp;
            }
        }
    });
    trBenchmark("Eigen tr(A*B)", left, right, ref,
    [](ComplexD &res, const MatLeft &a, const MatRight &b)
    { 
        res = (a*b).trace();
    });
    trBenchmark("Eigen row-wise dot", left, right, ref,
    [](ComplexD &res, const MatLeft &a, const MatRight &b)
    {
        res = 0.;
        parallel_for (unsigned int r = 0; r < a.rows(); ++r)
        {
            ComplexD tmp;

            tmp = a.row(r).conjugate().dot(b.col(r));
            parallel_critical
            {
                res += tmp;
            }
        }
    });
    trBenchmark("Eigen col-wise dot", left, right, ref,
    [](ComplexD &res, const MatLeft &a, const MatRight &b)
    {
        res = 0.;
        parallel_for (unsigned int c = 0; c < a.cols(); ++c)
        {
            ComplexD tmp;

            tmp = a.col(c).conjugate().dot(b.row(c));
            parallel_critical
            {
                res += tmp;
            }
        }
    });
    trBenchmark("Eigen Hadamard", left, right, ref,
    [](ComplexD &res, const MatLeft &a, const MatRight &b)
    { 
        res = a.cwiseProduct(b.transpose()).sum();
    });
#ifdef USE_MKL
    trBenchmark("MKL row-wise zdotu", left, right, ref,
    [](ComplexD &res, const MatLeft &a, const MatRight &b)
    {
        res = 0.;
        parallel_for (unsigned int r = 0; r < a.rows(); ++r)
        {
            ComplexD tmp;

            zdotuRow(tmp, r, a, b);
            parallel_critical
            {
                res += tmp;
            }
        }
    });
    trBenchmark("MKL col-wise zdotu", left, right, ref,
    [](ComplexD &res, const MatLeft &a, const MatRight &b)
    {
        res = 0.;
        parallel_for (unsigned int c = 0; c < a.cols(); ++c)
        {
            ComplexD tmp;

            zdotuCol(tmp, c, a, b);
            parallel_critical
            {
                res += tmp;
            }
        }
    });
#endif
    BARRIER();
    if (rank == 0)
    {
        std::cout << std::endl;
    }
}

template <typename Mat>
void fullMulBenchmark(const unsigned int ni, const unsigned int nj, const unsigned int nMat)
{
    std::vector<Mat> left, right;
    Mat              ref;
    int              rank, nMpi;

    left.resize(nMat, Mat::Random(ni, nj));
    right.resize(nMat, Mat::Random(nj, ni));
    GET_RANK(rank, nMpi);
    if (rank == 0)
    {
        std::cout << "==== A*B benchmarks" << std::endl;
        std::cout << "all matrices use ";
        if (Mat::Options == Eigen::RowMajor)
        {
            std::cout << "row-major ordering" << std::endl;
        }
        else if (Mat::Options == Eigen::ColMajor)
        {
            std::cout << "col-major ordering" << std::endl;
        }
        std::cout << std::endl;
    }
    BARRIER();
    ref = left.back()*right.back();
    mulBenchmark("Hadrons A2AContraction::mul", left, right, ref,
    [](Mat &res, const Mat &a, const Mat &b)
    { 
        A2AContraction::mul(res, a, b);
    });
    mulBenchmark("Eigen A*B", left, right, ref,
    [](Mat &res, const Mat &a, const Mat &b)
    { 
        res = a*b;
    });
#ifdef USE_MKL
    mulBenchmark("MKL A*B", left, right, ref,
    [](Mat &res, const Mat &a, const Mat &b)
    {
        const ComplexD one(1., 0.), zero(0., 0.);
        if (Mat::Options == Eigen::RowMajor)
        {
            cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(),
                        a.cols(), &one, a.data(), a.cols(), b.data(), b.cols(), &zero,
                        res.data(), res.cols());
        }
        else if (Mat::Options == Eigen::ColMajor)
        {
            cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(),
                        a.cols(), &one, a.data(), a.rows(), b.data(), b.rows(), &zero,
                        res.data(), res.rows());
        }
    });
#endif
    BARRIER();
    if (rank == 0)
    {
        std::cout << std::endl;
    }
}

int main(int argc, char *argv[])
{
    // parse command line
    Eigen::Index ni, nj, nMat;
    int          nMpi, rank;

    if (argc != 4)
    {
        std::cerr << "usage: " << argv[0] << " <Ni> <Nj> <#matrices>";
        std::cerr << std::endl;
        
        return EXIT_FAILURE;
    }
    ni   = std::stoi(argv[1]);
    nj   = std::stoi(argv[2]);
    nMat = std::stoi(argv[3]);

    INIT();
    GET_RANK(rank, nMpi);
    if (rank == 0)
    {
        std::cout << "\n*** ALL-TO-ALL MATRIX CONTRACTION BENCHMARK ***\n" << std::endl;
        std::cout << nMat << " couples of " << ni << "x" << nj << " matrices\n" << std::endl;

        std::cout << nMpi << " MPI processes" << std::endl;
#ifdef GRID_OMP
        #pragma omp parallel
        {
            #pragma omp single
            std::cout << omp_get_num_threads() << " threads\n" << std::endl; 
        }
#else
        std::cout << "Single-threaded\n" << std::endl; 
#endif

#ifdef EIGEN_USE_MKL_ALL
        std::cout << "Eigen uses the MKL" << std::endl;
#endif
        std::cout << "Eigen uses " << Eigen::nbThreads() << " threads" << std::endl;
#ifdef USE_MKL
        std::cout << "MKL   uses " << mkl_get_max_threads() << " threads" << std::endl;
#endif
        std::cout << std::endl;
    }

    fullTrBenchmark<A2AMatrix<ComplexD>, A2AMatrix<ComplexD>>(ni, nj, nMat);
    fullTrBenchmark<A2AMatrix<ComplexD>, A2AMatrixTr<ComplexD>>(ni, nj, nMat);
    fullTrBenchmark<A2AMatrixTr<ComplexD>, A2AMatrix<ComplexD>>(ni, nj, nMat);
    fullTrBenchmark<A2AMatrixTr<ComplexD>, A2AMatrixTr<ComplexD>>(ni, nj, nMat);
    fullMulBenchmark<A2AMatrix<ComplexD>>(ni, nj, nMat);
    fullMulBenchmark<A2AMatrixTr<ComplexD>>(ni, nj, nMat);
    FINALIZE();

    return EXIT_SUCCESS;
}
