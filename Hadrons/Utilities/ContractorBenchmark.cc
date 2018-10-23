#include <Hadrons/Global.hpp>
#include <Hadrons/DiskVector.hpp>
#ifdef USE_MKL
#include "mkl.h"
#include "mkl_cblas.h"
#endif

using namespace Grid;

#define EIGEN_ROW_MAJOR 1
#define EIGEN_COL_MAJOR 2
#ifndef EIGEN_ORDER 
#define EIGEN_ORDER EIGEN_ROW_MAJOR
#endif
#if (EIGEN_ORDER == EIGEN_ROW_MAJOR)
typedef Eigen::Matrix<ComplexD, -1, -1, Eigen::RowMajor> CMat;
#elif (EIGEN_ORDER == EIGEN_COL_MAJOR)
typedef Eigen::Matrix<ComplexD, -1, -1, Eigen::ColMajor> CMat;
#endif

typedef std::vector<std::vector<CMat>>                      CMatSet;

#pragma omp declare reduction(ComplexPlus: ComplexD: omp_out += omp_in) 

template <typename Function>
inline void trBenchmark(const std::string name, const CMatSet &mat, 
                        const ComplexD ref, Function fn)
{
    double       t, flops, bytes, n = mat[0][0].rows()*mat[0][0].cols();
    unsigned int nMat = mat[0].size();
    ComplexD     buf;

    t = 0.;
    for (unsigned int i = 0; i < nMat; ++i)
    {
        t -= usecond();
        fn(buf, mat[0][i], mat[1][i]);
        t += usecond();
    }
    flops = nMat*(6.*n + 2.*(n - 1.));
    bytes = nMat*(2.*n*sizeof(ComplexD));

    std::cout << std::setw(30) << name << ": diff= "
              << std::setw(12) << std::norm(buf-ref)
              << std::setw(10) << t/1.0e6 << " sec "
              << std::setw(10) << flops/t/1.0e3 << " GFlop/s " 
              << std::setw(10) << bytes/t*1.0e6/1024/1024/1024 << " GB/s "
              << std::endl;
    ::sleep(1);
}

template <typename Function>
inline void mulBenchmark(const std::string name, const CMatSet &mat, 
                         const CMat ref, Function fn)
{
    double       t, flops, bytes;
    double       nr = mat[0][0].rows(), nc = mat[0][0].cols(), n = nr*nc;
    unsigned int nMat = mat[0].size();
    CMat         buf(mat[0][0].rows(), mat[0][0].rows());

    t = 0.;
    for (unsigned int i = 0; i < nMat; ++i)
    {
        t -= usecond();
        fn(buf, mat[0][i], mat[1][i]);
        t += usecond();
    }
    flops = nMat*(n*(6.*nc + 2.*(nc - 1.)));
    bytes = nMat*((2.*n+nr*nr)*sizeof(ComplexD));

    std::cout << std::setw(30) << name << ": diff= "
              << std::setw(12) << (buf-ref).squaredNorm()
              << std::setw(10) << t/1.0e6 << " sec "
              << std::setw(10) << flops/t/1.0e3 << " GFlop/s " 
              << std::setw(10) << bytes/t*1.0e6/1024/1024/1024 << " GB/s "
              << std::endl;
    ::sleep(1);
}

int main(int argc, char *argv[])
{
    // parse command line
    Eigen::Index ni, nj, nMat;

    if (argc != 4)
    {
        std::cerr << "usage: " << argv[0] << " <Ni> <Nj> <#matrices>";
        std::cerr << std::endl;
        
        return EXIT_FAILURE;
    }
    ni   = std::stoi(argv[1]);
    nj   = std::stoi(argv[2]);
    nMat = std::stoi(argv[3]);

    std::cout << "\n*** ALL-TO-ALL MATRIX CONTRACTION BENCHMARK ***\n" << std::endl;
    std::cout << nMat << " " << ni << "x" << nj << " matrices\n" << std::endl;
#ifdef EIGEN_USE_MKL_ALL
    std::cout << "Eigen uses the MKL" << std::endl;
#endif
    std::cout << "Eigen uses ";
#if (EIGEN_ORDER == EIGEN_ROW_MAJOR)
    std::cout << "row-major ordering" << std::endl;
#elif (EIGEN_ORDER == EIGEN_COL_MAJOR)
    std::cout << "column-major ordering" << std::endl;
#endif
    std::cout << "Eigen uses " << Eigen::nbThreads() << " thread(s)" << std::endl;
#ifdef USE_MKL
    std::cout << "MKL   uses " << mkl_get_max_threads() << " thread(s)" << std::endl;
#endif
    std::cout << std::endl;

    CMatSet  mat(2);
    CMat     buf;
    ComplexD ref;

    mat[0].resize(nMat, Eigen::MatrixXcd::Random(ni, nj));
    mat[1].resize(nMat, Eigen::MatrixXcd::Random(nj, ni));

    std::cout << "==== tr(A*B) benchmarks" << std::endl;
    ref = (mat[0].back()*mat[1].back()).trace();
    trBenchmark("Naive loop rows first", mat, ref,
    [](ComplexD &res, const CMat &a, const CMat &b)
    { 
        res = 0.;
        #pragma omp parallel for schedule(static) reduction(ComplexPlus:res)
        for (unsigned int i = 0; i < a.rows(); ++i)
        for (unsigned int j = 0; j < a.cols(); ++j)
        {
            res += a(i, j)*b(j, i);
        }
    });
    trBenchmark("Naive loop cols first", mat, ref,
    [](ComplexD &res, const CMat &a, const CMat &b)
    { 
        res = 0.;
        #pragma omp parallel for schedule(static) reduction(ComplexPlus:res)
        for (unsigned int j = 0; j < a.cols(); ++j)
        for (unsigned int i = 0; i < a.rows(); ++i)
        {
            res += a(i, j)*b(j, i);
        }
    });
    trBenchmark("Eigen tr(A*B)", mat, ref,
    [](ComplexD &res, const CMat &a, const CMat &b)
    { 
        res = (a*b).trace();
    });
    trBenchmark("Eigen global dot", mat, ref,
    [&buf](ComplexD &res, const CMat &a, const CMat &b)
    { 
        buf = b.transpose();
        Eigen::Map<const Eigen::VectorXcd> av(a.data(), a.rows()*a.cols());
        Eigen::Map<const Eigen::VectorXcd> bv(buf.data(), b.rows()*b.cols());

        res = av.conjugate().dot(bv);
    });
    trBenchmark("Eigen row-wise dot", mat, ref,
    [](ComplexD &res, const CMat &a, const CMat &b)
    { 
        res = 0.;
        #pragma omp parallel for schedule(static) reduction(ComplexPlus:res)
        for (unsigned int r = 0; r < a.rows(); ++r)
        {
            res += a.row(r).conjugate().dot(b.col(r));
        }
    });
    trBenchmark("Eigen col-wise dot", mat, ref,
    [](ComplexD &res, const CMat &a, const CMat &b)
    { 
        res = 0.;
        #pragma omp parallel for schedule(static) reduction(ComplexPlus:res)
        for (unsigned int c = 0; c < a.cols(); ++c)
        {
            res += a.col(c).conjugate().dot(b.row(c));
        }
    });
    trBenchmark("Eigen Hadamard", mat, ref,
    [](ComplexD &res, const CMat &a, const CMat &b)
    { 
        res = a.cwiseProduct(b.transpose()).sum();
    });
#ifdef USE_MKL
    trBenchmark("MKL row-wise zdotu", mat, ref,
    [](ComplexD &res, const CMat &a, const CMat &b)
    { 
        res = 0.;
        #pragma omp parallel for schedule(static) reduction(ComplexPlus:res)
        for (unsigned int r = 0; r < a.rows(); ++r)
        {
            ComplexD tmp;

#if (EIGEN_ORDER == EIGEN_ROW_MAJOR)
            cblas_zdotu_sub(a.cols(), a.data() + r*a.cols(), 1, b.data() + r, b.cols(), &tmp);
#elif (EIGEN_ORDER == EIGEN_COL_MAJOR)
            cblas_zdotu_sub(a.cols(), a.data() + r, a.rows(), b.data() + r*b.rows(), 1, &tmp);
#endif
            res += tmp;
        }
    });
    trBenchmark("MKL col-wise zdotu", mat, ref,
    [](ComplexD &res, const CMat &a, const CMat &b)
    { 
        res = 0.;
        #pragma omp parallel for schedule(static) reduction(ComplexPlus:res)
        for (unsigned int c = 0; c < a.cols(); ++c)
        {
            ComplexD tmp;

#if (EIGEN_ORDER == EIGEN_ROW_MAJOR)
            cblas_zdotu_sub(a.rows(), a.data() + c, a.cols(), b.data() + c*b.cols(), 1, &tmp);
#elif (EIGEN_ORDER == EIGEN_COL_MAJOR)
            cblas_zdotu_sub(a.rows(), a.data() + c*a.rows(), 1, b.data() + c, b.rows(), &tmp);
#endif
            res += tmp;
        }
    });
#endif

    std::cout << std::endl;
    std::cout << "==== A*B benchmarks" << std::endl;
    buf = mat[0].back()*mat[1].back();
    mulBenchmark("Naive", mat, buf,
    [](CMat &res, const CMat &a, const CMat &b)
    { 
        unsigned int ni = a.rows(), nj = a.cols();

        #pragma omp parallel for collapse(2)
        for (unsigned int i = 0; i < ni; ++i)
        for (unsigned int k = 0; k < ni; ++k)
        {
            res(i, k) = a(i, 0)*b(0, k);
        }
        #pragma omp parallel for collapse(2)
        for (unsigned int i = 0; i < ni; ++i)
        for (unsigned int k = 0; k < ni; ++k)
        for (unsigned int j = 1; j < nj; ++j)
        {
            res(i, k) += a(i, j)*b(j, k);
        }
    });
    mulBenchmark("Eigen A*B", mat, buf,
    [](CMat &res, const CMat &a, const CMat &b)
    { 
        res = a*b;
    });
#ifdef USE_MKL
    mulBenchmark("MKL A*B", mat, buf,
    [](CMat &res, const CMat &a, const CMat &b)
    {
        const ComplexD one(1., 0.), zero(0., 0.);
#if (EIGEN_ORDER == EIGEN_ROW_MAJOR)
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(),
                    a.cols(), &one, a.data(), a.cols(), b.data(), b.cols(), &zero,
                    res.data(), res.cols());
#elif (EIGEN_ORDER == EIGEN_COL_MAJOR)
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(),
                    a.cols(), &one, a.data(), a.rows(), b.data(), b.rows(), &zero,
                    res.data(), res.rows());
#endif
    });
#endif

    std::cout << std::endl;

    return EXIT_SUCCESS;
}
