#include <Hadrons/Global.hpp>
#include <Hadrons/DiskVector.hpp>

using namespace Grid;

#ifndef EIGEN_ORDER 
#define EIGEN_ORDER RowMajor
#endif

typedef Eigen::Matrix<ComplexD, -1, -1, Eigen::EIGEN_ORDER> CMat;
typedef std::vector<std::vector<CMat>>                      CMatSet;

template <typename TwoMatFn>
inline void trBenchmark(const std::string name, const CMatSet &mat, TwoMatFn fn)
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

    std::cout << std::setw(30) << name << ": result= "
              << std::setw(10) << buf
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

    std::cout << "==== generating random matrices" << std::endl;    
    CMatSet mat(2);
    CMat    buf;

    mat[0].resize(nMat, Eigen::MatrixXcd::Random(ni, nj));
    mat[1].resize(nMat, Eigen::MatrixXcd::Random(nj, ni));

    std::cout << "==== tr(A*B) benchmark" << std::endl;
    trBenchmark("Naive loop rows first", mat,
    [](ComplexD &res, const CMat &a, const CMat &b)
    { 
        res = 0.;
        #pragma omp parallel for schedule(static) reduction(+:res)
        for (unsigned int i = 0; i < a.rows(); ++i)
        for (unsigned int j = 0; j < a.cols(); ++j)
        {
            res += a(i, j)*b(j, i);
        }
    });
    trBenchmark("Naive loop cols first", mat,
    [](ComplexD &res, const CMat &a, const CMat &b)
    { 
        res = 0.;
        for (unsigned int j = 0; j < a.cols(); ++j)
        for (unsigned int i = 0; i < a.rows(); ++i)
        {
            res += a(i, j)*b(j, i);
        }
    });
    trBenchmark("Eigen tr(A*B)", mat,
    [](ComplexD &res, const CMat &a, const CMat &b)
    { 
        res = (a*b).trace();
    });
    trBenchmark("Eigen global dot", mat,
    [&buf](ComplexD &res, const CMat &a, const CMat &b)
    { 
        buf = b.transpose();
        Eigen::Map<const Eigen::VectorXcd> av(a.data(), a.rows()*a.cols());
        Eigen::Map<const Eigen::VectorXcd> bv(buf.data(), b.rows()*b.cols());

        res = av.conjugate().dot(bv);
    });
    trBenchmark("Eigen row-wise dot", mat,
    [](ComplexD &res, const CMat &a, const CMat &b)
    { 
        res = 0.;
        #pragma omp parallel for schedule(static) reduction(+:res)
        for (unsigned int r = 0; r < a.rows(); ++r)
        {
            res += a.row(r).conjugate().dot(b.col(r));
        }
    });
    trBenchmark("Eigen col-wise dot", mat,
    [](ComplexD &res, const CMat &a, const CMat &b)
    { 
        res = 0.;
        #pragma omp parallel for schedule(static) reduction(+:res)
        for (unsigned int r = 0; r < a.cols(); ++r)
        {
            res += a.col(r).conjugate().dot(b.row(r));
        }
    });
    trBenchmark("Eigen Hadamard", mat,
    [](ComplexD &res, const CMat &a, const CMat &b)
    { 
        res = a.cwiseProduct(b.transpose()).sum();
    });

    return EXIT_SUCCESS;
}
