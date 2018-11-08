/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Utilities/Contractor.cc

Copyright (C) 2015-2018


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
#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DiskVector.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

#define TIME_MOD(t) (((t) + par.global.nt) % par.global.nt)

namespace Contractor
{
    class GlobalPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(GlobalPar,
                                        unsigned int, nt,
                                        std::string, diskVectorDir,
                                        std::string, output);
    };

    class A2AMatrixPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMatrixPar,
                                        std::string, file,
                                        std::string, dataset,
                                        unsigned int, cacheSize,
                                        std::string, name);
    };

    class ProductPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(ProductPar,
                                        std::string, terms,
                                        std::vector<std::string>, times,
                                        std::string, translations);
    };
}

struct ContractorPar
{
    Contractor::GlobalPar                  global;
    std::vector<Contractor::A2AMatrixPar>  a2aMatrix;
    std::vector<Contractor::ProductPar>    product;
};

void makeTimeSeq(std::vector<std::vector<unsigned int>> &timeSeq, 
                 const std::vector<std::set<unsigned int>> &times,
                 std::vector<unsigned int> &current,
                 const unsigned int depth)
{
    if (depth > 0)
    {
        for (auto t: times[times.size() - depth])
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


std::set<unsigned int> parseTimeRange(const std::string str, const unsigned int nt)
{
    std::regex               rex("([0-9]+)|(([0-9]+)\\.\\.([0-9]+))");
    std::smatch              sm;
    std::vector<std::string> rstr = strToVec<std::string>(str);
    std::set<unsigned int>   tSet;
    

    for (auto &s: rstr)
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
    ContractorPar par;
    unsigned int  nMat, nCont;
    XmlReader     reader(parFilename);

    read(reader, "global",    par.global);
    read(reader, "a2aMatrix", par.a2aMatrix);
    read(reader, "product",   par.product);
    nMat  = par.a2aMatrix.size();
    nCont = par.product.size();

    // create diskvectors
    std::map<std::string, EigenDiskVector<ComplexD>> a2aMat;
    unsigned int                                     cacheSize;

    for (auto &p: par.a2aMatrix)
    {
        std::string dirName = par.global.diskVectorDir + "/" + p.name;

        a2aMat.emplace(p.name, EigenDiskVector<ComplexD>(dirName, par.global.nt, p.cacheSize));
    }

    // load data
    for (unsigned int i = 0; i < a2aMat.size(); ++i)
    {
        auto   &p = par.a2aMatrix[i];
        double t, size;

        std::cout << "======== Loading '" << p.file << "'" << std::endl;

        A2AMatrixIo<HADRONS_A2AM_IO_TYPE> a2aIo(p.file, p.dataset, par.global.nt);

        a2aIo.load(a2aMat.at(p.name), &t);
        std::cout << "Read " << a2aIo.getSize() << " bytes in " << t << " usec, " << a2aIo.getSize()/t*1.0e6/1024/1024 << " MB/s" << std::endl;
    }

    // contract
    EigenDiskVector<ComplexD>::Matrix buf;

    for (auto &p: par.product)
    {
        std::vector<std::string>               term = strToVec<std::string>(p.terms);
        std::vector<std::set<unsigned int>>    times;
        std::vector<std::vector<unsigned int>> timeSeq;
        std::set<unsigned int>                 translations;
        std::vector<ComplexD>                  corr(par.global.nt);
        std::vector<A2AMatrixTr<ComplexD>>     lastTerm(par.global.nt);
        A2AMatrix<ComplexD>                    prod, buf, tmp;
        double                                 fusec, busec, flops, bytes;

        std::cout << "======== Product tr(";
        for (unsigned int g = 0; g < term.size(); ++g)
        {
            std::cout << term[g] << ((g == term.size() - 1) ? ')' : '*');
        }
        std::cout << std::endl;
        if (term.size() != p.times.size() + 1)
        {
            HADRONS_ERROR(Size, "number of terms (" + std::to_string(term.size()) 
                          + ") different from number of times (" 
                          + std::to_string(p.times.size() + 1) + ")");
        }
        for (auto &s: p.times)
        {
            times.push_back(parseTimeRange(s, par.global.nt));
        }
        translations = parseTimeRange(p.translations, par.global.nt);
        makeTimeSeq(timeSeq, times);
        std::cout << timeSeq.size()*translations.size()*(term.size() - 2) << " A*B, "
                  << timeSeq.size()*translations.size()*par.global.nt << " tr(A*B)"
                  << std::endl;

        std::cout << "-- caching transposed last term" << std::endl;
        for (unsigned int t = 0; t < par.global.nt; ++t)
        {
            const A2AMatrix<ComplexD> &ref = a2aMat.at(term.back())[t];

            lastTerm[t] = ref;
        }
        for (auto &t: timeSeq)
        {
            for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
            {
                corr[tLast] = 0.;
            }
            for (auto &dt: translations)
            {
                std::cout << "-- position " << t << ", translation " << dt << std::endl;
                if (term.size() > 2)
                {
                    std::cout << "*" << std::setw(12) << "products";
                }
                flops  = 0.;
                bytes  = 0.;
                fusec  = 0.;
                busec  = 0.;
                prod = a2aMat.at(term[0])[TIME_MOD(t[0] + dt)];
                for (unsigned int i = 1; i < term.size() - 1; ++i)
                {
                    const A2AMatrix<ComplexD> &ref = a2aMat.at(term[i])[TIME_MOD(t[i] + dt)];
                    fusec -= usecond();
                    busec -= usecond();
                    A2AContraction::mul(tmp, prod, ref);
                    fusec += usecond();
                    flops += A2AContraction::mulFlops(prod, ref);
                    prod   = tmp;
                    busec += usecond();
                    bytes += 3.*tmp.rows()*tmp.cols()*sizeof(ComplexD);
                }
                if (term.size() > 2)
                {
                    std::cout << Flops(flops, fusec) << " " << Bytes(bytes, busec) << std::endl;
                }
                std::cout << "*" << std::setw(12) << "traces";
                flops  = 0.;
                bytes  = 0.;
                fusec  = 0.;
                busec  = 0.;
                for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                {
                    fusec -= usecond();
                    busec -= usecond();
                    A2AContraction::accTrMul(corr[TIME_MOD(tLast - dt)], prod, lastTerm[tLast]);
                    fusec += usecond();
                    busec += usecond();
                    flops += A2AContraction::accTrMulFlops(prod, lastTerm[tLast]);
                    bytes += 2.*prod.rows()*prod.cols()*sizeof(ComplexD);
                }
                std::cout << Flops(flops, fusec) << " " << Bytes(bytes, busec) << std::endl;
            }
            for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
            {
                std::cout << tLast << " " << corr[tLast] << std::endl;
            }
        }
    }
    
    return EXIT_SUCCESS;
}
