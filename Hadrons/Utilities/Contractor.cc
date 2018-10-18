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
                                        std::vector<std::string>, timeRange);
    };
}

struct ContractorPar
{
    Contractor::GlobalPar                  global;
    std::vector<Contractor::A2AMatrixPar>  a2aMatrix;
    std::vector<Contractor::ProductPar>    product;
};

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

        std::cout << "-- Loading '" << p.file << "'..." << std::endl;

        A2AMatrixIo<HADRONS_A2AM_IO_TYPE> a2aIo(p.file, p.dataset, par.global.nt);

        a2aIo.load(a2aMat.at(p.name), &t);
        std::cout << "Read " << a2aIo.getSize() << " bytes in " << t << " usec, " << a2aIo.getSize()/t*1.0e6/1024/1024 << " MB/s" << std::endl;
    }

    // contract
    EigenDiskVector<ComplexD>::Matrix buf;

    for (auto &p: par.product)
    {
        std::vector<std::string>            term = strToVec<std::string>(p.terms);
        std::vector<std::set<unsigned int>> times;

        if (term.size() != p.timeRange.size())
        {
            HADRONS_ERROR(Size, "number of terms (" + std::to_string(term.size()) 
                          + ") different from number of time ranges (" 
                          + std::to_string(p.timeRange.size()) + ")");
        }
        for (auto &s: p.timeRange)
        {
            times.push_back(parseTimeRange(s, par.global.nt));
        }
        
        // test: just 2-pt function for now
        if (term.size() != 2)
        {
            HADRONS_ERROR(Implementation, "only 2-pt function implemented");
        }

        std::vector<ComplexD> corr(par.global.nt);
        double                tTrace = 0., flops, bytes;

#define TIME_MOD(t) (((t) + par.global.nt) % par.global.nt)
        for (unsigned int t = 0; t < par.global.nt; ++t)
        {
            corr[t] = ComplexD(0., 0.);
        }
        for (auto &trans: times[1])
        for (auto &t0: times[0])
        {
            unsigned int tt0 = TIME_MOD(t0 + trans);
            EigenDiskVector<ComplexD>::Matrix m0 = a2aMat.at(term[0])[tt0];
            EigenDiskVector<ComplexD>::Matrix m1 = a2aMat.at(term[1])[trans];

            std::cout << "tr(" << term[0] << "[" << tt0 << "]*"
                      << term[1] << "[" << trans << "])";
            tTrace  = -usecond();
            corr[t0] += (m0*m1).trace();
            tTrace += usecond();
            flops = 6.*m0.rows()*m0.cols() + 2.*(m0.rows()*m0.cols() - 1.);
            bytes = 2.*m0.rows()*m0.cols()*sizeof(ComplexD);
            std::cout << " -- Perf " << flops/tTrace/1.0e3 << " GFlop/s " 
                      << bytes/tTrace*1.0e6/1024/1024/1024 << " GB/s "<< std::endl;
        }
        for (unsigned int t = 0; t < par.global.nt; ++t)
        {
            std::cout << t << " " << corr[t] << std::endl;
        }
#undef TIME_MOD
   }

    
    return EXIT_SUCCESS;
}
