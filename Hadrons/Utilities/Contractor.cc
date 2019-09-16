/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Utilities/Contractor.cc

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>

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

#include <chrono>
#include <ctime>
#include <Hadrons/Utilities/Contractor.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DiskVector.hpp>
#include <Hadrons/TimerArray.hpp>

using namespace Grid;
using namespace Hadrons;

// Separator to be used between contraction terms only (underscores elsewhere)
std::string Separator{ "_" };

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

  makeTimeSeq(timeSeq, times, current, static_cast<unsigned int>(times.size()));
}

void saveCorrelator(const Contractor::CorrelatorResult &result, const std::string dir,
                    const unsigned int dt, const unsigned int traj)
{
    std::string              fileStem = "", filename;
    std::vector<std::string> terms = strToVec<std::string>(result.contraction.terms);

    for (unsigned int i = 0; i < terms.size() - 1; i++)
    {
        fileStem += terms[i] + Separator + std::to_string(result.times[i]) + Separator;
    }
    fileStem += terms.back();
    if (!result.contraction.translationAverage)
    {
        fileStem += Separator + "dt" + Separator + std::to_string(dt);
    }
    filename = dir + "/" + RESULT_FILE_NAME(fileStem, traj);
    std::cout << "Saving correlator to '" << filename << "'" << std::endl;
    makeFileDir(dir);
    ResultWriter writer(filename);
    write(writer, fileStem, result);
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
    // parse command line
    std::string parFilename;
    bool        bOnlyWriteUsedA2AMatrices{ false };
    int         ArgCount{ 0 };
    bool        bCmdLineError{ false };
    for( int i = 1; i < argc; i++ ) {
        if( argv[i][0] == '-' ) {
            // Switches
            bool bSwitchOK = false;
            switch( argv[i][1] ) {
                case 'a':
                    if( argv[i][2] == 0 ) {
                        bOnlyWriteUsedA2AMatrices = true;
                        bSwitchOK = true;
                        std::cout << "Only A2AMatrices used in each contraction will be written" << std::endl;
                    }
                    break;
                case 's':
                    if( argv[i][2] )
                        Separator = &argv[i][2];
                    else
                        Separator = ".";
                    bSwitchOK = true;
                    std::cout << "Using \"" << Separator << "\" as name separator" << std::endl;
                    break;
            }
            if( !bSwitchOK ) {
                std::cerr << "Urecognised switch \"" << argv[i] << "\"" << std::endl;
                bCmdLineError = true;
            }
        } else {
            // Arguments
            switch( ++ArgCount ) {
                case 1:
                    parFilename = argv[i];
                    break;
                default:
                    std::cerr << "Unused argument \"" << argv[i] << "\"" << std::endl;
                    break;
            }
        }
    }

    if (ArgCount != 1 or bCmdLineError)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file>"
                     "\n\t-a\tSimple Correlators (only describe A2AMatrices used for contraction)"
                     "\n\t-s[sep]\tSeparator \"sep\" used between name components."
                     "\n\t\tDefaults to \"_\", or \".\" if -s specified without sep"
                  << std::endl;
        
        return EXIT_FAILURE;
    }

    // Log what file we're processing and when we started
    const std::chrono::system_clock::time_point start{ std::chrono::system_clock::now() };
    std::time_t now = std::chrono::system_clock::to_time_t( start );
    std::cout << "Start " << parFilename << " " << std::ctime( &now );

    // parse parameter file
    Contractor::ContractorPar par;
    XmlReader     reader(parFilename);

    read(reader, "global",    par.global);
    read(reader, "a2aMatrix", par.a2aMatrix);
    read(reader, "product",   par.product);
    const unsigned int nMat  { static_cast<unsigned int>(par.a2aMatrix.size()) };
    const unsigned int nCont { static_cast<unsigned int>(par.product.size()) };

    // create diskvectors
    std::map<std::string, EigenDiskVector<ComplexD>> a2aMat;

    for (auto &p: par.a2aMatrix)
    {
        std::string dirName = par.global.diskVectorDir + "/" + p.name;

        a2aMat.emplace(p.name, EigenDiskVector<ComplexD>(dirName, par.global.nt, p.cacheSize));
    }

    // trajectory loop
    for (unsigned int traj = par.global.trajCounter.start; 
         traj < par.global.trajCounter.end; traj += par.global.trajCounter.step)
    {
        std::cout << ":::::::: Trajectory " << traj << std::endl;

        // load data
        int iSeq = 0;
        for (auto &p: par.a2aMatrix)
        {
            std::string filename = p.file;
            double      t;

            tokenReplace(filename, "traj", traj);
            std::cout << "======== Loading '" << filename << "'"
                      << "\nA2AMatrix " << ++iSeq << " of " << nMat << " = " << p.name << std::endl;

            A2AMatrixIo<HADRONS_A2AM_IO_TYPE> a2aIo(filename, p.dataset, par.global.nt);

            a2aIo.load(a2aMat.at(p.name), &t);
            std::cout << "Read " << a2aIo.getSize() << " bytes in " << t/1.0e6 
                    << " sec, " << a2aIo.getSize()/t*1.0e6/1024/1024 << " MB/s" << std::endl;
        }

        // contract
        EigenDiskVector<ComplexD>::Matrix buf;

        iSeq = 0;
        for (auto &p: par.product)
        {
            std::vector<std::string>               term = strToVec<std::string>(p.terms);
            std::vector<std::set<unsigned int>>    times;
            std::vector<std::vector<unsigned int>> timeSeq;
            std::set<unsigned int>                 translations;
            std::vector<A2AMatrixTr<ComplexD>>     lastTerm(par.global.nt);
            A2AMatrix<ComplexD>                    prod, buf, tmp;
            TimerArray                             tAr;
            double                                 fusec, busec, flops, bytes;
            Contractor::CorrelatorResult           result;             

            tAr.startTimer("Total");
            std::cout << "======== Contraction " << ++iSeq << " of " << nCont << " tr(";
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
            for (auto &m: par.a2aMatrix)
            {
                // For simple correlators, only include A2AMatrix info for correlators in this contraction
                if ( ( !bOnlyWriteUsedA2AMatrices or std::find( term.begin(), term.end(), m.name ) != term.end() )
                  and std::find(result.a2aMatrix.begin(), result.a2aMatrix.end(), m) == result.a2aMatrix.end())
                {
                    result.a2aMatrix.push_back(m);
                    tokenReplace(result.a2aMatrix.back().file, "traj", traj);
                }
            }
            result.contraction = p;
            result.correlator.resize(par.global.nt, 0.);

            translations = parseTimeRange(p.translations, par.global.nt);
            makeTimeSeq(timeSeq, times);
            std::cout << timeSeq.size()*translations.size()*(term.size() - 2) << " A*B, "
                    << timeSeq.size()*translations.size()*par.global.nt << " tr(A*B)"
                    << std::endl;

            std::cout << "* Caching transposed last term" << std::endl;
            for (unsigned int t = 0; t < par.global.nt; ++t)
            {
                tAr.startTimer("Disk vector overhead");
                const A2AMatrix<ComplexD> &ref = a2aMat.at(term.back())[t];
                tAr.stopTimer("Disk vector overhead");

                tAr.startTimer("Transpose caching");
                lastTerm[t].resize(ref.rows(), ref.cols());
                thread_for( j,ref.cols(),{
                    for (unsigned int i = 0; i < ref.rows(); ++i)
                    {
                        lastTerm[t](i, j) = ref(i, j);
                    }
                });
                tAr.stopTimer("Transpose caching");
            }
            bytes = par.global.nt*lastTerm[0].rows()*lastTerm[0].cols()*sizeof(ComplexD);
            std::cout << Sec(tAr.getDTimer("Transpose caching")) << " " 
                      << Bytes(bytes, tAr.getDTimer("Transpose caching")) << std::endl;
            for (unsigned int i = 0; i < timeSeq.size(); ++i)
            {
                unsigned int dti = 0;
                auto         &t = timeSeq[i];

                result.times = t;
                for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                {
                    result.correlator[tLast] = 0.;
                }
                for (auto &dt: translations)
                {
                    std::cout << "* Step " << i*translations.size() + dti + 1
                            << "/" << timeSeq.size()*translations.size()
                            << " -- positions= " << t << ", dt= " << dt << std::endl;
                    if (term.size() > 2)
                    {
                        std::cout << std::setw(8) << "products";
                    }
                    flops  = 0.;
                    bytes  = 0.;
                    fusec  = tAr.getDTimer("A*B algebra");
                    busec  = tAr.getDTimer("A*B total");
                    tAr.startTimer("Linear algebra");
                    tAr.startTimer("Disk vector overhead");
                    prod = a2aMat.at(term[0])[TIME_MOD(t[0] + dt)];
                    tAr.stopTimer("Disk vector overhead");
                    for (unsigned int j = 1; j < term.size() - 1; ++j)
                    {
                        tAr.startTimer("Disk vector overhead");
                        const A2AMatrix<ComplexD> &ref = a2aMat.at(term[j])[TIME_MOD(t[j] + dt)];
                        tAr.stopTimer("Disk vector overhead");
                        
                        tAr.startTimer("A*B total");
                        tAr.startTimer("A*B algebra");
                        A2AContraction::mul(tmp, prod, ref);
                        tAr.stopTimer("A*B algebra");
                        flops += A2AContraction::mulFlops(prod, ref);
                        prod   = tmp;
                        tAr.stopTimer("A*B total");
                        bytes += 3.*tmp.rows()*tmp.cols()*sizeof(ComplexD);
                    }
                    if (term.size() > 2)
                    {
                        std::cout << Sec(tAr.getDTimer("A*B total") - busec) << " "
                                << Flops(flops, tAr.getDTimer("A*B algebra") - fusec) << " " 
                                << Bytes(bytes, tAr.getDTimer("A*B total") - busec) << std::endl;
                    }
                    std::cout << std::setw(8) << "traces";
                    flops  = 0.;
                    bytes  = 0.;
                    fusec  = tAr.getDTimer("tr(A*B)");
                    busec  = tAr.getDTimer("tr(A*B)");
                    for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                    {
                        tAr.startTimer("tr(A*B)");
                        A2AContraction::accTrMul(result.correlator[TIME_MOD(tLast - dt)], prod, lastTerm[tLast]);
                        tAr.stopTimer("tr(A*B)");
                        flops += A2AContraction::accTrMulFlops(prod, lastTerm[tLast]);
                        bytes += 2.*prod.rows()*prod.cols()*sizeof(ComplexD);
                    }
                    tAr.stopTimer("Linear algebra");
                    std::cout << Sec(tAr.getDTimer("tr(A*B)") - busec) << " "
                            << Flops(flops, tAr.getDTimer("tr(A*B)") - fusec) << " " 
                            << Bytes(bytes, tAr.getDTimer("tr(A*B)") - busec) << std::endl;
                    if (!p.translationAverage)
                    {
                        saveCorrelator(result, par.global.output, dt, traj);
                        for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                        {
                            result.correlator[tLast] = 0.;
                        }
                    }
                    dti++;
                }
                if (p.translationAverage)
                {
                    for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                    {
                        result.correlator[tLast] /= translations.size();
                    }
                    saveCorrelator(result, par.global.output, 0, traj);
                }
            }
            tAr.stopTimer("Total");
            printTimeProfile(tAr.getTimings(), tAr.getTimer("Total"));
        }
    }

    // Mention that we're finished, what the time is and how long it took
    const std::chrono::system_clock::time_point stop{ std::chrono::system_clock::now() };
    now = std::chrono::system_clock::to_time_t( stop );
    const std::chrono::duration<double> duration_seconds = stop - start;
    const double hours{ ( duration_seconds.count() + 0.5 ) / 3600 };
    std::cout << "Stop " << parFilename << " " << std::ctime( &now )
              << "Total duration " << std::fixed << std::setprecision(1) << hours << " hours." << std::endl;
    return EXIT_SUCCESS;
}
