/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/EigenPack.hpp

Copyright (C) 2015-2018

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
#ifndef Hadrons_EigenPack_hpp_
#define Hadrons_EigenPack_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/algorithms/iterative/Deflation.h>
#include <Grid/algorithms/iterative/LocalCoherenceLanczos.h>

BEGIN_HADRONS_NAMESPACE

// Lanczos type
#ifndef HADRONS_DEFAULT_LANCZOS_NBASIS
#define HADRONS_DEFAULT_LANCZOS_NBASIS 60
#endif

template <typename F>
class EigenPack
{
public:
    typedef F Field;
    struct PackRecord: Serializable
    {
        GRID_SERIALIZABLE_CLASS_MEMBERS(PackRecord,
                                        std::string, operatorPar,
                                        std::string, solverPar);
    };
    struct VecRecord: Serializable
    {
        GRID_SERIALIZABLE_CLASS_MEMBERS(VecRecord,
                                        unsigned int, index,
                                        double,       eval);
        VecRecord(void): index(0), eval(0.) {}
    };
public:
    std::vector<RealD> eval;
    std::vector<F>     evec;
    PackRecord         record;
public:
    EigenPack(void)          = default;
    virtual ~EigenPack(void) = default;

    EigenPack(const size_t size, GridBase *grid)
    {
        resize(size, grid);
    }

    void resize(const size_t size, GridBase *grid)
    {
        eval.resize(size);
        evec.resize(size, grid);
    }

    virtual void read(const std::string fileStem, const int traj = -1)
    {
        std::string evecFilename, evalFilename;

        makeFilenames(evecFilename, evalFilename, fileStem, traj);
        XmlReader xmlReader(evalFilename);
        LOG(Message) << "Reading " << eval.size() << " eigenvalues from '" 
                     << evalFilename << "'" << std::endl;
        Grid::read(xmlReader, "evals", eval);
        basicRead(evec, evecFilename, evec.size());
    }

    virtual void write(const std::string fileStem, const int traj = -1)
    {
        std::string evecFilename, evalFilename;

        makeFilenames(evecFilename, evalFilename, fileStem, traj);
        XmlWriter xmlWriter(evalFilename);
        LOG(Message) << "Writing " << eval.size() << " eigenvalues to '" 
                     << evalFilename << "'" << std::endl;
        Grid::write(xmlWriter, "evals", eval);
        basicWrite(evecFilename, evec, evec.size());
    }
protected:
    void makeFilenames(std::string &evecFilename, std::string &evalFilename,
                       const std::string stem, const int traj = -1)
    {
        std::string t = (traj < 0) ? "" : ("." + std::to_string(traj));

        evecFilename = stem + "_evec" + t + ".bin";
        evalFilename = stem + "_eval" + t + ".xml";
    }

    template <typename T>
    void basicRead(std::vector<T> &evec, const std::string filename, 
                   const unsigned int size)
    {
        ScidacReader    binReader;

        binReader.open(filename);
        binReader.readScidacFileRecord(evec[0]._grid, record);
        for(int k = 0; k < size; ++k) 
        {
            VecRecord vecRecord;

            binReader.readScidacFieldRecord(evec[k], vecRecord);
            if (vecRecord.index != k)
            {
                HADRON_ERROR(Io, "Eigenvector " + std::to_string(k) + " has a"
                             + " wrong index (expected " + std::to_string(vecRecord.index) 
                             + ") in file '" + filename + "'");
            }
        }
        binReader.close();
    }

    template <typename T>
    void basicWrite(const std::string filename, std::vector<T> &evec, 
                    const unsigned int size)
    {
        ScidacWriter    binWriter(evec[0]._grid->IsBoss());

        binWriter.open(filename);
        binWriter.writeScidacFileRecord(evec[0]._grid, record);
        for(int k = 0; k < size; ++k) 
        {
            VecRecord vecRecord;

            vecRecord.index = k;
            vecRecord.eval  = eval[k];
            binWriter.writeScidacFieldRecord(evec[k], vecRecord);
        }
        binWriter.close();
    }
};

template <typename FineF, typename CoarseF>
class CoarseEigenPack: public EigenPack<FineF>
{
public:
    typedef CoarseF CoarseField;
public:
    std::vector<RealD>   evalCoarse;
    std::vector<CoarseF> evecCoarse;
public:
    CoarseEigenPack(void)          = default;
    virtual ~CoarseEigenPack(void) = default;

    CoarseEigenPack(const size_t sizeFine, const size_t sizeCoarse, 
                    GridBase *gridFine, GridBase *gridCoarse)
    {
        resize(sizeFine, sizeCoarse, gridFine, gridCoarse);
    }

    void resize(const size_t sizeFine, const size_t sizeCoarse, 
                GridBase *gridFine, GridBase *gridCoarse)
    {
        EigenPack<FineF>::resize(sizeFine, gridFine);
        evalCoarse.resize(sizeCoarse);
        evecCoarse.resize(sizeCoarse, gridCoarse);
    }

    void readFine(const std::string fileStem, const int traj = -1)
    {
        std::string evecFineFilename, evalFineFilename;
        std::string evecCoarseFilename, evalCoarseFilename;

        this->makeFilenames(evecFineFilename, evalFineFilename, 
                            fileStem + "_fine", traj);
        XmlReader xmlFineReader(evalFineFilename);
        LOG(Message) << "Reading " << this->eval.size() << " fine eigenvalues from '" 
                     << evalFineFilename << "'" << std::endl;
        Grid::read(xmlFineReader, "evals", this->eval);
        LOG(Message) << "Reading " << this->evec.size() << " fine eigenvectors from '" 
                     << evecFineFilename << "'" << std::endl;
        this->basicRead(this->evec, evecFineFilename, this->evec.size());
    }

    void readCoarse(const std::string fileStem, const int traj = -1)
    {
        std::string evecCoarseFilename, evalCoarseFilename;

        this->makeFilenames(evecCoarseFilename, evalCoarseFilename, 
                            fileStem + "_coarse", traj);
        XmlReader xmlCoarseReader(evalCoarseFilename);
        LOG(Message) << "Reading " << evalCoarse.size() << " coarse eigenvalues from '" 
                     << evalCoarseFilename << "'" << std::endl;
        Grid::read(xmlCoarseReader, "evals", evalCoarse);
        LOG(Message) << "Reading " << evecCoarse.size() << " coarse eigenvectors from '" 
                     << evecCoarseFilename << "'" << std::endl;
        this->basicRead(evecCoarse, evecCoarseFilename, evecCoarse.size());
    }

    virtual void read(const std::string fileStem, const int traj = -1)
    {
        readFine(fileStem, traj);
        readCoarse(fileStem, traj);
    }

    void writeFine(const std::string fileStem, const int traj = -1)
    {
        std::string evecFineFilename, evalFineFilename;

        this->makeFilenames(evecFineFilename, evalFineFilename, 
                            fileStem + "_fine", traj);
        XmlWriter xmlFineWriter(evalFineFilename);
        LOG(Message) << "Writing " << this->eval.size() << " fine eigenvalues to '" 
                     << evalFineFilename << "'" << std::endl;
        Grid::write(xmlFineWriter, "evals", this->eval);
        LOG(Message) << "Writing " << this->evec.size() << " fine eigenvectors to '" 
                     << evecFineFilename << "'" << std::endl;
        this->basicWrite(evecFineFilename, this->evec, this->evec.size());
    }

    void writeCoarse(const std::string fileStem, const int traj = -1)
    {
        std::string evecCoarseFilename, evalCoarseFilename;

        this->makeFilenames(evecCoarseFilename, evalCoarseFilename,
                            fileStem + "_coarse", traj);
        XmlWriter xmlCoarseWriter(evalCoarseFilename);
        LOG(Message) << "Writing " << evalCoarse.size() << " coarse eigenvalues to '" 
                     << evalCoarseFilename << "'" << std::endl;
        Grid::write(xmlCoarseWriter, "evals", evalCoarse);
        LOG(Message) << "Writing " << evecCoarse.size() << " coarse eigenvectors to '" 
                     << evecCoarseFilename << "'" << std::endl;
        this->basicWrite(evecCoarseFilename, evecCoarse, evecCoarse.size());
    }
    
    virtual void write(const std::string fileStem, const int traj = -1)
    {
        writeFine(fileStem, traj);
        writeCoarse(fileStem, traj);
    }
};

template <typename FImpl>
using FermionEigenPack = EigenPack<typename FImpl::FermionField>;

template <typename FImpl, int nBasis>
using CoarseFermionEigenPack = CoarseEigenPack<
    typename FImpl::FermionField,
    typename LocalCoherenceLanczos<typename FImpl::SiteSpinor, 
                                   typename FImpl::SiteComplex, 
                                   nBasis>::CoarseField>;

END_HADRONS_NAMESPACE

#endif // Hadrons_EigenPack_hpp_
