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
    struct PackRecord
    {
        std::string operatorXml, solverXml;
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

    virtual void read(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        if (multiFile)
        {
            for(int k = 0; k < evec.size(); ++k)
            {
                basicReadSingle(evec[k], eval[k], evecFilename(fileStem, k, traj), k);
            }
        }
        else
        {
            basicRead(evec, eval, evecFilename(fileStem, -1, traj), evec.size());
        }
    }

    virtual void write(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        if (multiFile)
        {
            for(int k = 0; k < evec.size(); ++k)
            {
                basicWriteSingle(evecFilename(fileStem, k, traj), evec[k], eval[k], k);
            }
        }
        else
        {
            basicWrite(evecFilename(fileStem, -1, traj), evec, eval, evec.size());
        }
    }
protected:
    std::string evecFilename(const std::string stem, const int vec, const int traj)
    {
        std::string t = (traj < 0) ? "" : ("." + std::to_string(traj));

        if (vec == -1)
        {
            return stem + t + ".bin";
        }
        else
        {
            return stem + t + "/v" + std::to_string(vec) + ".bin";
        }
    }

    template <typename T>
    void basicRead(std::vector<T> &evec, std::vector<double> &eval,
                   const std::string filename, const unsigned int size)
    {
        ScidacReader    binReader;

        binReader.open(filename);
        binReader.skipPastObjectRecord(SCIDAC_FILE_XML);
        for(int k = 0; k < size; ++k) 
        {
            VecRecord vecRecord;

            LOG(Message) << "Reading eigenvector " << k << std::endl;
            binReader.readScidacFieldRecord(evec[k], vecRecord);
            if (vecRecord.index != k)
            {
                HADRONS_ERROR(Io, "Eigenvector " + std::to_string(k) + " has a"
                              + " wrong index (expected " + std::to_string(vecRecord.index) 
                              + ") in file '" + filename + "'");
            }
            eval[k] = vecRecord.eval;
        }
        binReader.close();
    }

    template <typename T>
    void basicReadSingle(T &evec, double &eval, const std::string filename, 
                         const unsigned int index)
    {
        ScidacReader binReader;
        VecRecord    vecRecord;

        binReader.open(filename);
        binReader.skipPastObjectRecord(SCIDAC_FILE_XML);
        LOG(Message) << "Reading eigenvector " << index << std::endl;
        binReader.readScidacFieldRecord(evec, vecRecord);
        if (vecRecord.index != index)
        {
            HADRONS_ERROR(Io, "Eigenvector " + std::to_string(index) + " has a"
                          + " wrong index (expected " + std::to_string(vecRecord.index) 
                          + ") in file '" + filename + "'");
        }
        eval = vecRecord.eval;
        binReader.close();
    }

    template <typename T>
    void basicWrite(const std::string filename, std::vector<T> &evec, 
                    const std::vector<double> &eval, const unsigned int size)
    {
        ScidacWriter binWriter(evec[0]._grid->IsBoss());
        XmlWriter    xmlWriter("", "eigenPackPar");

        makeFileDir(filename, evec[0]._grid);
        xmlWriter.pushXmlString(record.operatorXml);
        xmlWriter.pushXmlString(record.solverXml);
        binWriter.open(filename);
        binWriter.writeLimeObject(1, 1, xmlWriter, "parameters", SCIDAC_FILE_XML);
        for(int k = 0; k < size; ++k) 
        {
            VecRecord vecRecord;

            vecRecord.index = k;
            vecRecord.eval  = eval[k];
            LOG(Message) << "Writing eigenvector " << k << std::endl;
            binWriter.writeScidacFieldRecord(evec[k], vecRecord, DEFAULT_ASCII_PREC);
        }
        binWriter.close();
    }

    template <typename T>
    void basicWriteSingle(const std::string filename, T &evec, 
                          const double eval, const unsigned int index)
    {
        ScidacWriter binWriter(evec._grid->IsBoss());
        XmlWriter    xmlWriter("", "eigenPackPar");
        VecRecord    vecRecord;

        makeFileDir(filename, evec._grid);
        xmlWriter.pushXmlString(record.operatorXml);
        xmlWriter.pushXmlString(record.solverXml);
        binWriter.open(filename);
        binWriter.writeLimeObject(1, 1, xmlWriter, "parameters", SCIDAC_FILE_XML);
        vecRecord.index = index;
        vecRecord.eval  = eval;
        LOG(Message) << "Writing eigenvector " << index << std::endl;
        binWriter.writeScidacFieldRecord(evec, vecRecord, DEFAULT_ASCII_PREC);
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

    void readFine(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        if (multiFile)
        {
            for(int k = 0; k < this->evec.size(); ++k)
            {
                this->basicReadSingle(this->evec[k], this->eval[k], this->evecFilename(fileStem + "_fine", k, traj), k);
            }
        }
        else
        {
            this->basicRead(this->evec, this->eval, this->evecFilename(fileStem + "_fine", -1, traj), this->evec.size());
        }
    }

    void readCoarse(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        if (multiFile)
        {
            for(int k = 0; k < evecCoarse.size(); ++k)
            {
                this->basicReadSingle(evecCoarse[k], evalCoarse[k], this->evecFilename(fileStem + "_coarse", k, traj), k);
            }
        }
        else
        {
            this->basicRead(evecCoarse, evalCoarse, this->evecFilename(fileStem + "_coarse", -1, traj), evecCoarse.size());
        }
    }

    virtual void read(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        readFine(fileStem, multiFile, traj);
        readCoarse(fileStem, multiFile, traj);
    }

    void writeFine(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        if (multiFile)
        {
            for(int k = 0; k < this->evec.size(); ++k)
            {
                this->basicWriteSingle(this->evecFilename(fileStem + "_fine", k, traj), this->evec[k], this->eval[k], k);
            }
        }
        else
        {
            this->basicWrite(this->evecFilename(fileStem + "_fine", -1, traj), this->evec, this->eval, this->evec.size());
        }
    }

    void writeCoarse(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        if (multiFile)
        {
            for(int k = 0; k < evecCoarse.size(); ++k)
            {
                this->basicWriteSingle(this->evecFilename(fileStem + "_coarse", k, traj), evecCoarse[k], evalCoarse[k], k);
            }
        }
        else
        {
            this->basicWrite(this->evecFilename(fileStem + "_coarse", -1, traj), evecCoarse, evalCoarse, evecCoarse.size());
        }
    }
    
    virtual void write(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        writeFine(fileStem, multiFile, traj);
        writeCoarse(fileStem, multiFile, traj);
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
