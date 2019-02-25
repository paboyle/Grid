/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/EigenPack.hpp

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
#ifndef Hadrons_EigenPack_hpp_
#define Hadrons_EigenPack_hpp_

#include <Hadrons/Global.hpp>
#include <Grid/algorithms/iterative/Deflation.h>
#include <Grid/algorithms/iterative/LocalCoherenceLanczos.h>

BEGIN_HADRONS_NAMESPACE

// Lanczos type
#ifndef HADRONS_DEFAULT_LANCZOS_NBASIS
#define HADRONS_DEFAULT_LANCZOS_NBASIS 60
#endif

#define HADRONS_DUMP_EP_METADATA(record) \
LOG(Message) << "Eigenpack metadata:" << std::endl;\
LOG(Message) << "* operator" << std::endl;\
LOG(Message) << (record).operatorXml << std::endl;\
LOG(Message) << "* solver" << std::endl;\
LOG(Message) << (record).solverXml << std::endl;

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

namespace EigenPackIo
{
    inline void readHeader(PackRecord &record, ScidacReader &binReader)
    {
        std::string recordXml;

        binReader.readLimeObject(recordXml, SCIDAC_FILE_XML);
        XmlReader xmlReader(recordXml, true, "eigenPackPar");
        xmlReader.push();
        xmlReader.readCurrentSubtree(record.operatorXml);
        xmlReader.nextElement();
        xmlReader.readCurrentSubtree(record.solverXml);
    }

    template <typename T, typename TIo = T>
    void readElement(T &evec, RealD &eval, const unsigned int index,
                     ScidacReader &binReader, TIo *ioBuf = nullptr)
    {
        VecRecord vecRecord;

        LOG(Message) << "Reading eigenvector " << index << std::endl;
        if (ioBuf == nullptr)
        {
            binReader.readScidacFieldRecord(evec, vecRecord);
        }
        else
        {
            binReader.readScidacFieldRecord(*ioBuf, vecRecord);
            precisionChange(evec, *ioBuf);
        }
        if (vecRecord.index != index)
        {
            HADRONS_ERROR(Io, "Eigenvector " + std::to_string(index) + " has a"
                            + " wrong index (expected " + std::to_string(vecRecord.index) 
                            + ")");
        }
        eval = vecRecord.eval;
    }

    template <typename T, typename TIo = T>
    static void readPack(std::vector<T> &evec, std::vector<RealD> &eval,
                         PackRecord &record, const std::string filename, 
                         const unsigned int size, bool multiFile, 
                         GridBase *gridIo = nullptr)
    {
        std::unique_ptr<TIo> ioBuf{nullptr};
        ScidacReader         binReader;

        if (typeHash<T>() != typeHash<TIo>())
        {
            if (gridIo == nullptr)
            {
                HADRONS_ERROR(Definition, 
                              "I/O type different from vector type but null I/O grid passed");
            }
            ioBuf.reset(new TIo(gridIo));
        }
        if (multiFile)
        {
            std::string fullFilename;

            for(int k = 0; k < size; ++k) 
            {
                fullFilename = filename + "/v" + std::to_string(k) + ".bin";
                binReader.open(fullFilename);
                readHeader(record, binReader);
                readElement(evec[k], eval[k], k, binReader, ioBuf.get());
                binReader.close();
            }
        }
        else
        {
            binReader.open(filename);
            readHeader(record, binReader);
            for(int k = 0; k < size; ++k) 
            {
                readElement(evec[k], eval[k], k, binReader, ioBuf.get());
            }
            binReader.close();
        }
    }

    inline void writeHeader(ScidacWriter &binWriter, PackRecord &record)
    {
        XmlWriter xmlWriter("", "eigenPackPar");

        xmlWriter.pushXmlString(record.operatorXml);
        xmlWriter.pushXmlString(record.solverXml);
        binWriter.writeLimeObject(1, 1, xmlWriter, "parameters", SCIDAC_FILE_XML);
    }

    template <typename T, typename TIo = T>
    void writeElement(ScidacWriter &binWriter, T &evec, RealD &eval, 
                      const unsigned int index, TIo *ioBuf, 
                      T *testBuf = nullptr)
    {
        VecRecord vecRecord;

        LOG(Message) << "Writing eigenvector " << index << std::endl;
        vecRecord.eval  = eval;
        vecRecord.index = index;
        if ((ioBuf == nullptr) || (testBuf == nullptr))
        {
            binWriter.writeScidacFieldRecord(evec, vecRecord, DEFAULT_ASCII_PREC);
        }
        else
        {
            precisionChange(*ioBuf, evec);
            precisionChange(*testBuf, *ioBuf);
            *testBuf -= evec;
            LOG(Message) << "Precision diff norm^2 " << norm2(*testBuf) << std::endl;
            binWriter.writeScidacFieldRecord(*ioBuf, vecRecord, DEFAULT_ASCII_PREC);
        }   
    }
    
    template <typename T, typename TIo = T>
    static void writePack(const std::string filename, std::vector<T> &evec, 
                          std::vector<RealD> &eval, PackRecord &record, 
                          const unsigned int size, bool multiFile, 
                          GridBase *gridIo = nullptr)
    {
        GridBase             *grid = evec[0]._grid;
        std::unique_ptr<TIo> ioBuf{nullptr}; 
        std::unique_ptr<T>   testBuf{nullptr};
        ScidacWriter         binWriter(grid->IsBoss());

        if (typeHash<T>() != typeHash<TIo>())
        {
            if (gridIo == nullptr)
            {
                HADRONS_ERROR(Definition, 
                              "I/O type different from vector type but null I/O grid passed");
            }
            ioBuf.reset(new TIo(gridIo));
            testBuf.reset(new T(grid));
        }
        if (multiFile)
        {
            std::string fullFilename;

            for(int k = 0; k < size; ++k) 
            {
                fullFilename = filename + "/v" + std::to_string(k) + ".bin";

                makeFileDir(fullFilename, grid);
                binWriter.open(fullFilename);
                writeHeader(binWriter, record);
                writeElement(binWriter, evec[k], eval[k], k, ioBuf.get(), testBuf.get());
                binWriter.close();
            }
        }
        else
        {
            makeFileDir(filename, grid);
            binWriter.open(filename);
            writeHeader(binWriter, record);
            for(int k = 0; k < size; ++k) 
            {
                writeElement(binWriter, evec[k], eval[k], k, ioBuf.get(), testBuf.get());
            }
            binWriter.close();
        }
    }
}

template <typename F>
class BaseEigenPack
{
public:
    typedef F Field;
public:
    std::vector<RealD> eval;
    std::vector<F>     evec;
    PackRecord         record;
public:
    BaseEigenPack(void)          = default;
    BaseEigenPack(const size_t size, GridBase *grid)
    {
        resize(size, grid);
    }
    virtual ~BaseEigenPack(void) = default;
    void resize(const size_t size, GridBase *grid)
    {
        eval.resize(size);
        evec.resize(size, grid);
    }
};

template <typename F, typename FIo = F>
class EigenPack: public BaseEigenPack<F>
{
public:
    typedef F   Field;
    typedef FIo FieldIo;
public:
    EigenPack(void)          = default;
    virtual ~EigenPack(void) = default;

    EigenPack(const size_t size, GridBase *grid, GridBase *gridIo = nullptr)
    : BaseEigenPack<F>(size, grid)
    {
        if (typeHash<F>() != typeHash<FIo>())
        {
            if (gridIo == nullptr)
            {
                HADRONS_ERROR(Definition, 
                              "I/O type different from vector type but null I/O grid passed");
            }
        }
        gridIo_ = gridIo;
    }

    virtual void read(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        EigenPackIo::readPack<F, FIo>(this->evec, this->eval, this->record, 
                                      evecFilename(fileStem, traj, multiFile), 
                                      this->evec.size(), multiFile, gridIo_);
        HADRONS_DUMP_EP_METADATA(this->record);
    }

    virtual void write(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        EigenPackIo::writePack<F, FIo>(evecFilename(fileStem, traj, multiFile), 
                                       this->evec, this->eval, this->record, 
                                       this->evec.size(), multiFile, gridIo_);
    }
protected:
    std::string evecFilename(const std::string stem, const int traj, const bool multiFile)
    {
        std::string t = (traj < 0) ? "" : ("." + std::to_string(traj));

        if (multiFile)
        {
            return stem + t;
        }
        else
        {
            return stem + t + ".bin";
        }
    }
protected:
    GridBase *gridIo_;
};

template <typename FineF, typename CoarseF, 
          typename FineFIo = FineF, typename CoarseFIo = CoarseF>
class CoarseEigenPack: public EigenPack<FineF, FineFIo>
{
public:
    typedef CoarseF   CoarseField;
    typedef CoarseFIo CoarseFieldIo;
public:      
    std::vector<CoarseF> evecCoarse;
    std::vector<RealD>   evalCoarse;
public:
    CoarseEigenPack(void)          = default;
    virtual ~CoarseEigenPack(void) = default;

    CoarseEigenPack(const size_t sizeFine, const size_t sizeCoarse, 
                    GridBase *gridFine, GridBase *gridCoarse,
                    GridBase *gridFineIo = nullptr, 
                    GridBase *gridCoarseIo = nullptr)
    {
        if (typeHash<FineF>() != typeHash<FineFIo>())
        {
            if (gridFineIo == nullptr)
            {
                HADRONS_ERROR(Definition, 
                              "Fine I/O type different from vector type but null fine I/O grid passed");
            }
        }
        if (typeHash<CoarseF>() != typeHash<CoarseFIo>())
        {
            if (gridCoarseIo == nullptr)
            {
                HADRONS_ERROR(Definition, 
                              "Coarse I/O type different from vector type but null coarse I/O grid passed");
            }
        }
        this->gridIo_ = gridFineIo;
        gridCoarseIo_ = gridCoarseIo;
        resize(sizeFine, sizeCoarse, gridFine, gridCoarse);
    }

    void resize(const size_t sizeFine, const size_t sizeCoarse, 
                GridBase *gridFine, GridBase *gridCoarse)
    {
        EigenPack<FineF, FineFIo>::resize(sizeFine, gridFine);
        evalCoarse.resize(sizeCoarse);
        evecCoarse.resize(sizeCoarse, gridCoarse);
    }

    void readFine(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        EigenPack<FineF, FineFIo>::read(fileStem + "_fine", multiFile, traj);
    }

    void readCoarse(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        PackRecord dummy;

        EigenPackIo::readPack<CoarseF, CoarseFIo>(evecCoarse, evalCoarse, dummy, 
                              this->evecFilename(fileStem + "_coarse", traj, multiFile), 
                              evecCoarse.size(), multiFile, gridCoarseIo_);
    }

    virtual void read(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        readFine(fileStem, multiFile, traj);
        readCoarse(fileStem, multiFile, traj);
    }

    void writeFine(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        EigenPack<FineF, FineFIo>::write(fileStem + "_fine", multiFile, traj);
    }

    void writeCoarse(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        EigenPackIo::writePack<CoarseF, CoarseFIo>(this->evecFilename(fileStem + "_coarse", traj, multiFile), 
                                                   evecCoarse, evalCoarse, this->record, 
                                                   evecCoarse.size(), multiFile, gridCoarseIo_);
    }
    
    virtual void write(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        writeFine(fileStem, multiFile, traj);
        writeCoarse(fileStem, multiFile, traj);
    }
private:
    GridBase *gridCoarseIo_;
};

template <typename FImpl>
using BaseFermionEigenPack = BaseEigenPack<typename FImpl::FermionField>;

template <typename FImpl, typename FImplIo = FImpl>
using FermionEigenPack = EigenPack<typename FImpl::FermionField, typename FImplIo::FermionField>;

template <typename FImpl, int nBasis, typename FImplIo = FImpl>
using CoarseFermionEigenPack = CoarseEigenPack<
    typename FImpl::FermionField,
    typename LocalCoherenceLanczos<typename FImpl::SiteSpinor, 
                                   typename FImpl::SiteComplex, 
                                   nBasis>::CoarseField,
    typename FImplIo::FermionField,
    typename LocalCoherenceLanczos<typename FImplIo::SiteSpinor, 
                                   typename FImplIo::SiteComplex, 
                                   nBasis>::CoarseField>;

#undef HADRONS_DUMP_EP_METADATA

END_HADRONS_NAMESPACE

#endif // Hadrons_EigenPack_hpp_
