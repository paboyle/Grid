/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/algorithms/multigrid/HDCGMultiGridParams.h

    Copyright (C) 2026

Author: Peter Boyle <pboyle@bnl.gov>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See the full license in the file "LICENSE" in the top level distribution
    directory
*************************************************************************************/
/*  END LEGAL */
#pragma once

NAMESPACE_BEGIN(Grid);

//////////////////////////////////////////////////////////////////////////////////////
// Parameters of the two-level mrhs HDCG (Hermitian, Schur-preconditioned)
// multigrid, HDCGMultiGrid.h.  Serialisable: read from XML, template written
// on a missing file.  Setup (blocking, batch, basis file, fine precision) is
// the same struct as the PVdagM chain's; Block2 is only used when a
// coarse-coarse level exists but must still divide the coarse lattice.
//////////////////////////////////////////////////////////////////////////////////////

GRID_SERIALIZABLE_ENUM(HDCGSubspaceMethod, undef, ChebyshevNew, 1, Chebyshev, 2, FineEvecs, 3);
GRID_SERIALIZABLE_ENUM(HDCGRefineMethod,   undef, None, 1, CG, 2, HDCG, 3);
GRID_SERIALIZABLE_ENUM(HDCGLanczosMethod,  undef, None, 1, Block, 2, Plain, 3);
GRID_SERIALIZABLE_ENUM(HDCGCoarseMethod,   undef, DeflatedCG, 1, Chebyshev, 2);
GRID_SERIALIZABLE_ENUM(HDCGSmootherType,   undef, CG, 1, Chebyshev, 2);
GRID_SERIALIZABLE_ENUM(HDCGOuterAlgorithm, undef, fPcg, 1, PrecBlockCGrQ, 2);

// Near-null basis.  ChebyshevNew filters noise through the Hermitian operator
// with the upper bound only; Chebyshev takes lo/hi/order; FineEvecs sums
// FineEvecSample consecutive fp32 fine eigenvector records per basis vector.
// A ChebyHi <= 0 means 1.1 x the power-method estimate.  FileControl is the
// scidac binary-IO control word of the basis files (0 for pre-2026 HDCG files).
struct HDCGSubspaceParams : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(HDCGSubspaceParams,
                                  HDCGSubspaceMethod, Method,
                                  RealD,       ChebyHi,
                                  RealD,       ChebyLo,
                                  int,         ChebyOrder,
                                  std::string, FineEvecFile,
                                  int,         FineEvecSample,
                                  int,         FileControl,
                                  HDCGRefineMethod, Refine,
                                  RealD,       RefineShift,         // CG: the shift of the refining CG; HDCG: of the refining operator
                                  RealD,       RefineTol,
                                  int,         RefineMaxIt,
                                  int,         RefineSmootherOrder, // HDCG refine only
                                  std::string, RefinedSubspaceFile);
  HDCGSubspaceParams()
    : Method(HDCGSubspaceMethod::ChebyshevNew), ChebyHi(0.0), ChebyLo(0.01), ChebyOrder(500),
      FineEvecFile(""), FineEvecSample(1),
      FileControl(BinaryIO::BINARYIO_LEXICOGRAPHIC|BinaryIO::BINARYIO_AGGREGATE),
      Refine(HDCGRefineMethod::None), RefineShift(1.0e-3), RefineTol(1.0e-3), RefineMaxIt(500),
      RefineSmootherOrder(11), RefinedSubspaceFile("") {};
};

// Coarse eigenvectors for the deflated coarse solve: block Lanczos at the
// solve's Nrhs, or plain Lanczos at Nrhs 1, on the same coarse operator.
// ChebyHi <= 0 means 1.1 x the power-method estimate on the coarse operator.
struct HDCGLanczosParams : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(HDCGLanczosParams,
                                  HDCGLanczosMethod, Method,
                                  int,   Nstop,
                                  int,   Nk,
                                  int,   Nm,
                                  RealD, ChebyLo,
                                  RealD, ChebyHi,
                                  int,   ChebyOrder,
                                  RealD, Tol,
                                  int,   MaxIt,
                                  std::string, EvecFile);
  HDCGLanczosParams()
    : Method(HDCGLanczosMethod::Block), Nstop(32), Nk(32), Nm(128),
      ChebyLo(0.01), ChebyHi(0.0), ChebyOrder(201), Tol(1.0e-4), MaxIt(100), EvecFile("") {};
};

struct HDCGCoarseSolverParams : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(HDCGCoarseSolverParams,
                                  HDCGCoarseMethod, Method,
                                  RealD, Tol,            // DeflatedCG
                                  int,   MaxIt,
                                  RealD, ChebyLo,        // Chebyshev polynomial inverse
                                  RealD, ChebyHi,
                                  int,   ChebyOrder);
  HDCGCoarseSolverParams()
    : Method(HDCGCoarseMethod::DeflatedCG), Tol(5.0e-2), MaxIt(3000),
      ChebyLo(1.0e-2), ChebyHi(40.0), ChebyOrder(120) {};
};

// CG: Order iterations of CG on the operator shifted by Shift.
// Chebyshev: order Order polynomial on [ChebyLo,ChebyHi] of the unshifted operator.
struct HDCGSmootherParams : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(HDCGSmootherParams,
                                  HDCGSmootherType, Type,
                                  RealD, Shift,
                                  int,   Order,
                                  RealD, ChebyLo,
                                  RealD, ChebyHi);
  HDCGSmootherParams() : Type(HDCGSmootherType::CG), Shift(2.0), Order(7), ChebyLo(2.0), ChebyHi(92.0) {};
};

struct HDCGOuterParams : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(HDCGOuterParams,
                                  RealD, Tol,
                                  int,   MaxIt,
                                  HDCGOuterAlgorithm, Algorithm);
  HDCGOuterParams() : Tol(1.0e-8), MaxIt(500), Algorithm(HDCGOuterAlgorithm::fPcg) {};
};

struct HDCGMultiGridParams : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(HDCGMultiGridParams,
                                  MGSetupParams,          Setup,
                                  HDCGSubspaceParams,     Subspace,
                                  HDCGLanczosParams,      Lanczos,
                                  HDCGCoarseSolverParams, CoarseSolver,
                                  HDCGSmootherParams,     Smoother,
                                  HDCGOuterParams,        Outer);
  HDCGMultiGridParams() {};
};

inline void CheckValidity(const HDCGMultiGridParams &P)
{
  GRID_ASSERT( P.Setup.Block1.size() == 4 );
  GRID_ASSERT( P.Setup.Block2.size() == 4 );
  GRID_ASSERT( P.Setup.CoarsenBatch  >= 1 );
  GRID_ASSERT( P.Setup.FinePrecision != MGPrecision::undef );
  GRID_ASSERT( P.Subspace.Method     != HDCGSubspaceMethod::undef );
  GRID_ASSERT( P.Subspace.Refine     != HDCGRefineMethod::undef );
  GRID_ASSERT( P.Subspace.FineEvecSample >= 1 );
  GRID_ASSERT( P.Lanczos.Method      != HDCGLanczosMethod::undef );
  GRID_ASSERT( P.Lanczos.Nstop <= P.Lanczos.Nk && P.Lanczos.Nk <= P.Lanczos.Nm );
  GRID_ASSERT( P.CoarseSolver.Method != HDCGCoarseMethod::undef );
  GRID_ASSERT( P.CoarseSolver.Tol     > 0.0 );
  GRID_ASSERT( P.Smoother.Type       != HDCGSmootherType::undef );
  GRID_ASSERT( P.Smoother.Order       > 0 );
  GRID_ASSERT( P.Outer.Tol            > 0.0 );
  GRID_ASSERT( P.Outer.Algorithm     != HDCGOuterAlgorithm::undef );
}

// As ReadPVdagMMultiGridParams: a missing file gets a template written next
// to it and false returned; an empty filename means defaults.
inline bool ReadHDCGMultiGridParams(HDCGMultiGridParams &P, const std::string &xmlfile)
{
  if ( xmlfile.length() ) {
    bool good;
    { std::ifstream f(xmlfile); good = f.good(); }
    if ( !good ) {
      std::cout << GridLogMessage << "HDCGMultiGridParams: " << xmlfile << " does not exist" << std::endl;
      if ( GlobalSharedMemory::WorldRank == 0 ) {
        XmlWriter WR(xmlfile+".templ");
        write(WR, "HDCGMultiGridParams", P);
        std::cout << GridLogMessage << "HDCGMultiGridParams: template written to "
                  << xmlfile << ".templ" << std::endl;
      }
      return false;
    }
    XmlReader RD(xmlfile);
    read(RD, "HDCGMultiGridParams", P);
  }
  CheckValidity(P);
  std::cout << GridLogMessage << "HDCGMultiGridParams ("
            << (xmlfile.length() ? xmlfile : std::string("defaults")) << "):" << std::endl;
  std::cout << P << std::endl;
  return true;
}

NAMESPACE_END(Grid);
