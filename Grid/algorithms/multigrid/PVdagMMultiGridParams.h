/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/algorithms/multigrid/PVdagMMultiGridParams.h

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
// Serialisable parameters for the mrhs PVdagM multigrid solver chain.
//
// Named for the multigrid FORM it configures, as WilsonMGParams is for the
// Wilson MG tests; an MrhsHDCGParams will compose the same generic MG*
// sub-structs for the Hermitian chain.  Read from XML (or JSON -- the
// serialisation macros give both) via ReadPVdagMMultiGridParams below and
// printed at startup by the macro's operator<<, so every log identifies its
// own run.  The constructor defaults ARE a tuned operating point, not
// arbitrary: an unconfigured run reproduces it, so changing them changes what
// an unconfigured run does.  Smoother mmax == nstep (full GCR history).
//
// There are NO environment-variable controls anywhere in this subsystem:
// parameters live here, library-internal constants are hard defaults in their
// classes, verbosity is the --log channels.
// Research instruments (GCR coefficient recording/replay, Chebyshev and
// stationary smoother variants, power iteration) are deliberately NOT part
// of this interface: they remain programmatic, for algorithmic studies,
// not consumer knobs.
//////////////////////////////////////////////////////////////////////////////////////

// Arithmetic precision of the fine level INSIDE the preconditioner -- the
// smoother and the V-cycle's own fine residuals.  The outer Krylov is
// always fp64: the fp64/fp32 seam sits at its preconditioner call.
// Orthogonal to FineSloppyComms, which is the halo WIRE format.
GRID_SERIALIZABLE_ENUM(MGPrecision, undef, fp64, 1, fp32, 2);

// The smoother is the adaptive shifted PGCR -- the one correct route for
// this non-Hermitian chain (stationary replay and Chebyshev were explored
// and did not win here; Chebyshev remains the HERMITIAN chain's smoother).
struct MGSmootherParams : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(MGSmootherParams,
                                  RealD, Shift,
                                  int,   Nstep,
                                  int,   Mmax);
  MGSmootherParams(RealD shift=0.1, int nstep=6, int mmax=6)
    : Shift(shift), Nstep(nstep), Mmax(mmax) {};
};

struct MGCoarseSolverParams : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(MGCoarseSolverParams,
                                  RealD, Tol,
                                  int,   Order,
                                  int,   Mmax);
  MGCoarseSolverParams() : Tol(0.05), Order(200), Mmax(16) {};
};

struct MGOuterParams : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(MGOuterParams,
                                  RealD, Tol,
                                  int,   MaxIterations,
                                  int,   Mmax,
                                  int,   Nstep);
  MGOuterParams() : Tol(1.0e-8), MaxIterations(1000), Mmax(4), Nstep(8) {};
};

struct MGDenseParams : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(MGDenseParams,
                                  int, LeafSpan);   // big-leaf span in blocks; see BlockCyclicSchurInverse
  MGDenseParams() : LeafSpan(9) {};
};

struct MGSetupParams : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(MGSetupParams,
                                  std::vector<int>, Block1,          // fine -> coarse blocking
                                  std::vector<int>, Block2,          // coarse -> coarse-coarse blocking
                                  int,              CoarsenBatch,
                                  std::string,      SubspaceFile,    // scidac; empty = create from noise, no I/O
                                  int,              FineSloppyComms, // reduced-precision halo WIRE inside the preconditioner only: fp32 on the fp64 operator, bf16 on the fp32 operator
                                  MGPrecision,      FinePrecision,   // fine-level ARITHMETIC inside the preconditioner
                                  int,              RetainSubspace); // keep the RAW basis after setup (nbasis fine vectors); needed ONLY for a fixed-basis rebuild on a changed gauge field (HMC).  A valence solve never rebuilds, so 0 frees it.
  MGSetupParams()
    : Block1({2,2,3,3}), Block2({4,4,2,4}), CoarsenBatch(9),
      SubspaceFile(""), FineSloppyComms(1), FinePrecision(MGPrecision::fp64),
      RetainSubspace(0) {};
};

struct PVdagMMultiGridParams : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(PVdagMMultiGridParams,
                                  MGSetupParams,        Setup,
                                  MGSmootherParams,     FineSmoother,
                                  MGSmootherParams,     CoarseSmoother,
                                  MGCoarseSolverParams, CoarseSolver,
                                  MGOuterParams,        Outer,
                                  MGDenseParams,        Dense);
  PVdagMMultiGridParams()
    : FineSmoother  (0.1, 6, 6),
      CoarseSmoother(0.1, 2, 2) {};
};

inline void CheckValidity(const PVdagMMultiGridParams &P)
{
  GRID_ASSERT( P.Setup.Block1.size() == 4 );
  GRID_ASSERT( P.Setup.Block2.size() == 4 );
  GRID_ASSERT( P.Setup.CoarsenBatch  >= 1 );
  GRID_ASSERT( P.Setup.FinePrecision != MGPrecision::undef );
  GRID_ASSERT( P.FineSmoother.Nstep   > 0 );
  GRID_ASSERT( P.CoarseSmoother.Nstep > 0 );
  GRID_ASSERT( P.CoarseSolver.Tol     > 0.0 );
  GRID_ASSERT( P.Outer.Tol            > 0.0 );
  GRID_ASSERT( P.Dense.LeafSpan      >= 1 );
}

//////////////////////////////////////////////////////////////////////////////////////
// Read from an XML file; on a MISSING file the boss writes <file>.templ from
// the defaults and returns false so the caller can Grid_finalize() and exit
// -- the produced template is the documentation of every parameter.  An
// empty filename means "defaults": no I/O, returns true.  Always prints the
// effective parameters, so the log describes the run.
//
// Drivers fetch the filename themselves, e.g.
//   std::string pfile("");
//   if( GridCmdOptionExists(argv,argv+argc,"--multigrid-params") )
//     pfile = GridCmdOptionPayload(argv,argv+argc,"--multigrid-params");
//   if( !ReadPVdagMMultiGridParams(Params, pfile) ) { Grid_finalize(); return 0; }
//////////////////////////////////////////////////////////////////////////////////////
inline bool ReadPVdagMMultiGridParams(PVdagMMultiGridParams &P, const std::string &xmlfile)
{
  if ( xmlfile.length() ) {
    bool good;
    { std::ifstream f(xmlfile); good = f.good(); }
    if ( !good ) {
      std::cout << GridLogMessage << "PVdagMMultiGridParams: " << xmlfile << " does not exist" << std::endl;
      if ( GlobalSharedMemory::WorldRank == 0 ) {
        XmlWriter WR(xmlfile+".templ");
        write(WR, "PVdagMMultiGridParams", P);
        std::cout << GridLogMessage << "PVdagMMultiGridParams: template written to "
                  << xmlfile << ".templ" << std::endl;
      }
      return false;
    }
    XmlReader RD(xmlfile);
    read(RD, "PVdagMMultiGridParams", P);
  }
  CheckValidity(P);
  std::cout << GridLogMessage << "PVdagMMultiGridParams ("
            << (xmlfile.length() ? xmlfile : std::string("defaults")) << "):" << std::endl;
  std::cout << P << std::endl;
  return true;
}

NAMESPACE_END(Grid);
