/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/algorithms/multigrid/MultiGridIO.h

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

//////////////////////////////////////////////////////////////////////
// Multigrid setup I/O shared by the PVdagM and HDCG chains: a set of
// bare-vector scidac records (a near-null basis, or a set of
// eigenvectors with their eigenvalues beside them in XML).
//
// The coarse OPERATOR is not written: coarsening is cheap relative
// to producing the basis, and is redone.
//
// `control` is the scidac binary-IO control word.  The library default
// is lexicographic + aggregated; pass 0 to read files written by the
// pre-2026 HDCG drivers, which used control 0.
//////////////////////////////////////////////////////////////////////
static const int MultiGridIODefaultControl = BinaryIO::BINARYIO_LEXICOGRAPHIC|BinaryIO::BINARYIO_AGGREGATE;

template <class Field>
void saveSubspace(std::vector<Field> &subspace, std::string const fname,
		  int control = MultiGridIODefaultControl)
{
#ifdef HAVE_LIME
  Grid::emptyUserRecord record;
  Grid::ScidacWriter SW(subspace[0].Grid()->IsBoss());
  SW.open(fname);
  for (int k = 0; k < (int)subspace.size(); k++) {
    SW.writeScidacFieldRecord(subspace[k], record, 0, control);
  }
  SW.close();
#endif
}
template <class Field>
void loadSubspace(std::vector<Field> &subspace, std::string const fname,
		  int control = MultiGridIODefaultControl)
{
#ifdef HAVE_LIME
  Grid::emptyUserRecord record;
  Grid::ScidacReader SR;
  SR.open(fname);
  for (int k = 0; k < (int)subspace.size(); k++) {
    SR.readScidacFieldRecord(subspace[k], record, control);
  }
  SR.close();
#endif
}

//////////////////////////////////////////////////////////////////////
// Eigenpairs: the vectors as a subspace file, the values in
// <fname>.evals.xml.  evec.size() fixes how many are read.
//////////////////////////////////////////////////////////////////////
template <class Field>
void saveEigenpairs(std::vector<Field> &evec, std::vector<RealD> &eval, std::string const fname,
		    int control = MultiGridIODefaultControl)
{
  GRID_ASSERT(evec.size()==eval.size());
  saveSubspace(evec, fname, control);
  if ( evec[0].Grid()->IsBoss() ) {
    XmlWriter WR(fname+".evals.xml");
    write(WR,"evals",eval);
  }
}
// The eigenvalue file sets the count; the vectors are allocated on grid.
template <class Field>
void loadEigenpairs(std::vector<Field> &evec, std::vector<RealD> &eval, std::string const fname,
		    GridBase *grid, int control = MultiGridIODefaultControl)
{
  XmlReader RD(fname+".evals.xml");
  read(RD,"evals",eval);
  evec.clear();
  evec.resize(eval.size(), Field(grid));
  loadSubspace(evec, fname, control);
}

NAMESPACE_END(Grid);
