/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: Grid/algorithms/blas/A2ASpatialSum.h

    Copyright (C) 2025

Author: Peter Boyle <pboyle@bnl.gov>

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
#pragma once

NAMESPACE_BEGIN(Grid);

/*
  A2ASpatialSum

  Replaces the scalar spatial accumulation loop in A2A extended meson field
  contractions with a batched GEMM over local time slices, enabling GPU offload.

  Given:
    leftv[N_i][osite]    — conjugated left SpinColourVectors (SIMD-packed)
    loopRight[N_j][osite]— type-contracted right SpinColourVectors (SIMD-packed)

  Computes:
    EMF[i,j,t] = sum_{x,s,c} leftv[i][x,t,s,c] * loopRight[j][x,t,s,c]

  via batched GEMM over nt local time slices, then GlobalSumVector across MPI.

  Memory layout (all C row-major):
    W_buf  [nt][N_i][nxyz*Nsc]  — W[t][i][x*Nsc+sc]  = leftv[i] at (x,t)
    LR_buf [nt][N_j][nxyz*Nsc]  — LR[t][j][x*Nsc+sc] = loopRight[j] at (x,t)
    EMF_buf[nt][N_j][N_i]       — column-major result; EMF[i,j,t] = EMF_buf[t][j][i]

  BLAS call (column-major, OP_T on A so A is read as W[i][k]):
    C = A^T * B  where A=W[N_i×K C-row], B=LR[N_j×K C-row], C=[N_j×N_i C-row]
    → C[i,j] = sum_k W[i][k] * LR[j][k] = EMF[i,j] ✓
*/
template<class vobj>
class A2ASpatialSum
{
public:
  typedef typename vobj::scalar_type   scalar;
  typedef typename vobj::scalar_object sobj;

  GridBase *grid;
  int N_i, N_j;
  int nt, nxyz, Nsc;

  deviceVector<scalar>   W_buf;
  deviceVector<scalar>   LR_buf;
  deviceVector<scalar>   EMF_buf;
  deviceVector<scalar *> W_ptrs;
  deviceVector<scalar *> LR_ptrs;
  deviceVector<scalar *> EMF_ptrs;

  A2ASpatialSum() : grid(nullptr), N_i(0), N_j(0), nt(0), nxyz(0), Nsc(0) {}

  void Allocate(int _N_i, int _N_j, GridBase *_grid)
  {
    grid = _grid;
    N_i  = _N_i;
    N_j  = _N_j;
    Coordinate ldims = grid->LocalDimensions();
    nt   = ldims[grid->Nd() - 1];
    nxyz = grid->lSites() / nt;
    Nsc  = sizeof(sobj) / sizeof(scalar);

    W_buf.resize(nt * N_i * nxyz * Nsc);
    LR_buf.resize(nt * N_j * nxyz * Nsc);
    EMF_buf.resize(nt * N_j * N_i);

    // Build persistent batch pointer arrays
    W_ptrs.resize(nt);
    LR_ptrs.resize(nt);
    EMF_ptrs.resize(nt);
    scalar *Wh   = &W_buf[0];
    scalar *LRh  = &LR_buf[0];
    scalar *EMFh = &EMF_buf[0];
    int lN_i = N_i, lN_j = N_j, lnxyz = nxyz, lNsc = Nsc;
    for (int t = 0; t < nt; t++) {
      acceleratorPut(W_ptrs[t],   Wh   + t * lN_i * lnxyz * lNsc);
      acceleratorPut(LR_ptrs[t],  LRh  + t * lN_j * lnxyz * lNsc);
      acceleratorPut(EMF_ptrs[t], EMFh + t * lN_j * lN_i);
    }
  }

  void PackLeft(const std::vector<Lattice<vobj>> &leftv)
  {
    GRID_ASSERT((int)leftv.size() == N_i);
    PackVectors(leftv, &W_buf[0], N_i);
  }

  void PackRight(const std::vector<Lattice<vobj>> &loopRight)
  {
    GRID_ASSERT((int)loopRight.size() == N_j);
    PackVectors(loopRight, &LR_buf[0], N_j);
  }

private:
  // Pack vecs[N] lattice fields into buf[nt][N][nxyz*Nsc], extracting all SIMD lanes.
  void PackVectors(const std::vector<Lattice<vobj>> &vecs, scalar *buf, int N)
  {
    int nd     = grid->_ndimension;
    int osites = grid->oSites();
    int Nsimd  = vobj::Nsimd();
    int lN     = N;
    int lNsc   = Nsc;
    int lnxyz  = nxyz;
    Coordinate rdimensions = grid->_rdimensions;
    Coordinate ldims       = grid->LocalDimensions();
    Coordinate simd        = grid->_simd_layout;

    for (int n = 0; n < N; n++) {
      autoView(src_v, vecs[n], AcceleratorRead);
      accelerator_for(sf, osites, Nsimd, {
#ifdef GRID_SIMT
        {
          int lane = acceleratorSIMTlane(Nsimd);
#else
          for (int lane = 0; lane < Nsimd; lane++) {
#endif
          Coordinate icoor(nd), ocoor(nd), lcoor(nd);
          Lexicographic::CoorFromIndex(icoor, lane, simd);
          Lexicographic::CoorFromIndex(ocoor, sf, rdimensions);
          for (int d = 0; d < nd; d++)
            lcoor[d] = rdimensions[d] * icoor[d] + ocoor[d];

          int     l_t = lcoor[nd - 1];
          Coordinate xyz_coor = lcoor;
          xyz_coor[nd - 1] = 0;
          int64_t l_xyz;
          Lexicographic::IndexFromCoor(xyz_coor, l_xyz, ldims);

          sobj    data   = extractLane(lane, src_v[sf]);
          scalar *data_s = (scalar *)&data;

          int64_t base = (int64_t)l_t * lN * lnxyz * lNsc
                       + (int64_t)n   * lnxyz * lNsc
                       + l_xyz * lNsc;
          for (int sc = 0; sc < lNsc; sc++)
            buf[base + sc] = data_s[sc];
        }
      });
    }
  }

public:

  // Batched GEMM + MPI reduction → result[nt_global][N_i][N_j]
  //
  // BLAS (column-major, OP_T on A):
  //   C[N_j×N_i] = A^T[N_i×K] * B[N_j×K]    with K=nxyz*Nsc
  //   reading A as C row-major [N_i][K] and B as C row-major [N_j][K]
  //   → C[i,j] = sum_k W[i,k] * LR[j,k] = EMF[i,j] ✓
  void Sum(Eigen::Tensor<ComplexD, 3> &result)
  {
    GridBLAS BLAS;

    int K = nxyz * Nsc;
    BLAS.gemmBatched(GridBLAS_OP_T, GridBLAS_OP_N,
                     N_i, N_j, K,
                     scalar(1.0),
                     W_ptrs,
                     LR_ptrs,
                     scalar(0.0),
                     EMF_ptrs);
    BLAS.synchronise();

    // Copy from device and distribute into global-t layout
    int nt_global = result.dimension(0);
    int nd        = grid->Nd();
    int lt_start  = grid->LocalStarts()[nd - 1];

    std::vector<scalar> host_emf(nt * N_j * N_i);
    acceleratorCopyFromDevice(&EMF_buf[0], host_emf.data(),
                              nt * N_j * N_i * sizeof(scalar));

    // EMF_buf[t][j*N_i + i] = EMF[i,j] for local t
    std::vector<scalar> global_emf(nt_global * N_i * N_j, scalar(0.0));
    for (int lt = 0; lt < nt; lt++) {
      int gt = lt + lt_start;
      for (int i = 0; i < N_i; i++)
      for (int j = 0; j < N_j; j++)
        global_emf[gt * N_i * N_j + i * N_j + j] = host_emf[lt * N_j * N_i + j * N_i + i];
    }
    grid->GlobalSumVector(global_emf.data(), nt_global * N_i * N_j);

    for (int gt = 0; gt < nt_global; gt++)
    for (int i  = 0; i  < N_i;      i++)
    for (int j  = 0; j  < N_j;      j++)
      result(gt, i, j) = global_emf[gt * N_i * N_j + i * N_j + j];
  }
};

NAMESPACE_END(Grid);
