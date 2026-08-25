/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/communicator/RingAllReduce.h

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

/////////////////////////////////////////////////////////////////////////////
// Vector all-reduce on point-to-point only (SendToRecvFrom), no MPI
// collectives: reduce-scatter ring followed by all-gather ring.  A
// communicator-level primitive: needs only CartesianCommunicator.
//
//   RingAllReduce(comm, buf, n)           flat ring over all P ranks
//   CartesianRingAllReduce(comm, buf, n)  ring along each processor dimension
//                                         in turn (P_d ranks per ring)
//
// Why: Cray MPICH device-buffer MPI_Allreduce aborts above ~8 MB (MPI_FLOAT,
// measured 4.4 MB pass / 13.3 MB fail) and delivers 5.8 GB/s where P2P rings
// deliver ~17 GB/s.  A ring has no size cliff -- every step is one symmetric
// SendToRecvFrom of one chunk.
//
// Cost model, N bytes per rank:
//   flat      : 2(P-1) steps, 2N(P-1)/P bytes per rank -- bandwidth-optimal,
//               latency 2(P-1) x step  (P=288: 574 steps)
//   cartesian : sum_d 2(P_d-1) steps, ~2N per dimension -- few steps, each on
//               a neighbour link (3.6.4.4: 26 steps, ~8N bytes)
// Rule of thumb: cartesian for tens of MB, flat for hundreds of MB.
//
// Buffer memory space: whatever SendToRecvFrom accepts on this build (device
// under ACCELERATOR_AWARE_MPI); the reduction runs as accelerator_for on a
// deviceVector working copy, so the caller's buffer is only ever memcpy'd.
// Deterministic: summation order is fixed by rank and chunk index, so the
// result is bitwise reproducible run to run (MPI_Allreduce need not be),
// though it differs from MPI's order at rounding level.
//
// n elements of T; T must support + on the accelerator (RealF/RealD,
// ComplexF/ComplexD).  sizeof(T)*chunk must be a multiple of 4 bytes
// (SendToRecvFrom counts int32 words); true for all supported T.
/////////////////////////////////////////////////////////////////////////////

// Ring all-reduce among P ranks with given next/prev neighbours; `me` is this
// rank's position in the ring.  work has P*c elements, scratch has c.
template<class T>
void RingAllReduceCore(CartesianCommunicator *comm,
		       T *work, T *scratch, uint64_t c, int P, int me, int next, int prev)
{
  if ( P==1 ) return;
  uint64_t bytes = c*sizeof(T);
  // reduce-scatter: after P-1 steps rank me owns fully reduced chunk (me+1)%P
  for(int s=0;s<P-1;s++){
    int sendc = (me - s     + 2*P) % P;
    int recvc = (me - s - 1 + 2*P) % P;
    comm->SendToRecvFrom((void *)&work[sendc*c], next, (void *)scratch, prev, bytes);
    T *dst = &work[recvc*c];
    accelerator_for(i, c, 1, { dst[i] = dst[i] + scratch[i]; });
  }
  // all-gather: circulate the reduced chunks
  for(int s=0;s<P-1;s++){
    int sendc = (me - s + 1 + 2*P) % P;
    int recvc = (me - s     + 2*P) % P;
    comm->SendToRecvFrom((void *)&work[sendc*c], next, (void *)&work[recvc*c], prev, bytes);
  }
}

template<class T>
void RingAllReduce(CartesianCommunicator *comm, T *buf, uint64_t n)
{
  int P  = comm->ProcessorCount();
  int me = comm->ThisRank();
  if ( P==1 || n==0 ) return;
  uint64_t c = (n + P - 1)/P;
  deviceVector<T> work(c*P);
  deviceVector<T> scratch(c);
  T *w = &work[0];
  accelerator_for(i, c*P, 1, { w[i] = T(0.0); });
  acceleratorCopyDeviceToDevice((void *)buf, (void *)w, n*sizeof(T));
  RingAllReduceCore(comm, w, &scratch[0], c, P, me, (me+1)%P, (me+P-1)%P);
  acceleratorCopyDeviceToDevice((void *)w, (void *)buf, n*sizeof(T));
}

template<class T>
void CartesianRingAllReduce(CartesianCommunicator *comm, T *buf, uint64_t n)
{
  if ( comm->ProcessorCount()==1 || n==0 ) return;
  int Nd = comm->_ndimension;
  for(int d=0;d<Nd;d++){
    int P = comm->_processors[d];
    if ( P==1 ) continue;
    int me = comm->_processor_coor[d];
    int next, prev;
    comm->ShiftedRanks(d, 1, prev, next);   // (dim, shift, source, dest)
    uint64_t c = (n + P - 1)/P;
    deviceVector<T> work(c*P);
    deviceVector<T> scratch(c);
    T *w = &work[0];
    accelerator_for(i, c*P, 1, { w[i] = T(0.0); });
    acceleratorCopyDeviceToDevice((void *)buf, (void *)w, n*sizeof(T));
    RingAllReduceCore(comm, w, &scratch[0], c, P, me, next, prev);
    acceleratorCopyDeviceToDevice((void *)w, (void *)buf, n*sizeof(T));
  }
}

NAMESPACE_END(Grid);
