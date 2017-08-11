    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/ImplicitlyRestartedLanczos.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Chulwoo Jung <chulwoo@bnl.gov>

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
#ifndef GRID_LANC_H
#define GRID_LANC_H

#include <string.h>		//memset

#ifdef USE_LAPACK
#ifdef USE_MKL
#include<mkl_lapack.h>
#else
void LAPACK_dstegr (char *jobz, char *range, int *n, double *d, double *e,
		    double *vl, double *vu, int *il, int *iu, double *abstol,
		    int *m, double *w, double *z, int *ldz, int *isuppz,
		    double *work, int *lwork, int *iwork, int *liwork,
		    int *info);
//#include <lapacke/lapacke.h>
#endif
#endif

#include <Grid/algorithms/densematrix/DenseMatrix.h>
//#include <Grid/algorithms/iterative/EigenSort.h>

// eliminate temorary vector in calc()
#define MEM_SAVE

namespace Grid
{

  struct Bisection
  {

#if 0
    static void get_eig2 (int row_num, std::vector < RealD > &ALPHA,
			  std::vector < RealD > &BETA,
			  std::vector < RealD > &eig)
    {
      int i, j;
        std::vector < RealD > evec1 (row_num + 3);
        std::vector < RealD > evec2 (row_num + 3);
      RealD eps2;
        ALPHA[1] = 0.;
        BETHA[1] = 0.;
      for (i = 0; i < row_num - 1; i++)
	{
	  ALPHA[i + 1] = A[i * (row_num + 1)].real ();
	  BETHA[i + 2] = A[i * (row_num + 1) + 1].real ();
	}
      ALPHA[row_num] = A[(row_num - 1) * (row_num + 1)].real ();
        bisec (ALPHA, BETHA, row_num, 1, row_num, 1e-10, 1e-10, evec1, eps2);
        bisec (ALPHA, BETHA, row_num, 1, row_num, 1e-16, 1e-16, evec2, eps2);

      // Do we really need to sort here?
      int begin = 1;
      int end = row_num;
      int swapped = 1;
      while (swapped)
	{
	  swapped = 0;
	  for (i = begin; i < end; i++)
	    {
	      if (mag (evec2[i]) > mag (evec2[i + 1]))
		{
		  swap (evec2 + i, evec2 + i + 1);
		  swapped = 1;
		}
	    }
	  end--;
	  for (i = end - 1; i >= begin; i--)
	    {
	      if (mag (evec2[i]) > mag (evec2[i + 1]))
		{
		  swap (evec2 + i, evec2 + i + 1);
		  swapped = 1;
		}
	    }
	  begin++;
	}

      for (i = 0; i < row_num; i++)
	{
	  for (j = 0; j < row_num; j++)
	    {
	      if (i == j)
		H[i * row_num + j] = evec2[i + 1];
	      else
		H[i * row_num + j] = 0.;
	    }
	}
    }
#endif

    static void bisec (std::vector < RealD > &c,
		       std::vector < RealD > &b,
		       int n,
		       int m1,
		       int m2,
		       RealD eps1,
		       RealD relfeh, std::vector < RealD > &x, RealD & eps2)
    {
      std::vector < RealD > wu (n + 2);

      RealD h, q, x1, xu, x0, xmin, xmax;
      int i, a, k;

      b[1] = 0.0;
      xmin = c[n] - fabs (b[n]);
      xmax = c[n] + fabs (b[n]);
      for (i = 1; i < n; i++)
	{
	  h = fabs (b[i]) + fabs (b[i + 1]);
	  if (c[i] + h > xmax)
	    xmax = c[i] + h;
	  if (c[i] - h < xmin)
	    xmin = c[i] - h;
	}
      xmax *= 2.;

      eps2 = relfeh * ((xmin + xmax) > 0.0 ? xmax : -xmin);
      if (eps1 <= 0.0)
	eps1 = eps2;
      eps2 = 0.5 * eps1 + 7.0 * (eps2);
      x0 = xmax;
      for (i = m1; i <= m2; i++)
	{
	  x[i] = xmax;
	  wu[i] = xmin;
	}

      for (k = m2; k >= m1; k--)
	{
	  xu = xmin;
	  i = k;
	  do
	    {
	      if (xu < wu[i])
		{
		  xu = wu[i];
		  i = m1 - 1;
		}
	      i--;
	    }
	  while (i >= m1);
	  if (x0 > x[k])
	    x0 = x[k];
	  while ((x0 - xu) > 2 * relfeh * (fabs (xu) + fabs (x0)) + eps1)
	    {
	      x1 = (xu + x0) / 2;

	      a = 0;
	      q = 1.0;
	      for (i = 1; i <= n; i++)
		{
		  q =
		    c[i] - x1 -
		    ((q != 0.0) ? b[i] * b[i] / q : fabs (b[i]) / relfeh);
		  if (q < 0)
		    a++;
		}
//      printf("x1=%0.14e a=%d\n",x1,a);
	      if (a < k)
		{
		  if (a < m1)
		    {
		      xu = x1;
		      wu[m1] = x1;
		    }
		  else
		    {
		      xu = x1;
		      wu[a + 1] = x1;
		      if (x[a] > x1)
			x[a] = x1;
		    }
		}
	      else
		x0 = x1;
	    }
	  printf ("x0=%0.14e xu=%0.14e k=%d\n", x0, xu, k);
	  x[k] = (x0 + xu) / 2;
	}
    }
  };

/////////////////////////////////////////////////////////////
// Implicitly restarted lanczos
/////////////////////////////////////////////////////////////


  template < class Field > class SimpleLanczos
  {

    const RealD small = 1.0e-16;
  public:
    int lock;
    int get;
    int Niter;
    int converged;

    int Nstop;			// Number of evecs checked for convergence
    int Nk;			// Number of converged sought
    int Np;			// Np -- Number of spare vecs in kryloc space
    int Nm;			// Nm -- total number of vectors


    RealD OrthoTime;

    RealD eresid;

    SortEigen < Field > _sort;

    LinearOperatorBase < Field > &_Linop;

    OperatorFunction < Field > &_poly;

    /////////////////////////
    // Constructor
    /////////////////////////
    void init (void)
    {
    };
    void Abort (int ff, DenseVector < RealD > &evals,
		DenseVector < DenseVector < RealD > >&evecs);

    SimpleLanczos (LinearOperatorBase < Field > &Linop,	// op
		   OperatorFunction < Field > &poly,	// polynmial
		   int _Nstop,	// sought vecs
		   int _Nk,	// sought vecs
		   int _Nm,	// spare vecs
		   RealD _eresid,	// resid in lmdue deficit 
		   int _Niter):	// Max iterations
     
      _Linop (Linop),
      _poly (poly),
      Nstop (_Nstop), Nk (_Nk), Nm (_Nm), eresid (_eresid), Niter (_Niter)
    {
      Np = Nm - Nk;
      assert (Np > 0);
    };

    /////////////////////////
    // Sanity checked this routine (step) against Saad.
    /////////////////////////
    void RitzMatrix (DenseVector < Field > &evec, int k)
    {

      if (1)
	return;

      GridBase *grid = evec[0]._grid;
      Field w (grid);
      std::cout << GridLogMessage << "RitzMatrix " << std::endl;
      for (int i = 0; i < k; i++)
	{
	  _Linop.HermOp (evec[i], w);
//      _poly(_Linop,evec[i],w);
	  std::cout << GridLogMessage << "[" << i << "] ";
	  for (int j = 0; j < k; j++)
	    {
	      ComplexD in = innerProduct (evec[j], w);
	      if (fabs ((double) i - j) > 1)
		{
		  if (abs (in) > 1.0e-9)
		    {
		      std::cout << GridLogMessage << "oops" << std::endl;
		      abort ();
		    }
		  else
		    std::cout << GridLogMessage << " 0 ";
		}
	      else
		{
		  std::cout << GridLogMessage << " " << in << " ";
		}
	    }
	  std::cout << GridLogMessage << std::endl;
	}
    }

    void step (DenseVector < RealD > &lmd,
	       DenseVector < RealD > &lme,
	       Field & last, Field & current, Field & next, uint64_t k)
    {
      if (lmd.size () <= k)
	lmd.resize (k + Nm);
      if (lme.size () <= k)
	lme.resize (k + Nm);


//      _poly(_Linop,current,next );   // 3. wk:=Avk−βkv_{k−1}
      _Linop.HermOp (current, next);	// 3. wk:=Avk−βkv_{k−1}
      if (k > 0)
	{
	  next -= lme[k - 1] * last;
	}
//      std::cout<<GridLogMessage << "<last|next>" << innerProduct(last,next) <<std::endl;

      ComplexD zalph = innerProduct (current, next);	// 4. αk:=(wk,vk)
      RealD alph = real (zalph);

      next = next - alph * current;	// 5. wk:=wk−αkvk
//      std::cout<<GridLogMessage << "<current|next>" << innerProduct(current,next) <<std::endl;

      RealD beta = normalise (next);	// 6. βk+1 := ∥wk∥2. If βk+1 = 0 then Stop
      // 7. vk+1 := wk/βk+1
//       norm=beta;

      int interval = Nm / 100 + 1;
      if ((k % interval) == 0)
	std::
	  cout << GridLogMessage << k << " : alpha = " << zalph << " beta " <<
	  beta << std::endl;
      const RealD tiny = 1.0e-20;
      if (beta < tiny)
	{
	  std::cout << GridLogMessage << " beta is tiny " << beta << std::
	    endl;
	}
      lmd[k] = alph;
      lme[k] = beta;

    }

    void qr_decomp (DenseVector < RealD > &lmd,
		    DenseVector < RealD > &lme,
		    int Nk,
		    int Nm,
		    DenseVector < RealD > &Qt, RealD Dsh, int kmin, int kmax)
    {
      int k = kmin - 1;
      RealD x;

      RealD Fden = 1.0 / hypot (lmd[k] - Dsh, lme[k]);
      RealD c = (lmd[k] - Dsh) * Fden;
      RealD s = -lme[k] * Fden;

      RealD tmpa1 = lmd[k];
      RealD tmpa2 = lmd[k + 1];
      RealD tmpb = lme[k];

      lmd[k] = c * c * tmpa1 + s * s * tmpa2 - 2.0 * c * s * tmpb;
      lmd[k + 1] = s * s * tmpa1 + c * c * tmpa2 + 2.0 * c * s * tmpb;
      lme[k] = c * s * (tmpa1 - tmpa2) + (c * c - s * s) * tmpb;
      x = -s * lme[k + 1];
      lme[k + 1] = c * lme[k + 1];

      for (int i = 0; i < Nk; ++i)
	{
	  RealD Qtmp1 = Qt[i + Nm * k];
	  RealD Qtmp2 = Qt[i + Nm * (k + 1)];
	  Qt[i + Nm * k] = c * Qtmp1 - s * Qtmp2;
	  Qt[i + Nm * (k + 1)] = s * Qtmp1 + c * Qtmp2;
	}

      // Givens transformations
      for (int k = kmin; k < kmax - 1; ++k)
	{

	  RealD Fden = 1.0 / hypot (x, lme[k - 1]);
	  RealD c = lme[k - 1] * Fden;
	  RealD s = -x * Fden;

	  RealD tmpa1 = lmd[k];
	  RealD tmpa2 = lmd[k + 1];
	  RealD tmpb = lme[k];

	  lmd[k] = c * c * tmpa1 + s * s * tmpa2 - 2.0 * c * s * tmpb;
	  lmd[k + 1] = s * s * tmpa1 + c * c * tmpa2 + 2.0 * c * s * tmpb;
	  lme[k] = c * s * (tmpa1 - tmpa2) + (c * c - s * s) * tmpb;
	  lme[k - 1] = c * lme[k - 1] - s * x;

	  if (k != kmax - 2)
	    {
	      x = -s * lme[k + 1];
	      lme[k + 1] = c * lme[k + 1];
	    }

	  for (int i = 0; i < Nk; ++i)
	    {
	      RealD Qtmp1 = Qt[i + Nm * k];
	      RealD Qtmp2 = Qt[i + Nm * (k + 1)];
	      Qt[i + Nm * k] = c * Qtmp1 - s * Qtmp2;
	      Qt[i + Nm * (k + 1)] = s * Qtmp1 + c * Qtmp2;
	    }
	}
    }

#ifdef USE_LAPACK
#ifdef USE_MKL
#define LAPACK_INT MKL_INT
#else
#define LAPACK_INT long long
#endif
    void diagonalize_lapack (DenseVector < RealD > &lmd, DenseVector < RealD > &lme, int N1,	// all
			     int N2,	// get
			     GridBase * grid)
    {
      const int size = Nm;
      LAPACK_INT NN = N1;
      double evals_tmp[NN];
      double DD[NN];
      double EE[NN];
      for (int i = 0; i < NN; i++)
	for (int j = i - 1; j <= i + 1; j++)
	  if (j < NN && j >= 0)
	    {
	      if (i == j)
		DD[i] = lmd[i];
	      if (i == j)
		evals_tmp[i] = lmd[i];
	      if (j == (i - 1))
		EE[j] = lme[j];
	    }
      LAPACK_INT evals_found;
      LAPACK_INT lwork =
	((18 * NN) >
	 (1 + 4 * NN + NN * NN) ? (18 * NN) : (1 + 4 * NN + NN * NN));
      LAPACK_INT liwork = 3 + NN * 10;
      LAPACK_INT iwork[liwork];
      double work[lwork];
      LAPACK_INT isuppz[2 * NN];
      char jobz = 'N';		// calculate evals only
      char range = 'I';		// calculate il-th to iu-th evals
      //    char range = 'A'; // calculate all evals
      char uplo = 'U';		// refer to upper half of original matrix
      char compz = 'I';		// Compute eigenvectors of tridiagonal matrix
      int ifail[NN];
      LAPACK_INT info;
//  int total = QMP_get_number_of_nodes();
//  int node = QMP_get_node_number();
//  GridBase *grid = evec[0]._grid;
      int total = grid->_Nprocessors;
      int node = grid->_processor;
      int interval = (NN / total) + 1;
      double vl = 0.0, vu = 0.0;
      LAPACK_INT il = interval * node + 1, iu = interval * (node + 1);
      if (iu > NN)
	iu = NN;
      double tol = 0.0;
      if (1)
	{
	  memset (evals_tmp, 0, sizeof (double) * NN);
	  if (il <= NN)
	    {
	      printf ("total=%d node=%d il=%d iu=%d\n", total, node, il, iu);
#ifdef USE_MKL
	      dstegr (&jobz, &range, &NN,
#else
	      LAPACK_dstegr (&jobz, &range, &NN,
#endif
			     (double *) DD, (double *) EE, &vl, &vu, &il, &iu,	// these four are ignored if second parameteris 'A'
			     &tol,	// tolerance
			     &evals_found, evals_tmp, (double *) NULL, &NN,
			     isuppz, work, &lwork, iwork, &liwork, &info);
	      for (int i = iu - 1; i >= il - 1; i--)
		{
		  printf ("node=%d evals_found=%d evals_tmp[%d] = %g\n", node,
			  evals_found, i - (il - 1), evals_tmp[i - (il - 1)]);
		  evals_tmp[i] = evals_tmp[i - (il - 1)];
		  if (il > 1)
		    evals_tmp[i - (il - 1)] = 0.;
		}
	    }
	  {
	    grid->GlobalSumVector (evals_tmp, NN);
	  }
	}
// cheating a bit. It is better to sort instead of just reversing it, but the document of the routine says evals are sorted in increasing order. qr gives evals in decreasing order.
    }
#undef LAPACK_INT
#endif


    void diagonalize (DenseVector < RealD > &lmd,
		      DenseVector < RealD > &lme,
		      int N2, int N1, GridBase * grid)
    {

#ifdef USE_LAPACK
      const int check_lapack = 0;	// just use lapack if 0, check against lapack if 1

      if (!check_lapack)
	return diagonalize_lapack (lmd, lme, N2, N1, grid);

//      diagonalize_lapack(lmd2,lme2,Nm2,Nm,Qt,grid);
#endif
    }

#if 1
    static RealD normalise (Field & v)
    {
      RealD nn = norm2 (v);
      nn = sqrt (nn);
      v = v * (1.0 / nn);
      return nn;
    }

    void orthogonalize (Field & w, DenseVector < Field > &evec, int k)
    {
      double t0 = -usecond () / 1e6;
      typedef typename Field::scalar_type MyComplex;
      MyComplex ip;

      if (0)
	{
	  for (int j = 0; j < k; ++j)
	    {
	      normalise (evec[j]);
	      for (int i = 0; i < j; i++)
		{
		  ip = innerProduct (evec[i], evec[j]);	// are the evecs normalised? ; this assumes so.
		  evec[j] = evec[j] - ip * evec[i];
		}
	    }
	}

      for (int j = 0; j < k; ++j)
	{
	  ip = innerProduct (evec[j], w);	// are the evecs normalised? ; this assumes so.
	  w = w - ip * evec[j];
	}
      normalise (w);
      t0 += usecond () / 1e6;
      OrthoTime += t0;
    }

    void setUnit_Qt (int Nm, DenseVector < RealD > &Qt)
    {
      for (int i = 0; i < Qt.size (); ++i)
	Qt[i] = 0.0;
      for (int k = 0; k < Nm; ++k)
	Qt[k + k * Nm] = 1.0;
    }


    void calc (DenseVector < RealD > &eval, const Field & src, int &Nconv)
    {

      GridBase *grid = src._grid;
//      assert(grid == src._grid);

      std::
	cout << GridLogMessage << " -- Nk = " << Nk << " Np = " << Np << std::
	endl;
      std::cout << GridLogMessage << " -- Nm = " << Nm << std::endl;
      std::cout << GridLogMessage << " -- size of eval   = " << eval.
	size () << std::endl;

//      assert(c.size() && Nm == eval.size());

      DenseVector < RealD > lme (Nm);
      DenseVector < RealD > lmd (Nm);


      Field current (grid);
      Field last (grid);
      Field next (grid);

      Nconv = 0;

      RealD beta_k;

      // Set initial vector
      // (uniform vector) Why not src??
      //      evec[0] = 1.0;
      current = src;
      std::cout << GridLogMessage << "norm2(src)= " << norm2 (src) << std::
	endl;
      normalise (current);
      std::
	cout << GridLogMessage << "norm2(evec[0])= " << norm2 (current) <<
	std::endl;

      // Initial Nk steps
      OrthoTime = 0.;
      double t0 = usecond () / 1e6;
      RealD norm;		// sqrt norm of last vector

      uint64_t iter = 0;

      bool initted = false;
      std::vector < RealD > low (Nstop * 10);
      std::vector < RealD > high (Nstop * 10);
      RealD cont = 0.;
      while (1) {
	  cont = 0.;
	  std::vector < RealD > lme2 (Nm);
	  std::vector < RealD > lmd2 (Nm);
	  for (uint64_t k = 0; k < Nm; ++k, iter++) {
	      step (lmd, lme, last, current, next, iter);
	      last = current;
	      current = next;
	    }
	  double t1 = usecond () / 1e6;
	  std::cout << GridLogMessage << "IRL::Initial steps: " << t1 -
	    t0 << "seconds" << std::endl;
	  t0 = t1;
	  std::
	    cout << GridLogMessage << "IRL::Initial steps:OrthoTime " <<
	    OrthoTime << "seconds" << std::endl;

	  // getting eigenvalues
	  lmd2.resize (iter + 2);
	  lme2.resize (iter + 2);
	  for (uint64_t k = 0; k < iter; ++k) {
	      lmd2[k + 1] = lmd[k];
	      lme2[k + 2] = lme[k];
	    }
	  t1 = usecond () / 1e6;
	  std::cout << GridLogMessage << "IRL:: copy: " << t1 -
	    t0 << "seconds" << std::endl;
	  t0 = t1;
	  {
	    int total = grid->_Nprocessors;
	    int node = grid->_processor;
	    int interval = (Nstop / total) + 1;
	    int iu = (iter + 1) - (interval * node + 1);
	    int il = (iter + 1) - (interval * (node + 1));
	    std::vector < RealD > eval2 (iter + 3);
	    RealD eps2;
	    Bisection::bisec (lmd2, lme2, iter, il, iu, 1e-16, 1e-10, eval2,
			      eps2);
//        diagonalize(eval2,lme2,iter,Nk,grid);
	    RealD diff = 0.;
	    for (int i = il; i <= iu; i++) {
		if (initted)
		  diff =
		    fabs (eval2[i] - high[iu-i]) / (fabs (eval2[i]) +
						      fabs (high[iu-i]));
		if (initted && (diff > eresid))
		  cont = 1.;
		if (initted)
		  printf ("eval[%d]=%0.14e %0.14e, %0.14e\n", i, eval2[i],
			  high[iu-i], diff);
		high[iu-i] = eval2[i];
	      }
	    il = (interval * node + 1);
	    iu = (interval * (node + 1));
	    Bisection::bisec (lmd2, lme2, iter, il, iu, 1e-16, 1e-10, eval2,
			      eps2);
	    for (int i = il; i <= iu; i++) {
		if (initted)
		  diff =
		    fabs (eval2[i] - low[i]) / (fabs (eval2[i]) +
						fabs (low[i]));
		if (initted && (diff > eresid))
		  cont = 1.;
		if (initted)
		  printf ("eval[%d]=%0.14e %0.14e, %0.14e\n", i, eval2[i],
			  low[i], diff);
		low[i] = eval2[i];
	      }
	    t1 = usecond () / 1e6;
	    std::cout << GridLogMessage << "IRL:: diagonalize: " << t1 -
	      t0 << "seconds" << std::endl;
	    t0 = t1;
	  }

	  for (uint64_t k = 0; k < Nk; ++k) {
//          eval[k] = eval2[k];
	    }
	  if (initted)
	    {
	      grid->GlobalSumVector (&cont, 1);
	      if (cont < 1.) return;
	    }
	  initted = true;
	}

    }






/**
   There is some matrix Q such that for any vector y
   Q.e_1 = y and Q is unitary.
**/
    template < class T >
      static T orthQ (DenseMatrix < T > &Q, DenseVector < T > y)
    {
      int N = y.size ();	//Matrix Size
      Fill (Q, 0.0);
      T tau;
      for (int i = 0; i < N; i++)
	{
	  Q[i][0] = y[i];
	}
      T sig = conj (y[0]) * y[0];
      T tau0 = fabs (sqrt (sig));

      for (int j = 1; j < N; j++)
	{
	  sig += conj (y[j]) * y[j];
	  tau = abs (sqrt (sig));

	  if (abs (tau0) > 0.0)
	    {

	      T gam = conj ((y[j] / tau) / tau0);
	      for (int k = 0; k <= j - 1; k++)
		{
		  Q[k][j] = -gam * y[k];
		}
	      Q[j][j] = tau0 / tau;
	    }
	  else
	    {
	      Q[j - 1][j] = 1.0;
	    }
	  tau0 = tau;
	}
      return tau;
    }

/**
	There is some matrix Q such that for any vector y
	Q.e_k = y and Q is unitary.
**/
    template < class T >
      static T orthU (DenseMatrix < T > &Q, DenseVector < T > y)
    {
      T tau = orthQ (Q, y);
      SL (Q);
      return tau;
    }


/**
	Wind up with a matrix with the first con rows untouched

say con = 2
	Q is such that Qdag H Q has {x, x, val, 0, 0, 0, 0, ...} as 1st colum
	and the matrix is upper hessenberg
	and with f and Q appropriately modidied with Q is the arnoldi factorization

**/

    template < class T > static void Lock (DenseMatrix < T > &H,	///Hess mtx     
					   DenseMatrix < T > &Q,	///Lock Transform
					   T val,	///value to be locked
					   int con,	///number already locked
					   RealD small, int dfg, bool herm)
    {
      //ForceTridiagonal(H);

      int M = H.dim;
      DenseVector < T > vec;
      Resize (vec, M - con);

      DenseMatrix < T > AH;
      Resize (AH, M - con, M - con);
      AH = GetSubMtx (H, con, M, con, M);

      DenseMatrix < T > QQ;
      Resize (QQ, M - con, M - con);

      Unity (Q);
      Unity (QQ);

      DenseVector < T > evals;
      Resize (evals, M - con);
      DenseMatrix < T > evecs;
      Resize (evecs, M - con, M - con);

      Wilkinson < T > (AH, evals, evecs, small);

      int k = 0;
      RealD cold = abs (val - evals[k]);
      for (int i = 1; i < M - con; i++)
	{
	  RealD cnew = abs (val - evals[i]);
	  if (cnew < cold)
	    {
	      k = i;
	      cold = cnew;
	    }
	}
      vec = evecs[k];

      ComplexD tau;
      orthQ (QQ, vec);
      //orthQM(QQ,AH,vec);

      AH = Hermitian (QQ) * AH;
      AH = AH * QQ;

      for (int i = con; i < M; i++)
	{
	  for (int j = con; j < M; j++)
	    {
	      Q[i][j] = QQ[i - con][j - con];
	      H[i][j] = AH[i - con][j - con];
	    }
	}

      for (int j = M - 1; j > con + 2; j--)
	{

	  DenseMatrix < T > U;
	  Resize (U, j - 1 - con, j - 1 - con);
	  DenseVector < T > z;
	  Resize (z, j - 1 - con);
	  T nm = norm (z);
	  for (int k = con + 0; k < j - 1; k++)
	    {
	      z[k - con] = conj (H (j, k + 1));
	    }
	  normalise (z);

	  RealD tmp = 0;
	  for (int i = 0; i < z.size () - 1; i++)
	    {
	      tmp = tmp + abs (z[i]);
	    }

	  if (tmp < small / ((RealD) z.size () - 1.0))
	    {
	      continue;
	    }

	  tau = orthU (U, z);

	  DenseMatrix < T > Hb;
	  Resize (Hb, j - 1 - con, M);

	  for (int a = 0; a < M; a++)
	    {
	      for (int b = 0; b < j - 1 - con; b++)
		{
		  T sum = 0;
		  for (int c = 0; c < j - 1 - con; c++)
		    {
		      sum += H[a][con + 1 + c] * U[c][b];
		    }		//sum += H(a,con+1+c)*U(c,b);}
		  Hb[b][a] = sum;
		}
	    }

	  for (int k = con + 1; k < j; k++)
	    {
	      for (int l = 0; l < M; l++)
		{
		  H[l][k] = Hb[k - 1 - con][l];
		}
	    }			//H(Hb[k-1-con][l] , l,k);}}

	  DenseMatrix < T > Qb;
	  Resize (Qb, M, M);

	  for (int a = 0; a < M; a++)
	    {
	      for (int b = 0; b < j - 1 - con; b++)
		{
		  T sum = 0;
		  for (int c = 0; c < j - 1 - con; c++)
		    {
		      sum += Q[a][con + 1 + c] * U[c][b];
		    }		//sum += Q(a,con+1+c)*U(c,b);}
		  Qb[b][a] = sum;
		}
	    }

	  for (int k = con + 1; k < j; k++)
	    {
	      for (int l = 0; l < M; l++)
		{
		  Q[l][k] = Qb[k - 1 - con][l];
		}
	    }			//Q(Qb[k-1-con][l] , l,k);}}

	  DenseMatrix < T > Hc;
	  Resize (Hc, M, M);

	  for (int a = 0; a < j - 1 - con; a++)
	    {
	      for (int b = 0; b < M; b++)
		{
		  T sum = 0;
		  for (int c = 0; c < j - 1 - con; c++)
		    {
		      sum += conj (U[c][a]) * H[con + 1 + c][b];
		    }		//sum += conj( U(c,a) )*H(con+1+c,b);}
		  Hc[b][a] = sum;
		}
	    }

	  for (int k = 0; k < M; k++)
	    {
	      for (int l = con + 1; l < j; l++)
		{
		  H[l][k] = Hc[k][l - 1 - con];
		}
	    }			//H(Hc[k][l-1-con] , l,k);}}

	}
    }
#endif


  };

}
#endif
