/* -*- Mode: C; comment-column: 22; fill-column: 79; compile-command: "gcc -o zolotarev zolotarev.c -ansi -pedantic -lm -DTEST"; -*- */
#define VERSION Source Time-stamp: <2015-05-18 16:32:08 neo>

/* These C routines evalute the optimal rational approximation to the signum
 * function for epsilon < |x| < 1 using Zolotarev's theorem.
 *
 * To obtain reliable results for high degree approximations (large n) it is
 * necessary to compute using sufficiently high precision arithmetic. To this
 * end the code has been parameterised to work with the preprocessor names
 * INTERNAL_PRECISION and PRECISION set to float, double, or long double as
 * appropriate. INTERNAL_PRECISION is used in computing the Zolotarev
 * coefficients, which are converted to PRECISION before being returned to the
 * caller. Presumably even higher precision could be obtained using GMP or
 * similar package, but bear in mind that rounding errors might also be
 * significant in evaluating the resulting polynomial. The convergence criteria
 * have been written in a precision-independent form. */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#ifndef INTERNAL_PRECISION
#define INTERNAL_PRECISION double
#endif

#include "Zolotarev.h"
#define ZOLOTAREV_INTERNAL
#undef ZOLOTAREV_DATA
#define ZOLOTAREV_DATA izd
#undef ZPRECISION
#define ZPRECISION INTERNAL_PRECISION
#include "Zolotarev.h"
#undef ZOLOTAREV_INTERNAL

/* The ANSI standard appears not to know what pi is */

#ifndef M_PI
#define M_PI ((INTERNAL_PRECISION) 3.141592653589793238462643383279502884197\
169399375105820974944592307816406286208998628034825342117068)
#endif

#define ZERO ((INTERNAL_PRECISION) 0)
#define ONE ((INTERNAL_PRECISION) 1)
#define TWO ((INTERNAL_PRECISION) 2)
#define THREE ((INTERNAL_PRECISION) 3)
#define FOUR ((INTERNAL_PRECISION) 4)
#define HALF (ONE/TWO)

/* The following obscenity seems to be the simplest (?) way to coerce the C
 * preprocessor to convert the value of a preprocessor token into a string. */

#define PP2(x) #x
#define PP1(a,b,c) a ## b(c)
#define STRINGIFY(name) PP1(PP,2,name)

/* Compute the partial fraction expansion coefficients (alpha) from the
 * factored form */
namespace Grid {
namespace Approx {

static void construct_partfrac(izd *z) {
  int dn = z -> dn, dd = z -> dd, type = z -> type;
  int j, k, da = dd + 1 + type;
  INTERNAL_PRECISION A = z -> A, *a = z -> a, *ap = z -> ap, *alpha;
  alpha = (INTERNAL_PRECISION*) malloc(da * sizeof(INTERNAL_PRECISION));
  for (j = 0; j < dd; j++)
    for (k = 0, alpha[j] = A; k < dd; k++)
      alpha[j] *=
	(k < dn ? ap[j] - a[k] : ONE) / (k == j ? ONE : ap[j] - ap[k]);
  if(type == 1)	      /* implicit pole at zero? */
    for (k = 0, alpha[dd] = A * (dn > dd ? - a[dd] : ONE); k < dd; k++) {
      alpha[dd] *= a[k] / ap[k];
      alpha[k] *= (dn > dd ? ap[k] - a[dd] : ONE) / ap[k];
    }
  alpha[da-1] = dn == da - 1 ? A : ZERO;
  z -> alpha = alpha;
  z -> da = da;
  return;
}

/* Convert factored polynomial into dense polynomial. The input is the overall
 * factor A and the roots a[i], such that p = A product(x - a[i], i = 1..d) */

static INTERNAL_PRECISION *poly_factored_to_dense(INTERNAL_PRECISION A, 
						  INTERNAL_PRECISION *a,
						  int d) {
  INTERNAL_PRECISION *p;
  int i, j;
  p = (INTERNAL_PRECISION *) malloc((d + 2) * sizeof(INTERNAL_PRECISION));
  p[0] = A;
  for (i = 0; i < d; i++) {
    p[i+1] = p[i];
    for (j = i; j > 0; j--) p[j] = p[j-1] - a[i]*p[j];
    p[0] *= - a[i];
  }
  return p;
}

/* Convert a rational function of the form R0(x) = x p(x^2)/q(x^2) (type 0) or
 * R1(x) = p(x^2)/[x q(x^2)] (type 1) into its continued fraction
 * representation. We assume that 0 <= deg(q) - deg(p) <= 1 for type 0 and 0 <=
 * deg(p) - deg(q) <= 1 for type 1. On input p and q are in factored form, and
 * deg(q) = dq, deg(p) = dp.  The output is the continued fraction coefficients
 * beta, where R(x) = beta[0] x + 1/(beta[1] x + 1/(...)).
 *
 * The method used is as follows. There are four cases to consider:
 *
 * 0.i.  Type 0, deg p = deg q
 *
 * 0.ii. Type 0, deg p = deg q - 1
 *
 * 1.i.  Type 1, deg p = deg q
 *
 * 1.ii. Type 1, deg p = deg q + 1
 *
 * and these are connected by two transformations:
 *
 * A. To obtain a continued fraction expansion of type 1 we use a single-step
 * polynomial division we find beta and r(x) such that p(x) = beta x q(x) +
 * r(x), with deg(r) = deg(q). This implies that p(x^2) = beta x^2 q(x^2) +
 * r(x^2), and thus R1(x) = x beta + r(x^2)/(x q(x^2)) = x beta + 1/R0(x)
 * with R0(x) = x q(x^2)/r(x^2).
 *
 * B. A continued fraction expansion of type 0 is obtained in a similar, but
 * not identical, manner. We use the polynomial division algorithm to compute
 * the quotient beta and the remainder r that satisfy p(x) = beta q(x) + r(x)
 * with deg(r) = deg(q) - 1. We thus have x p(x^2) = x beta q(x^2) + x r(x^2),
 * so R0(x) = x beta + x r(x^2)/q(x^2) = x beta + 1/R1(x) with R1(x) = q(x^2) /
 * (x r(x^2)).
 *
 * Note that the deg(r) must be exactly deg(q) for (A) and deg(q) - 1 for (B)
 * because p and q have disjoint roots all of multiplicity 1. This means that
 * the division algorithm requires only a single polynomial subtraction step.
 *
 * The transformations between the cases form the following finite state
 * automaton:
 *
 * +------+            +------+            +------+            +------+
 * |      |            |      | ---(A)---> |      |            |      |
 * | 0.ii | ---(B)---> | 1.ii |            | 0.i  | <---(A)--- | 1.i  |
 * |      |            |      | <---(B)--- |      |            |      |
 * +------+            +------+            +------+            +------+
 */

static INTERNAL_PRECISION *contfrac_A(INTERNAL_PRECISION *,
				      INTERNAL_PRECISION *,
				      INTERNAL_PRECISION *,
				      INTERNAL_PRECISION *, int, int);

static INTERNAL_PRECISION *contfrac_B(INTERNAL_PRECISION *,
				      INTERNAL_PRECISION *,
				      INTERNAL_PRECISION *,
				      INTERNAL_PRECISION *, int, int);

static void construct_contfrac(izd *z){
  INTERNAL_PRECISION *r, A = z -> A, *p = z -> a, *q = z -> ap;
  int dp = z -> dn, dq = z -> dd, type = z -> type;

  z -> db = 2 * dq + 1 + type;
  z -> beta = (INTERNAL_PRECISION *)
    malloc(z -> db * sizeof(INTERNAL_PRECISION));
  p = poly_factored_to_dense(A, p, dp);
  q = poly_factored_to_dense(ONE, q, dq);
  r = (INTERNAL_PRECISION *) malloc((MAX(dp,dq) + 1) *
				    sizeof(INTERNAL_PRECISION));
  if (type == 0) (void) contfrac_B(z -> beta, p, q, r, dp, dq);
  else (void) contfrac_A(z -> beta, p, q, r, dp, dq);
  free(p); free(q); free(r);
  return;
}

static INTERNAL_PRECISION *contfrac_A(INTERNAL_PRECISION *beta,
				      INTERNAL_PRECISION *p,
				      INTERNAL_PRECISION *q,
				      INTERNAL_PRECISION *r, int dp, int dq) {
  INTERNAL_PRECISION quot, *rb;
  int j;

  /* p(x) = x beta q(x) + r(x); dp = dq or dp = dq + 1 */

  quot = dp == dq ? ZERO : p[dp] / q[dq];
  r[0] = p[0];
  for (j = 1; j <= dp; j++) r[j] = p[j] - quot * q[j-1];
#ifdef DEBUG
  printf("%s: Continued Fraction form: deg p = %2d, deg q = %2d, beta = %g\n",
	 __FUNCTION__, dp, dq, (float) quot);
  for (j = 0; j <= dq + 1; j++)
    printf("\tp[%2d] = %14.6g\tq[%2d] = %14.6g\tr[%2d] = %14.6g\n",
	   j, (float) (j > dp ? ZERO : p[j]),
	   j, (float) (j == 0 ? ZERO : q[j-1]),
	   j, (float) (j == dp ? ZERO : r[j]));
#endif /* DEBUG */
  *(rb = contfrac_B(beta, q, r, p, dq, dq)) = quot;
  return rb + 1;
}

static INTERNAL_PRECISION *contfrac_B(INTERNAL_PRECISION *beta,
				      INTERNAL_PRECISION *p,
				      INTERNAL_PRECISION *q,
				      INTERNAL_PRECISION *r, int dp, int dq) {
  INTERNAL_PRECISION quot, *rb;
  int j;

  /* p(x) = beta q(x) + r(x); dp = dq or dp = dq - 1 */

  quot = dp == dq ? p[dp] / q[dq] : ZERO;
  for (j = 0; j < dq; j++) r[j] = p[j] - quot * q[j];
#ifdef DEBUG
  printf("%s: Continued Fraction form: deg p = %2d, deg q = %2d, beta = %g\n",
	 __FUNCTION__, dp, dq, (float) quot);
  for (j = 0; j <= dq; j++)
    printf("\tp[%2d] = %14.6g\tq[%2d] = %14.6g\tr[%2d] = %14.6g\n",
	   j, (float) (j > dp ? ZERO : p[j]),
	   j, (float) q[j],
	   j, (float) (j == dq ? ZERO : r[j]));
#endif /* DEBUG */
  *(rb = dq > 0 ? contfrac_A(beta, q, r, p, dq, dq-1) : beta) = quot;
  return rb + 1;
}

/* The global variable U is used to hold the argument u throughout the AGM
 * recursion. The global variables F and K are set in the innermost
 * instantiation of the recursive function AGM to the values of the elliptic
 * integrals F(u,k) and K(k) respectively. They must be made thread local to
 * make this code thread-safe in a multithreaded environment. */

static INTERNAL_PRECISION U, F, K;	/* THREAD LOCAL */

/* Recursive implementation of Gauss' arithmetico-geometric mean, which is the
 * kernel of the method used to compute the Jacobian elliptic functions
 * sn(u,k), cn(u,k), and dn(u,k) with parameter k (where 0 < k < 1), as well
 * as the elliptic integral F(s,k) satisfying F(sn(u,k)) = u and the complete
 * elliptic integral K(k).
 *
 * The algorithm used is a recursive implementation of the Gauss (Landen)
 * transformation.
 *
 * The function returns the value of sn(u,k'), where k' is the dual parameter,
 * and also sets the values of the global variables F and K.  The latter is
 * used to determine the sign of cn(u,k').
 *
 * The algorithm is deemed to have converged when b ceases to increase. This
 * works whatever INTERNAL_PRECISION is specified. */

static INTERNAL_PRECISION AGM(INTERNAL_PRECISION a,
			      INTERNAL_PRECISION b,
			      INTERNAL_PRECISION s) {
  static INTERNAL_PRECISION pb = -ONE;
  INTERNAL_PRECISION c, d, xi;

  if (b <= pb) {
    pb = -ONE;
    F = asin(s) / a;		/* Here, a is the AGM */
    K = M_PI / (TWO * a);
    return sin(U * a);
  }
  pb = b;
  c = a - b;
  d = a + b;
  xi = AGM(HALF*d, sqrt(a*b), ONE + c*c == ONE ?
	   HALF*s*d/a : (a - sqrt(a*a - s*s*c*d))/(c*s));
  return 2*a*xi / (d + c*xi*xi);
}

/* Computes sn(u,k), cn(u,k), dn(u,k), F(u,k), and K(k). It is essentially a
 * wrapper for the routine AGM. The sign of cn(u,k) is defined to be -1 if
 * K(k) < u < 3*K(k) and +1 otherwise, and thus sign is computed by some quite
 * unnecessarily obfuscated bit manipulations. */

static void sncndnFK(INTERNAL_PRECISION u, INTERNAL_PRECISION k,
		     INTERNAL_PRECISION* sn, INTERNAL_PRECISION* cn,
		     INTERNAL_PRECISION* dn, INTERNAL_PRECISION* elF,
		     INTERNAL_PRECISION* elK) {
  int sgn;
  U = u;
  *sn = AGM(ONE, sqrt(ONE - k*k), u);
  sgn = ((int) (fabs(u) / K)) % 4; /* sgn = 0, 1, 2, 3 */
  sgn ^= sgn >> 1;    /* (sgn & 1) = 0, 1, 1, 0 */
  sgn = 1 - ((sgn & 1) << 1);	/* sgn = 1, -1, -1, 1 */
  *cn = ((INTERNAL_PRECISION) sgn) * sqrt(ONE - *sn * *sn);
  *dn = sqrt(ONE - k*k* *sn * *sn);
  *elF = F;
  *elK = K;
}

/* Compute the coefficients for the optimal rational approximation R(x) to
 * sgn(x) of degree n over the interval epsilon < |x| < 1 using Zolotarev's
 * formula. 
 *
 * Set type = 0 for the Zolotarev approximation, which is zero at x = 0, and
 * type = 1 for the approximation which is infinite at x = 0. */

zolotarev_data* zolotarev(PRECISION epsilon, int n, int type) {
  INTERNAL_PRECISION A, c, cp, kp, ksq, sn, cn, dn, Kp, Kj, z, z0, t, M, F,
    l, invlambda, xi, xisq, *tv, s, opl;
  int m, czero, ts;
  zolotarev_data *zd;
  izd *d = (izd*) malloc(sizeof(izd));

  d -> type = type;
  d -> epsilon = (INTERNAL_PRECISION) epsilon;
  d -> n = n;
  d -> dd = n / 2;
  d -> dn = d -> dd - 1 + n % 2; /* n even: dn = dd - 1, n odd: dn = dd */
  d -> deg_denom = 2 * d -> dd;
  d -> deg_num = 2 * d -> dn + 1;

  d -> a = (INTERNAL_PRECISION*) malloc((1 + d -> dn) *
					sizeof(INTERNAL_PRECISION));
  d -> ap = (INTERNAL_PRECISION*) malloc(d -> dd *
					 sizeof(INTERNAL_PRECISION));
  ksq = d -> epsilon * d -> epsilon;
  kp = sqrt(ONE - ksq);
  sncndnFK(ZERO, kp, &sn, &cn, &dn, &F, &Kp); /* compute Kp = K(kp) */
  z0 = TWO * Kp / (INTERNAL_PRECISION) n;
  M = ONE;
  A = ONE / d -> epsilon;

  sncndnFK(HALF * z0, kp, &sn, &cn, &dn, &F, &Kj); /* compute xi */
  xi = ONE / dn;
  xisq = xi * xi;
  invlambda = xi;

  for (m = 0; m < d -> dd; m++) {
    czero = 2 * (m + 1) == n; /* n even and m = dd -1 */

    z = z0 * ((INTERNAL_PRECISION) m + ONE);
    sncndnFK(z, kp, &sn, &cn, &dn, &F, &Kj);
    t = cn / sn;
    c = - t*t;
    if (!czero) (d -> a)[d -> dn - 1 - m] = ksq / c;

    z = z0 * ((INTERNAL_PRECISION) m + HALF);
    sncndnFK(z, kp, &sn, &cn, &dn, &F, &Kj);
    t = cn / sn;
    cp = - t*t;
    (d -> ap)[d -> dd - 1 - m] = ksq / cp;

    M *= (ONE - c) / (ONE - cp);
    A *= (czero ? -ksq : c) * (ONE - cp) / (cp * (ONE - c));
    invlambda *= (ONE - c*xisq) / (ONE - cp*xisq);
  }
  invlambda /= M;
  d -> A = TWO / (ONE + invlambda) * A;
  d -> Delta = (invlambda - ONE) / (invlambda + ONE);

  d -> gamma = (INTERNAL_PRECISION*) malloc((1 + d -> n) *
					    sizeof(INTERNAL_PRECISION));
  l = ONE / invlambda;
  opl = ONE + l;
  sncndnFK(sqrt( d -> type == 1
		   ? (THREE + l) / (FOUR * opl)
		   : (ONE + THREE*l) / (opl*opl*opl)
	       ), sqrt(ONE - l*l), &sn, &cn, &dn, &F, &Kj);
  s = M * F;
  for (m = 0; m < d -> n; m++) {
    sncndnFK(s + TWO*Kp*m/n, kp, &sn, &cn, &dn, &F, &Kj);
    d -> gamma[m] = d -> epsilon / dn;
  }

  /* If R(x) is a Zolotarev rational approximation of degree (n,m) with maximum
   * error Delta, then (1 - Delta^2) / R(x) is also an optimal Chebyshev
   * approximation of degree (m,n) */

  if (d -> type == 1) {
    d -> A = (ONE - d -> Delta * d -> Delta) / d -> A;
    tv = d -> a; d -> a = d -> ap; d -> ap = tv;
    ts = d -> dn; d -> dn = d -> dd; d -> dd = ts;
    ts = d -> deg_num; d -> deg_num = d -> deg_denom; d -> deg_denom = ts;
  }

  construct_partfrac(d);
  construct_contfrac(d);

  /* Converting everything to PRECISION for external use only */

  zd = (zolotarev_data*) malloc(sizeof(zolotarev_data));
  zd -> A = (PRECISION) d -> A;
  zd -> Delta = (PRECISION) d -> Delta;
  zd -> epsilon = (PRECISION) d -> epsilon;
  zd -> n = d -> n;
  zd -> type = d -> type;
  zd -> dn = d -> dn;
  zd -> dd = d -> dd;
  zd -> da = d -> da;
  zd -> db = d -> db;
  zd -> deg_num = d -> deg_num;
  zd -> deg_denom = d -> deg_denom;

  zd -> a = (PRECISION*) malloc(zd -> dn * sizeof(PRECISION));
  for (m = 0; m < zd -> dn; m++) zd -> a[m] = (PRECISION) d -> a[m];
  free(d -> a);

  zd -> ap = (PRECISION*) malloc(zd -> dd * sizeof(PRECISION));
  for (m = 0; m < zd -> dd; m++) zd -> ap[m] = (PRECISION) d -> ap[m];
  free(d -> ap);

  zd -> alpha = (PRECISION*) malloc(zd -> da * sizeof(PRECISION));
  for (m = 0; m < zd -> da; m++) zd -> alpha[m] = (PRECISION) d -> alpha[m];
  free(d -> alpha);

  zd -> beta = (PRECISION*) malloc(zd -> db * sizeof(PRECISION));
  for (m = 0; m < zd -> db; m++) zd -> beta[m] = (PRECISION) d -> beta[m];
  free(d -> beta);

  zd -> gamma = (PRECISION*) malloc(zd -> n * sizeof(PRECISION));
  for (m = 0; m < zd -> n; m++) zd -> gamma[m] = (PRECISION) d -> gamma[m];
  free(d -> gamma);

  free(d);
  return zd;
}


void zolotarev_free(zolotarev_data *zdata)
{
    free(zdata -> a);
    free(zdata -> ap);
    free(zdata -> alpha);
    free(zdata -> beta);
    free(zdata -> gamma);
    free(zdata);
}


zolotarev_data* higham(PRECISION epsilon, int n) {
  INTERNAL_PRECISION A, M, c, cp, z, z0, t, epssq;
  int m, czero;
  zolotarev_data *zd;
  izd *d = (izd*) malloc(sizeof(izd));

  d -> type = 0;
  d -> epsilon = (INTERNAL_PRECISION) epsilon;
  d -> n = n;
  d -> dd = n / 2;
  d -> dn = d -> dd - 1 + n % 2; /* n even: dn = dd - 1, n odd: dn = dd */
  d -> deg_denom = 2 * d -> dd;
  d -> deg_num = 2 * d -> dn + 1;

  d -> a = (INTERNAL_PRECISION*) malloc((1 + d -> dn) *
					sizeof(INTERNAL_PRECISION));
  d -> ap = (INTERNAL_PRECISION*) malloc(d -> dd *
					 sizeof(INTERNAL_PRECISION));
  A = (INTERNAL_PRECISION) n;
  z0 = M_PI / A;
  A = n % 2 == 0 ? A : ONE / A;
  M = d -> epsilon * A;
  epssq = d -> epsilon * d -> epsilon;

  for (m = 0; m < d -> dd; m++) {
    czero = 2 * (m + 1) == n; /* n even and m = dd - 1*/

    if (!czero) {
      z = z0 * ((INTERNAL_PRECISION) m + ONE);
      t = tan(z);
      c = - t*t;
      (d -> a)[d -> dn - 1 - m] = c;
      M *= epssq - c;
    }

    z = z0 * ((INTERNAL_PRECISION) m + HALF);
    t = tan(z);
    cp = - t*t;
    (d -> ap)[d -> dd - 1 - m] = cp;
    M /= epssq - cp;
  }

  d -> gamma = (INTERNAL_PRECISION*) malloc((1 + d -> n) *
					    sizeof(INTERNAL_PRECISION));
  for (m = 0; m < d -> n; m++) d -> gamma[m] = ONE;

  d -> A = A;
  d -> Delta = ONE - M;

  construct_partfrac(d);
  construct_contfrac(d);

  /* Converting everything to PRECISION for external use only */

  zd = (zolotarev_data*) malloc(sizeof(zolotarev_data));
  zd -> A = (PRECISION) d -> A;
  zd -> Delta = (PRECISION) d -> Delta;
  zd -> epsilon = (PRECISION) d -> epsilon;
  zd -> n = d -> n;
  zd -> type = d -> type;
  zd -> dn = d -> dn;
  zd -> dd = d -> dd;
  zd -> da = d -> da;
  zd -> db = d -> db;
  zd -> deg_num = d -> deg_num;
  zd -> deg_denom = d -> deg_denom;

  zd -> a = (PRECISION*) malloc(zd -> dn * sizeof(PRECISION));
  for (m = 0; m < zd -> dn; m++) zd -> a[m] = (PRECISION) d -> a[m];
  free(d -> a);

  zd -> ap = (PRECISION*) malloc(zd -> dd * sizeof(PRECISION));
  for (m = 0; m < zd -> dd; m++) zd -> ap[m] = (PRECISION) d -> ap[m];
  free(d -> ap);

  zd -> alpha = (PRECISION*) malloc(zd -> da * sizeof(PRECISION));
  for (m = 0; m < zd -> da; m++) zd -> alpha[m] = (PRECISION) d -> alpha[m];
  free(d -> alpha);

  zd -> beta = (PRECISION*) malloc(zd -> db * sizeof(PRECISION));
  for (m = 0; m < zd -> db; m++) zd -> beta[m] = (PRECISION) d -> beta[m];
  free(d -> beta);

  zd -> gamma = (PRECISION*) malloc(zd -> n * sizeof(PRECISION));
  for (m = 0; m < zd -> n; m++) zd -> gamma[m] = (PRECISION) d -> gamma[m];
  free(d -> gamma);

  free(d);
  return zd;
}
}}

#ifdef TEST

#undef ZERO
#define ZERO ((PRECISION) 0)
#undef ONE
#define ONE ((PRECISION) 1)
#undef TWO
#define TWO ((PRECISION) 2)

/* Evaluate the rational approximation R(x) using the factored form */

static PRECISION zolotarev_eval(PRECISION x, zolotarev_data* rdata) {
  int m;
  PRECISION R;

  if (rdata -> type == 0) {
    R = rdata -> A * x;
    for (m = 0; m < rdata -> deg_denom/2; m++)
      R *= (2*(m+1) > rdata -> deg_num ? ONE : x*x - rdata -> a[m]) /
	(x*x - rdata -> ap[m]);
  } else {
    R = rdata -> A / x;
    for (m = 0; m < rdata -> deg_num/2; m++)
      R *= (x*x - rdata -> a[m]) /
	(2*(m+1) > rdata -> deg_denom ? ONE : x*x - rdata -> ap[m]);
  }
  return R;
}

/* Evaluate the rational approximation R(x) using the partial fraction form */

static PRECISION zolotarev_partfrac_eval(PRECISION x, zolotarev_data* rdata) {
  int m;
  PRECISION R = rdata -> alpha[rdata -> da - 1];
  for (m = 0; m < rdata -> dd; m++)
    R += rdata -> alpha[m] / (x * x - rdata -> ap[m]);
  if (rdata -> type == 1) R += rdata -> alpha[rdata -> dd] / (x * x);
  return R * x;
}    

/* Evaluate the rational approximation R(x) using continued fraction form. 
 *
 * If x = 0 and type = 1 then the result should be INF, whereas if x = 0 and
 * type = 0 then the result should be 0, but division by zero will occur at
 * intermediate stages of the evaluation. For IEEE implementations with
 * non-signalling overflow this will work correctly since 1/(1/0) = 1/INF = 0,
 * but with signalling overflow you will get an error message. */

static PRECISION zolotarev_contfrac_eval(PRECISION x, zolotarev_data* rdata) {
  int m;
  PRECISION R = rdata -> beta[0] * x;
  for (m = 1; m < rdata -> db; m++) R = rdata -> beta[m] * x + ONE / R;
  return R;
}    

/* Evaluate the rational approximation R(x) using Cayley form */

static PRECISION zolotarev_cayley_eval(PRECISION x, zolotarev_data* rdata) {
  int m;
  PRECISION T;

  T = rdata -> type == 0 ? ONE : -ONE;
  for (m = 0; m < rdata -> n; m++)
    T *= (rdata -> gamma[m] - x) / (rdata -> gamma[m] + x);
  return (ONE - T) / (ONE + T);
}

/* Test program. Apart from printing out the parameters for R(x) it produces
 * the following data files for plotting (unless NPLOT is defined):
 *
 * zolotarev-fn is a plot of R(x) for |x| < 1.2. This should look like sgn(x).
 *
 * zolotarev-err is a plot of the error |R(x) - sgn(x)| scaled by 1/Delta. This
 * should oscillate deg_num + deg_denom + 2 times between +1 and -1 over the
 * domain epsilon <= |x| <= 1.
 *
 * If ALLPLOTS is defined then zolotarev-partfrac (zolotarev-contfrac) is a
 * plot of the difference between the values of R(x) computed using the
 * factored form and the partial fraction (continued fraction) form, scaled by
 * 1/Delta. It should be zero everywhere. */

int main(int argc, char** argv) {
  
  int m, n, plotpts = 5000, type = 0;
  float eps, x, ypferr, ycferr, ycaylerr, maxypferr, maxycferr, maxycaylerr;
  zolotarev_data *rdata;
  PRECISION y;
  FILE *plot_function, *plot_error, 
    *plot_partfrac, *plot_contfrac, *plot_cayley;

  if (argc < 3 || argc > 4) {
    fprintf(stderr, "Usage: %s epsilon n [type]\n", *argv);
    exit(EXIT_FAILURE);
  }
  sscanf(argv[1], "%g", &eps);	/* First argument is epsilon */
  sscanf(argv[2], "%d", &n);	/* Second argument is n */
  if (argc == 4) sscanf(argv[3], "%d", &type); /* Third argument is type */

  if (type < 0 || type > 2) {
    fprintf(stderr, "%s: type must be 0 (Zolotarev R(0) = 0),\n"
	    "\t\t1 (Zolotarev R(0) = Inf, or 2 (Higham)\n", *argv);
    exit(EXIT_FAILURE);
  }

  rdata = type == 2 
    ? higham((PRECISION) eps, n) 
    : zolotarev((PRECISION) eps, n, type);

  printf("Zolotarev Test: R(epsilon = %g, n = %d, type = %d)\n\t" 
	 STRINGIFY(VERSION) "\n\t" STRINGIFY(HVERSION)
	 "\n\tINTERNAL_PRECISION = " STRINGIFY(INTERNAL_PRECISION)
	 "\tPRECISION = " STRINGIFY(PRECISION)
	 "\n\n\tRational approximation of degree (%d,%d), %s at x = 0\n"
	 "\tDelta = %g (maximum error)\n\n"
	 "\tA = %g (overall factor)\n",
	 (float) rdata -> epsilon, rdata -> n, type,
	 rdata -> deg_num, rdata -> deg_denom,
	 rdata -> type == 1 ? "infinite" : "zero",
	 (float) rdata -> Delta, (float) rdata -> A);
  for (m = 0; m < MIN(rdata -> dd, rdata -> dn); m++)
    printf("\ta[%2d] = %14.8g\t\ta'[%2d] = %14.8g\n",
	   m + 1, (float) rdata -> a[m], m + 1, (float) rdata -> ap[m]);
  if (rdata -> dd > rdata -> dn)
    printf("\t\t\t\t\ta'[%2d] = %14.8g\n",
	   rdata -> dn + 1, (float) rdata -> ap[rdata -> dn]);
  if (rdata -> dd < rdata -> dn)
    printf("\ta[%2d] = %14.8g\n",
	   rdata -> dd + 1, (float) rdata -> a[rdata -> dd]);

  printf("\n\tPartial fraction coefficients\n");
  printf("\talpha[ 0] = %14.8g\n",
	 (float) rdata -> alpha[rdata -> da - 1]);
  for (m = 0; m < rdata -> dd; m++)
    printf("\talpha[%2d] = %14.8g\ta'[%2d] = %14.8g\n",
	   m + 1, (float) rdata -> alpha[m], m + 1, (float) rdata -> ap[m]);
  if (rdata -> type == 1)
    printf("\talpha[%2d] = %14.8g\ta'[%2d] = %14.8g\n",
	   rdata -> dd + 1, (float) rdata -> alpha[rdata -> dd],
	   rdata -> dd + 1, (float) ZERO);

  printf("\n\tContinued fraction coefficients\n");
  for (m = 0; m < rdata -> db; m++)
    printf("\tbeta[%2d] = %14.8g\n", m, (float) rdata -> beta[m]);

  printf("\n\tCayley transform coefficients\n");
  for (m = 0; m < rdata -> n; m++)
    printf("\tgamma[%2d] = %14.8g\n", m, (float) rdata -> gamma[m]);

#ifndef NPLOT
  plot_function = fopen("zolotarev-fn.dat", "w");
  plot_error = fopen("zolotarev-err.dat", "w");
#ifdef ALLPLOTS
  plot_partfrac = fopen("zolotarev-partfrac.dat", "w");
  plot_contfrac = fopen("zolotarev-contfrac.dat", "w");
  plot_cayley = fopen("zolotarev-cayley.dat", "w");
#endif /* ALLPLOTS */
  for (m = 0, maxypferr = maxycferr = maxycaylerr = 0.0; m <= plotpts; m++) {
    x = 2.4 * (float) m / plotpts - 1.2;
    if (rdata -> type == 0 || fabs(x) * (float) plotpts > 1.0) {
      /* skip x = 0 for type 1, as R(0) is singular */
      y = zolotarev_eval((PRECISION) x, rdata);
      fprintf(plot_function, "%g %g\n", x, (float) y);
      fprintf(plot_error, "%g %g\n",
	      x, (float)((y - ((x > 0.0 ? ONE : -ONE))) / rdata -> Delta));
      ypferr = (float)((zolotarev_partfrac_eval((PRECISION) x, rdata) - y)
		       / rdata -> Delta);
      ycferr = (float)((zolotarev_contfrac_eval((PRECISION) x, rdata) - y)
		       / rdata -> Delta);
      ycaylerr = (float)((zolotarev_cayley_eval((PRECISION) x, rdata) - y)
		       / rdata -> Delta);
      if (fabs(x) < 1.0 && fabs(x) > rdata -> epsilon) {
	maxypferr = MAX(maxypferr, fabs(ypferr));
	maxycferr = MAX(maxycferr, fabs(ycferr));
	maxycaylerr = MAX(maxycaylerr, fabs(ycaylerr));
      }
#ifdef ALLPLOTS
      fprintf(plot_partfrac, "%g %g\n", x, ypferr);
      fprintf(plot_contfrac, "%g %g\n", x, ycferr);
      fprintf(plot_cayley, "%g %g\n", x, ycaylerr);
#endif /* ALLPLOTS */
    }
  }
#ifdef ALLPLOTS
  fclose(plot_cayley);
  fclose(plot_contfrac);
  fclose(plot_partfrac);
#endif /* ALLPLOTS */
  fclose(plot_error);
  fclose(plot_function);

  printf("\n\tMaximum PF error = %g (relative to Delta)\n", maxypferr);
  printf("\tMaximum CF error = %g (relative to Delta)\n", maxycferr);
  printf("\tMaximum Cayley error = %g (relative to Delta)\n", maxycaylerr);
#endif /* NPLOT */

  free(rdata -> a);
  free(rdata -> ap);
  free(rdata -> alpha);
  free(rdata -> beta);
  free(rdata -> gamma);
  free(rdata);

  return EXIT_SUCCESS;
}


#endif /* TEST */
