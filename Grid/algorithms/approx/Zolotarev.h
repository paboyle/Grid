/* -*- Mode: C; comment-column: 22; fill-column: 79; -*- */

#ifdef __cplusplus
namespace Grid {
namespace Approx {
#endif

#define HVERSION Header Time-stamp: <14-OCT-2004 09:26:51.00 adk@MISSCONTRARY>


#ifndef ZOLOTAREV_INTERNAL
#ifndef PRECISION
#define PRECISION double
#endif
#define ZPRECISION PRECISION
#define ZOLOTAREV_DATA zolotarev_data
#endif

/* This struct contains the coefficients which parameterise an optimal rational
 * approximation to the signum function.
 *
 * The parameterisations used are:
 *
 * Factored form for type 0 (R0(0) = 0)
 *
 * R0(x) = A * x * prod(x^2 - a[j], j = 0 .. dn-1) / prod(x^2 - ap[j], j = 0
 * .. dd-1),
 *
 * where deg_num = 2*dn + 1 and deg_denom = 2*dd.
 *
 * Factored form for type 1 (R1(0) = infinity)
 *
 * R1(x) = (A / x) * prod(x^2 - a[j], j = 0 .. dn-1) / prod(x^2 - ap[j], j = 0
 * .. dd-1),
 *
 * where deg_num = 2*dn and deg_denom = 2*dd + 1. 
 *
 * Partial fraction form
 *
 * R(x) = alpha[da] * x + sum(alpha[j] * x / (x^2 - ap[j]), j = 0 .. da-1)
 *
 * where da = dd for type 0 and da = dd + 1 with ap[dd] = 0 for type 1.
 *
 * Continued fraction form 
 *
 * R(x) = beta[db-1] * x + 1 / (beta[db-2] * x + 1 / (beta[db-3] * x + ...))
 *
 * with the final coefficient being beta[0], with d' = 2 * dd + 1 for type 0
 * and db = 2 * dd + 2 for type 1.
 *
 * Cayley form (Chiu's domain wall formulation)
 *
 * R(x) = (1 - T(x)) / (1 + T(x))
 *
 * where T(x) = prod((x - gamma[j]) / (x + gamma[j]), j = 0 .. n-1)
 */

typedef struct {
  ZPRECISION *a,      /* zeros of numerator, a[0 .. dn-1] */
    *ap,	      /* poles (zeros of denominator), ap[0 .. dd-1] */
    A,		      /* overall factor */
    *alpha,	      /* coefficients of partial fraction, alpha[0 .. da-1] */
    *beta,	      /* coefficients of continued fraction, beta[0 .. db-1] */
    *gamma,	      /* zeros of numerator of T in Cayley form */
    Delta,	      /* maximum error, |R(x) - sgn(x)| <= Delta */
    epsilon;	      /* minimum x value, epsilon < |x| < 1 */
  int n,	      /* approximation degree */
    type,	      /* 0: R(0) = 0, 1: R(0) = infinity */
    dn, dd, da, db,   /* number of elements of a, ap, alpha, and beta */
    deg_num,	      /* degree of numerator = deg_denom +/- 1 */
    deg_denom;	      /* degree of denominator */
} ZOLOTAREV_DATA;

#ifndef ZOLOTAREV_INTERNAL

/* zolotarev(epsilon, n, type) returns a pointer to an initialised
 * zolotarev_data structure. The arguments must satisfy the constraints that
 * epsilon > 0, n > 0, and type = 0 or 1. */

ZOLOTAREV_DATA* higham(PRECISION epsilon, int n) ;
ZOLOTAREV_DATA* zolotarev(PRECISION epsilon, int n, int type);
void zolotarev_free(zolotarev_data *zdata);
#endif

#ifdef __cplusplus
}}
#endif
