#!/usr/bin/env python3
"""
Reconstruct the GCR smoother polynomial from PGCR coefficient logs
(SmootherCoeffLog=1 in Example_pvdagm_3level_DenseCoarseMatrix, or
LogCoefficients(1) on any PrecGeneralisedConjugateResidualNonHermitian).

GCR with z = r (trivial preconditioner) builds  p_k = P_k(A) r0,
r_k = R_k(A) r0:

    R_0 = P_0 = 1
    R_{k+1}(l) = R_k(l) - a_k l P_k(l)
    P_{k+1}(l) = R_{k+1}(l) + sum_j b_{k,j} P_{k-j}(l)

so after m steps  x_m = S(A) r0  with  S = sum_k a_k P_k  and residual
polynomial  R_m = 1 - l S(l),  R_m(0) = 1.

Usage
  gcr_polynomial.py LOG [--name Fsmoother] [--lmin 0] [--lmax 3] [--n 600]
      writes <name>_poly.dat : l  |R_mean(l)|  Re R_mean  Im R_mean  |R|min  |R|max  |S_mean(l)|
      where mean = polynomial built from the per-step MEAN coefficients over
      all calls, and min/max are the envelope over the individual calls.
      Prints the mean coefficients and their relative spread per step.
  gcr_polynomial.py --selftest
      runs GCR in numpy on a random non-normal matrix, logs a,b the same way,
      and checks  R_m(A) r0 == the actual residual  to rounding.

gnuplot:
  set logscale y; plot 'Fsmoother_poly.dat' u 1:2 w l t '|R_m| mean', \
                       '' u 1:5 w l t 'min', '' u 1:6 w l t 'max'
"""
import re, sys, argparse
import numpy as np
from numpy.polynomial import polynomial as Poly   # coefficient arrays, low->high

A_RE = re.compile(r'(\S+)\s+coeff\[(\d+)\]\s+a=\(([^,]+),([^)]+)\)')
B_RE = re.compile(r'(\S+)\s+coeff\[(\d+)\]((?:\s+b\[\d+\]=\([^)]*\))+)')
BJ_RE = re.compile(r'b\[(\d+)\]=\(([^,]+),([^)]+)\)')

def parse(path, name):
    """-> list of calls; each call = list of (a_k, [b_k0, b_k1, ...]) in step order."""
    calls, cur = [], None
    for line in open(path, errors='replace'):
        m = A_RE.search(line)
        if m and m.group(1) == name:
            k = int(m.group(2)); a = complex(float(m.group(3)), float(m.group(4)))
            if k == 0:
                if cur: calls.append(cur)
                cur = []
            if cur is not None:
                cur.append([a, []])
            continue
        m = B_RE.search(line)
        if m and m.group(1) == name and cur is not None:
            k = int(m.group(2))
            bs = [complex(float(x[1]), float(x[2])) for x in BJ_RE.findall(m.group(3))]
            if k < len(cur): cur[k][1] = bs
    if cur: calls.append(cur)
    return calls

def polynomials(coeffs):
    """coeffs: list of (a_k, [b_kj]) -> (S, R) as low->high complex coefficient arrays."""
    R = [np.array([1.0+0j])]; P = [np.array([1.0+0j])]
    S = np.array([0.0+0j])
    for k, (a, bs) in enumerate(coeffs):
        S = Poly.polyadd(S, a*P[k])
        Rn = Poly.polysub(R[k], a*Poly.polymulx(P[k]))          # R_{k+1} = R_k - a l P_k
        Pn = Rn.copy()
        for j, b in enumerate(bs):                              # P_{k+1} = R_{k+1} + sum b_kj P_{k-j}
            Pn = Poly.polyadd(Pn, b*P[k-j])
        R.append(Rn); P.append(Pn)
    return S, R[-1]

def mean_coeffs(calls):
    m = min(len(c) for c in calls)
    out, spread = [], []
    for k in range(m):
        a = np.array([c[k][0] for c in calls])
        nb = min(len(c[k][1]) for c in calls)
        bs = [np.array([c[k][1][j] for c in calls]) for j in range(nb)]
        out.append((a.mean(), [b.mean() for b in bs]))
        sa = a.std()/max(abs(a.mean()),1e-300)
        sb = [b.std()/max(abs(b.mean()),1e-300) for b in bs]
        spread.append((sa, sb))
    return out, spread

def selftest(n=40, m=6, mmax=6, seed=1):
    rng = np.random.default_rng(seed)
    A = np.eye(n)*1.5 + 0.4*rng.standard_normal((n,n)) + 0.2j*rng.standard_normal((n,n))  # non-normal
    r0 = rng.standard_normal(n) + 1j*rng.standard_normal(n)
    r = r0.copy(); x = np.zeros(n, complex)
    p = [r.copy()]; q = [A@r]; qq = [np.vdot(q[0],q[0]).real]
    log = []
    for k in range(m):
        a = np.vdot(q[k], r)/qq[k]                               # <q,r>/<q,q>
        x = x + a*p[k]; r = r - a*q[k]
        z = r; Az = A@z; pn = z.copy(); qn = Az.copy(); bs = []
        for back in range(min(k+1, mmax-1)):
            b = -(np.vdot(q[k-back], Az)/qq[k-back]).real        # real part, as Grid does
            pn = pn + b*p[k-back]; qn = qn + b*q[k-back]; bs.append(complex(b))
        p.append(pn); q.append(qn); qq.append(np.vdot(qn,qn).real)
        log.append((a, bs))
    S, R = polynomials(log)
    def apply(poly, v):
        out = np.zeros_like(v); Ak = v.copy()
        for c in poly: out = out + c*Ak; Ak = A@Ak
        return out
    err_r = np.linalg.norm(apply(R, r0) - r)/np.linalg.norm(r)
    err_x = np.linalg.norm(apply(S, r0) - x)/np.linalg.norm(x)
    print(f"selftest: |R_m(A)r0 - r_m|/|r_m| = {err_r:.2e}   |S(A)r0 - x_m|/|x_m| = {err_x:.2e}")
    return err_r < 1e-10 and err_x < 1e-10

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("log", nargs="?")
    ap.add_argument("--name", default="Fsmoother")
    ap.add_argument("--lmin", type=float, default=0.0)
    ap.add_argument("--lmax", type=float, default=3.0)
    ap.add_argument("--n", type=int, default=600)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        sys.exit(0 if selftest() else 1)
    calls = parse(args.log, args.name)
    if not calls:
        sys.exit(f"no '{args.name} coeff[...]' lines found in {args.log}")
    mean, spread = mean_coeffs(calls)
    print(f"{args.name}: {len(calls)} calls, {len(mean)} steps")
    for k,(a,bs) in enumerate(mean):
        print(f"  step {k}: a = ({a.real:+.5f},{a.imag:+.5f})  rel spread {spread[k][0]:.3f}   "
              + "  ".join(f"b[{j}]={b.real:+.5f} (spread {spread[k][1][j]:.3f})" for j,b in enumerate(bs)))
    S, R = polynomials(mean)
    print("  S coefficients (low->high):", np.array2string(S, precision=5))
    lam = np.linspace(args.lmin, args.lmax, args.n)
    Rm = Poly.polyval(lam, R); Sm = Poly.polyval(lam, S)
    per = np.array([np.abs(Poly.polyval(lam, polynomials(c)[1])) for c in calls])
    out = f"{args.name}_poly.dat"
    with open(out, "w") as f:
        f.write("# lambda |R_mean| ReR ImR |R|min |R|max |S_mean|\n")
        for i,l in enumerate(lam):
            f.write(f"{l:.6f} {abs(Rm[i]):.6e} {Rm[i].real:.6e} {Rm[i].imag:.6e} "
                    f"{per[:,i].min():.6e} {per[:,i].max():.6e} {abs(Sm[i]):.6e}\n")
    print(f"  wrote {out}")
