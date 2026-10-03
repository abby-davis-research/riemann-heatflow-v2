#!/usr/bin/env python3
"""
ht_flow.py — fixed-node tanh-sinh quadrature core for de Bruijn-Newman heat flows.

H_t(z) = sum_j c_j * exp(t*u_j^2) * exp(i*u_j*z),   c_j = w_j * Phi(u_j)

The complex-exponential form is the analytic continuation; for real z the
imaginary part is ~1e-50 (rounding). All derivatives analytic in z.
"""
import os
import mpmath
from mpmath import mp


def _gl_compute(J, dps):
    """Gauss-Legendre nodes/weights on [-1,1] at dps digits. numpy seed + mp Newton."""
    import numpy as np
    mp.dps = dps
    xs0, ws0 = np.polynomial.legendre.leggauss(J)
    xs, ws = [], []
    for x0 in xs0:
        x = mp.mpf(str(x0))
        for _ in range(12):
            p = mp.legendre(J, x)
            pm1 = mp.legendre(J - 1, x)
            dp = J / (x * x - 1) * (x * p - pm1)
            dx = p / dp
            x -= dx
            if abs(dx) < mp.mpf(10) ** (-(dps - 5)):
                break
        p = mp.legendre(J, x)
        pm1 = mp.legendre(J - 1, x)
        dp = J / (x * x - 1) * (x * p - pm1)
        xs.append(x)
        ws.append(2 / ((1 - x * x) * dp * dp))
    return xs, ws


def gl_grid(J=256, umax=1.0, dps=50, cache_dir='/tmp'):
    """Fixed Gauss-Legendre nodes/weights on [0, umax], cached on disk."""
    import os
    mp.dps = dps
    fn = os.path.join(cache_dir, 'gl_%d_umax%s_%d.json' % (J, str(umax), dps))
    if os.path.exists(fn):
        import json
        d = json.load(open(fn))
        us = [mp.mpf(s) for s in d['us']]
        ws = [mp.mpf(s) for s in d['ws']]
        return us, ws
    xs, ws0 = _gl_compute(J, dps)
    umax = mp.mpf(umax)
    us = [umax / 2 * (x + 1) for x in xs]
    ws = [w * umax / 2 for w in ws0]
    import json
    json.dump({'us': [str(u) for u in us], 'ws': [str(w) for w in ws]},
              open(fn, 'w'))
    return us, ws


def tanh_sinh_grid(K=160, h=0.25, umax=2.0, dps=50):
    """Fixed tanh-sinh nodes/weights on [0, umax]. Returns (us, ws) as mpf lists."""
    mp.dps = dps
    h = mp.mpf(h)
    umax = mp.mpf(umax)
    us, ws = [], []
    half_pi = mp.pi / 2
    for k in range(-K, K + 1):
        kh = h * k
        sinh_kh = mp.sinh(kh)
        x = mp.tanh(half_pi * sinh_kh)          # in (-1, 1)
        w = half_pi * mp.cosh(kh) / (mp.cosh(half_pi * sinh_kh) ** 2) * h
        u = umax / 2 * (1 + x)
        us.append(u)
        ws.append(umax / 2 * w)
    return us, ws


def phi_zeta(u):
    """Riemann Phi(u) = sum_{n>=1} (2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) e^{-pi n^2 e^{4u}}."""
    s = mp.mpf(0)
    e4u = mp.e ** (4 * u)
    for n in range(1, 9):
        n2 = n * n
        s += (2 * mp.pi ** 2 * n2 * n2 * mp.e ** (9 * u)
              - 3 * mp.pi * n2 * mp.e ** (5 * u)) * mp.e ** (-mp.pi * n2 * e4u)
    return s


def dh_kappa(dps=50):
    mp.dps = dps
    return (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)


def phi_dh_wide(v, kappa):
    """Wide-frame DH Phi: 8 e^{3v} sum_{n>=1} n a_n e^{-pi n^2 e^{4v}/5},
    a = (1, kappa, -kappa, -1, 0) period 5."""
    a = [mp.mpf(1), kappa, -kappa, mp.mpf(-1), mp.mpf(0)]
    s = mp.mpf(0)
    e4v = mp.e ** (4 * v)
    for n in range(1, 31):
        an = a[(n - 1) % 5]
        if an:
            s += n * an * mp.e ** (-mp.pi * n * n * e4v / 5)
    return 8 * mp.e ** (3 * v) * s


class HeatFlow:
    """H_t(z) and derivatives via precomputed quadrature. phi: u -> Phi(u)."""

    def __init__(self, phi, J=256, umax=1.0, dps=50, grid='gl'):
        mp.dps = dps
        self.dps = dps
        if grid == 'gl':
            us, ws = gl_grid(J=J, umax=umax, dps=dps)
        else:
            us, ws = tanh_sinh_grid(K=J, umax=umax, dps=dps)
        self.u = us
        self.u2 = [x * x for x in us]
        # drop negligible weights to speed up (keep 1e-60 relative)
        self.c = []
        self.uu = []
        self.u2u = []
        for w, u, u2 in zip(ws, us, self.u2):
            c = w * phi(u)
            if abs(c) > mp.mpf(10) ** (-(dps - 10)):
                self.c.append(c)
                self.uu.append(u)
                self.u2u.append(u2)
        self.J = len(self.c)

    def H(self, z, t):
        z = mp.mpc(z)
        t = mp.mpf(t)
        s = mp.mpc(0)
        for c, u, u2 in zip(self.c, self.uu, self.u2u):
            s += c * mp.e ** (t * u2) * mp.cos(u * z)
        return s

    def Hp(self, z, t):
        z = mp.mpc(z)
        t = mp.mpf(t)
        s = mp.mpc(0)
        for c, u, u2 in zip(self.c, self.uu, self.u2u):
            s += c * mp.e ** (t * u2) * (-u) * mp.sin(u * z)
        return s

    def Hpp(self, z, t):
        z = mp.mpc(z)
        t = mp.mpf(t)
        s = mp.mpc(0)
        for c, u, u2 in zip(self.c, self.uu, self.u2u):
            s += c * mp.e ** (t * u2) * (-u * u) * mp.cos(u * z)
        return s

    def H0_derivs(self, z, K):
        """[H_0^{(k)}(z)] for k=0..K via the quadrature (t=0)."""
        z = mp.mpc(z)
        out = [mp.mpc(0)] * (K + 1)
        for c, u, u2 in zip(self.c, self.uu, self.u2u):
            uz = u * z
            co = mp.cos(uz); si = mp.sin(uz)
            pw = mp.mpc(1)  # u^k
            for k in range(K + 1):
                r = k % 4
                if r == 0:
                    out[k] += c * pw * co
                elif r == 1:
                    out[k] += c * pw * (-si)
                elif r == 2:
                    out[k] += c * pw * (-co)
                else:
                    out[k] += c * pw * si
                pw *= u
        return out

    def pde_residual(self, z, t, dt=1e-6):
        """|dH/dt + d2H/dz2| with dH/dt by central difference."""
        dHdt = (self.H(z, t + dt) - self.H(z, t - dt)) / (2 * dt)
        return abs(dHdt + self.Hpp(z, t))


# ---- Davenport-Heilbronn via Hurwitz (verification only) ----

def dh_f(s, kappa):
    """f(s) = 5^{-s} sum_{a=1..4} b(a) zeta_H(s, a/5), b=(1,k,-k,-1)."""
    b = [mp.mpf(1), kappa, -kappa, mp.mpf(-1)]
    return mp.power(5, -s) * sum(
        b[a - 1] * mp.hurwitz(s, mp.mpf(a) / 5) for a in range(1, 5))


def dh_Lambda(s, kappa):
    return (mp.mpf(5) / mp.pi) ** ((s + 1) / 2) * mp.gamma((s + 1) / 2) * dh_f(s, kappa)
