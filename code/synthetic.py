#!/usr/bin/env python3
"""
synthetic.py — RH-false synthetic flows via R(z) surgery (prereg v2.0, §4).

    H_0^syn(z) = H_0(z) * R(z),
    R(z) = [(1 - z^2/w^2)(1 - z^2/wbar^2)] / (1 - z^2/z_1^2),  w = z_1 + i*delta.

Surgically replaces the real pair ±z_1 with the quadruplet ±w, ±w̄; R is real
on R, R(0)=1, R→1 at ∞, and H_0*R is entire (denominator pole cancels against
H_0's simple zero). N_off(0) = 4 by construction.

H_t^syn(z) = e^{-t d^2/dz^2}[H_0*R] via Taylor-in-t with the Leibniz rule,
truncation M=48 (prereg §4; gated against M=96 in V4).

CORRECTION vs v1: v1's _assemble looped j over [0, min(k,4)], dropping the
h^{(j)} terms with j in (4, k] for k > 4 — i.e. the Taylor coefficients
F^{(2m)} for m >= 3 were wrong. Here j ranges over [max(0,k-4), k], the full
Leibniz range for N of degree 4, and h is computed to order 2M.

Performance: h-derivatives and the F array are t-independent, so they are
computed once per z and cached; H/Hp/Hpp at any t are then cheap sums.
"""
from mpmath import mp


class SyntheticFlow:
    def __init__(self, base_flow, z1, delta, M=48):
        self.base = base_flow
        self.z1 = mp.mpf(z1)
        self.delta = mp.mpf(delta)
        self.M = M
        self.w = self.z1 + mp.j * self.delta
        self.C0 = -self.z1 ** 2
        a = 1 / self.w ** 2
        b = 1 / mp.conj(self.w) ** 2
        self._A = a + b          # N(z) = 1 - A z^2 + B z^4
        self._B = a * b
        self._Fcache = {}        # str(z) -> [F^{(0)}(z), ..., F^{(2M)}(z)]

    # ---- N(z) and derivatives (deg 4 polynomial; exact) ----
    def _N_deriv(self, z, k):
        A, B = self._A, self._B
        if k == 0:
            return 1 - A * z ** 2 + B * z ** 4
        if k == 1:
            return -2 * A * z + 4 * B * z ** 3
        if k == 2:
            return -2 * A + 12 * B * z ** 2
        if k == 3:
            return 24 * B * z
        if k == 4:
            return 24 * B
        return mp.mpc(0)

    # ---- h(z) = H_0(z)/(z^2 - z_1^2), entire ----
    def _h_derivs(self, z, K):
        """h^{(0..K)}(z) via the recurrence from h*D = H_0, D = z^2 - z_1^2:
        h^{(k)} = [H_0^{(k)} - 2k z h^{(k-1)} - k(k-1) h^{(k-2)}] / D(z)."""
        H0d = self.base.H0_derivs(z, K)
        D = z * z - self.z1 * self.z1
        h = [mp.mpc(0)] * (K + 1)
        h[0] = H0d[0] / D
        if K >= 1:
            h[1] = (H0d[1] - 2 * z * h[0]) / D
        for k in range(2, K + 1):
            h[k] = (H0d[k] - 2 * k * z * h[k - 1]
                    - k * (k - 1) * h[k - 2]) / D
        return h

    def _F_at(self, z):
        """[F^{(k)}(z)] for k = 0..2M, F = C0 * h * N. Cached per z."""
        key = str(z)
        hit = self._Fcache.get(key)
        if hit is not None:
            return hit
        M = self.M
        K = 2 * M
        h = self._h_derivs(z, K)
        F = []
        for k in range(K + 1):
            s = mp.mpc(0)
            # N has degree 4 => N^{(k-j)} != 0 iff k-j <= 4 iff j >= k-4
            for j in range(max(0, k - 4), k + 1):
                s += mp.binomial(k, j) * h[j] * self._N_deriv(z, k - j)
            F.append(self.C0 * s)
        self._Fcache[key] = F
        return F

    def _taylor(self, z, t, deriv):
        """Sum_{m} (-t)^m/m! * F^{(2m+deriv)}(z), deriv = 0,1,2 for H,Hp,Hpp."""
        F = self._F_at(z)
        M = self.M
        t = mp.mpf(t)
        s = mp.mpc(0)
        c = mp.mpf(1)  # (-t)^m / m!, updated iteratively
        for m in range(M + 1):
            k = 2 * m + deriv
            if k <= 2 * M:
                s += c * F[k]
            c *= -t / (m + 1)
        return s

    def H(self, z, t):
        return self._taylor(z, t, 0)

    def Hp(self, z, t):
        return self._taylor(z, t, 1)

    def Hpp(self, z, t):
        return self._taylor(z, t, 2)

    def pde_residual(self, z, t, dt=1e-6):
        """|dH/dt + d2H/dz2| with dH/dt by central difference (for G-c, V4)."""
        dHdt = (self.H(z, t + dt) - self.H(z, t - dt)) / (2 * dt)
        return abs(dHdt + self.Hpp(z, t))

    def construction_check(self):
        """|H_0^syn(w)| — must be < 1e-25 (V4 construction check, prereg §7)."""
        return abs(self.H(self.w, 0))

    def n_off_zero(self):
        """N_off(0) by construction: the injected quadruplet."""
        return 4
