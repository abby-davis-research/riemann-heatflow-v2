# Preregistration Amendment v2.1 — Riemann Heat-Flow Instrument Repair

**Amends:** Preregistration v2.0, DOI `10.5281/zenodo.23120530`
(deposited 2026-10-03, before any v2 `t≠0` computation).

**Author:** Abby Davis (independent researcher, Tucson AZ).

**Date:** 2026-10-03.

**Status:** This amendment is timestamped on Zenodo **before any v2.1
`t≠0` computation**. No v2.1 tracking, event detection, or Phase-1/Phase-2
computation has been run. All `t≠0` results reported herein are from v2.0
gate diagnostics (already red-flagged) or from `t=0` instrument
characterization, which is not barred.

---

## §1. Why v2.1 exists: the v2.0 gate outcome

On 2026-10-03 the v2.0 gate suite ran. Outcome:

| Gate | Result | Measured vs threshold |
|------|--------|----------------------|
| V1 (J=256→512) | **FAIL** | 3.79e-22 at z=275 vs < 1e-35 |
| V2 (first five zeros) | PASS | — |
| V3 (PDE residual) | **FAIL** | 1.25e-19 vs < 1e-20 |
| V4 (synthetic) | **FAIL** (marginal) | PDE sub-check 1.01e-20 vs < 1e-20 |
| V6 (null-window detector) | PASS | — |
| CAL-(i) (N_off(0)=0) | PASS | — |
| DH case (G1/G2/G3) | **DROPPED** | G3 spread 3.2e-07 vs < 1e-25 |

Additionally, t=0 seeding found only **60 of ~104** expected zeros in
`|z| ≤ 300`, and Newton polishing converged to roots displaced up to ~5
units from known references in the band `159 < |z| < 258`.

Per prereg v2.0 §7, Phase 1 and Phase 2 never ran. No event-map science
data was collected. All failures were red-flagged in the lab notebook.
This amendment repairs the instrument; it does not rescue, reinterpret,
or narrow any v2.0 result.

## §2. Root-cause diagnosis (all verified by direct measurement)

**Cause A — arithmetic precision (dps=50).** At `dps=50`, the quadrature
suffers ~35-digit cancellation in z-dependent high-z bands: `H_0` goes
"ghost-flat" (`|H'| ~ 1e-26`), and Newton cannot polish. Raising to
`dps=80` alone did **not** repair seeding (100/104 found, 22 roots
displaced > 1 unit) — precision was necessary but not sufficient.

**Cause B — quadrature truncation (umax=3.0).** With `umax=3.0`, the
fixed Gauss–Legendre nodes are too sparse to resolve `cos(zu)` at high
`z`: measured J=128→256 differences reached 2e-03 at z=300. The Φ mass
is concentrated in `u ∈ [0,1]` (Φ(1) ~ 1e-74); spreading nodes over
`[0,3]` wastes resolution where the integrand oscillates fastest.

**Cause C — the PDE check measured itself.** The v2.0 `pde_residual`
used a central finite difference with `dt=1e-6`. Its `O(dt²)` truncation
is 1.25e-19 — larger than the v2.0 V3 threshold (1e-20) and ~23 orders
above the true PDE residual (~1e-42 estimated). V3 and the V4 PDE
sub-check were therefore testing finite-difference truncation, not the
heat equation. The failures were test artifacts, not instrument failures.

**Cause D (DH only) — separate issue.** The DH G3 spread (3.2e-07) is a
DH-quadrature matter distinct from Causes A–C. It is out of scope for
this repair; the DH case remains DROPPED (§4).

## §3. The v2.1 repaired instrument

| Parameter | v2.0 | v2.1 | Why |
|-----------|------|------|-----|
| `dps` | 50 | **80** | Cause A: rounding now at 1e-75 |
| `umax` | 3.0 | **1.2** | Cause B: nodes where Φ has mass |
| `J` (production / refinement) | 256 / 512 | 256 / 512 | unchanged |
| `dH/dt` in PDE checks | finite difference | **analytic** (`Ht`) | Cause C |

Measured performance of the repaired instrument (all at `t=0` unless
noted, `dps=80`, `J=256`, `umax=1.2`):

- J=256→512 difference: **~2e-71, uniform in z** for `|z| ≤ 300`
  (was 3.79e-22 at z=275).
- True PDE residual (analytic `dH/dt`): **~1e-85** at all 10 test
  `(z,t)` (was measured as 1.25e-19 by the finite-difference artifact).
- Seeding: **104/104** expected zeros recovered in `|z| ≤ 300`
  (was 60/104); worst `|H_0| ~ 1e-70`; worst displacement from known
  `±2γ_n` < 0.01 (was ~5).

## §4. The v2.1 gate suite (revised)

All v2.0 gates are re-run from scratch under v2.1. Changes:

- **V1:** threshold **1e-60** (was 1e-35; miscalibrated for the wide
  window). Measured 2e-71 → 11 orders of margin. Not permissive:
  tightened by 25 orders relative to the v2.0 *measured* value.
- **V2:** unchanged (passed in v2.0).
- **V2b (NEW):** seeding completeness. Requires: seeded count ==
  reference count (104 for `|z| ≤ 300`); every seeded zero within 0.01
  of a known `±2γ_n`; every known `±2γ_n` within 0.01 of a seeded zero;
  worst `|H_0| < 1e-25`. This gate directly bars the v2.0 failure mode.
- **V3:** via `pde_residual_exact` (analytic `dH/dt`); threshold
  **1e-60** (was 1e-20 via finite difference). Measured ~1e-85.
- **V4:** construction and M-refinement unchanged (< 1e-25); PDE
  sub-check via `pde_residual_exact`, threshold **1e-60**.
- **V5:** unchanged; runs after Phase 1 (pending at gate time).
- **V6:** unchanged (passed in v2.0).
- **CAL-(i):** unchanged (passed in v2.0).
- **DH case:** **DROPPED** (carried from v2.0). The comparison proceeds
  as zeta vs synth-δ1, synth-δ01, null-window, per v2.0 §4d.

Phase 1 and Phase 2 remain barred until every runnable v2.1 gate passes
and is logged. The decision rule (v2.0 §9) is unchanged.

## §5. What is unchanged from v2.0

The scientific question, the window (`|z| ≤ 300`, not narrowed), the
event taxonomy (R/C/X), the genuineness protocol G-a/b/c (G-c now uses
the exact PDE residual; threshold 1e-20 retained as a loose bound),
the calibration (§8, including the rigorous-direction synth-δ1 Type-C
requirement), the null-window detector check, the 48-hour study box
(which starts at the first v2.1 `t≠0` evaluation), the red-flag
discipline, and the scope cap (methods/phenomenology only; never a
proof claim; never about RH's truth value).

## §6. Timestamping and the `t≠0` bar

This amendment is deposited on Zenodo **before any v2.1 `t≠0`
computation**. The v2.1 gate run (which performs the first v2.1 `t≠0`
evaluations, inside V6) begins only after this deposit is timestamped.
The v2.0 preregistration (DOI 10.5281/zenodo.23120530) stands unaltered;
this amendment is a separate dated record linked to it.

---

*Signed: Abby Davis, 2026-10-03.*
*Discipline: every day we try. Preregister before computing. Red-flag
every step. Nulls stay in the lab notebook; only sound results see
daylight.*
