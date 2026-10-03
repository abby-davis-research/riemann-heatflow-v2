# CODE_NOTES.md — v2 heat-flow study code (preregistration v2.0)

**Status (2026-10-03): code complete. NO v2 t≠0 computation has been run.**
`ht_flow.py` is the byte-identical verified v1 core (PDE residual ≤ 3.5e−22);
its quadrature was NOT modified. All tolerances/thresholds below come
literally from prereg §§4–8; deviations and judgment calls are flagged.

## Modules

- **ht_flow.py** — frozen. `HeatFlow(phi, J, umax, dps)`: fixed-node
  Gauss–Legendre quadrature of H_t(z) = ∫Φe^{tu²}cos(zu)du, plus `H/Hp/Hpp`,
  `H0_derivs`, `pde_residual`, and the DH Hurwitz helpers (`dh_f`,
  `dh_Lambda`, `phi_dh_wide`, `dh_kappa`). v2 uses `J=256, umax=3.0, dps=50`
  (prereg §4) via constructor args — no code change needed.
- **track.py** — wide-window Newton tracker (§5) + event detector +
  genuineness protocol (§6) + mandatory trajectory serialization (§11).
  - `newton(flow, z0, t)`: tol |dz|<1e−30, max 25 iters, require |H|<1e−25.
  - `seed_zeros(flow, zmax=300)`: grid-seeded complex Newton at t=0,
    dedup by pairwise distance; `label_zeros` assigns deterministic
    `r±{i}` / `c±{i}` labels.
  - `track_leg(flow, seeds, t_start, t_end)`: adaptive Δt (0.01; halve if
    iters>12 floor 0.001; double if iters<4 cap 0.02; reject |dz|>0.5×gap,
    halve and retry); per-step event detection (R/C/X); returns full
    trajectories `(t, z, |H|, iters, dt)` per label — MANDATORY, a case
    without them is VOID.
  - `fit_fold(pts)`: g²/im²-vs-t linear fit, ≥5 pts, R²≥0.999, window
    [1e−7, 3e−2].
  - `validate_event(flow256, flow512, event, traj)`: G-a (≥10 smooth steps
    each side), G-b (J=512 refit, |Δt*|<3σ, |Δz*|<0.05), G-c (PDE<1e−20).
    Verdicts: genuine / suspect / artifact / logged (X).
  - `save_case_json(...)`: prereg §11 case file schema.
- **synthetic.py** — `SyntheticFlow(base, z1, delta, M=48)`: R(z) surgery
  (§4), Taylor-in-t via Leibniz rule, `H/Hp/Hpp`, `pde_residual` (for G-c/V4),
  `construction_check()` (|H_0^syn(w)|, V4). F-coefficients cached per z
  (t-independent) for performance.
- **dh.py** — DH gates G1/G2/G3 with JSON-LOGGED results, seeded RNG (7),
  45-minute wall-clock box (`BoxExpired` → DROPPED, logged). `make_dh_flow`,
  `find_offaxis_zeros` (genuinely-off-axis only: |Re(s)−1/2|>1e−3).
- **analyze.py** — `load_case_json`, `validated_events`, `ks_2samp`
  (two-sample KS, α=0.01), `n_off_curve` (|Im|>1e−6), `artifact_tally`,
  gates `gate_V1`…`gate_V6` + `log_gates`, `check_calibration` (§8),
  `check_landing` (§8(iv): 0<t*_inj≤2.0, R²≥0.999, N_off drop == 4),
  `separation_stats` (§9 features), `v5_stability`, `compare_event_maps`
  (V5), `confirm_V4_tracked`, `write_events_csv`, `write_figure_data`
  (fig1–fig4 data JSONs).

## What remains (not yet written/run)

1. **Runner scripts**: `run_gates.py` (V1–V6 → gate log), `run_phase1.py`
   (ζ characterization → event map + artifact tally), `run_phase2.py`
   (synth/DH/null-window comparison → decision rule). The task scoped code
   only; runners are the next step.
2. **Gate runs** (all t≠0 work is barred until the runner stage).
3. **make_figures.py** (prereg §11): figure *data* writers exist in
   analyze.py; the plotting script itself is not written.
4. `data/` and `figures/` dirs exist and are empty.

## Judgment calls & flagged items (prereg underspecification)

1. **DH Φ sum**: prereg §4 writes n=1..20; frozen `ht_flow.phi_dh_wide`
   sums n=1..30. More terms can only help; G3 is the arbiter. Logged in
   `dh.py`, not silently changed.
2. **Dedup distance** for t=0 seeding: 1e−3 (prereg: "pairwise distance",
   no value). Converged Newton zeros agree to ~1e−25, so 1e−3 is safe.
3. **Adaptive-dt "iters"**: uses the MAX Newton iterations over labels at
   the step (prereg doesn't say max vs mean; max is conservative).
4. **V4 "at all tracked (z,t)"** is circular pre-tracking (gates run before
   Phase 2). Implemented as: gate on 10 test (z,t) points pre-Phase-2, PLUS
   `confirm_V4_tracked` re-run on sampled tracked points post-tracking.
5. **>3σ for count features** (|E|, N_off(−1)): prereg doesn't define σ.
   Implemented as Poisson z = |a−b|/√(a+b) > 3. **CONFIRMED by Abby
   2026-10-03 ("it's good") — proceed to Phase 2 on this reading.**
6. **Newton-exhaustion mid-leg** (dt at floor, still failing): the failed
   labels are dropped with a logged TRACK_FAIL; the leg continues. (Prereg
   specifies halve-and-retry but not the terminal case.)
7. **R→C double count**: a nascent complex pair from a just-detected R event
   can trigger the C detector while |Im|<0.03. It will fail G-a (<10 steps
   on the complex side) → SUSPECT/artifact — noisy but prereg-consistent.
   No suppression rule was added (prereg doesn't authorize one).
8. **v1 synthetic Leibniz bug FIXED**: v1 looped j∈[0,min(k,4)], dropping
   h^{(j)}, j∈(4,k] for k>4 (wrong Taylor coefficients for m≥3). v2 uses
   the full j∈[max(0,k−4),k] range with h to order 2M.
9. **Taylor-in-t radius**: convergence at larger t is an empirical question;
   gate V4 (M=48 vs M=96 <1e−25) is the arbiter per prereg. If V4 fails,
   that's a gate failure (STOP), not a tuning exercise.
10. **Injected labels**: `check_landing` takes `injected_labels`
    (default `('w1','w2','w3','w4')`); the Phase-2 runner must label the
    quadruplet accordingly (or pass its own tuple).
11. **X events**: logged with coordinates, never fitted, verdict 'logged';
    they count as detector triggers for V6 (any trigger in (0.2,0.5] fails).
12. **V2 "14-digit agreement"**: implemented as relative |z−2γ_n|/|2γ_n| <
    1e−14 after Newton polish.
13. Recommended (not preregistered) diagnostic for the runner: compare the
    t=0 seed count against Riemann–von Mangoldt N(150)≈53 → ~106 zeros in
    |z|≤300, as a seeding-completeness sanity check (UNGATED).
14. Inherited from frozen ht_flow.py: the constructor drops quadrature nodes
    with |w·Φ(u)| ≤ 1e−40, so J=256 on [0,3] keeps ~92 effective nodes (Φ
    mass beyond u≈1 is ~1e−71; v1's audit). The V1 J=256→J=512 comparison
    exercises the same code path v1 used — the "J" in gates refers to the
    requested node count, as in v1.
15. Post-build addition (2026-10-03): `track_leg` accepts `box_seconds`
    (prereg §10 wall-clock box); on expiry the leg returns the partial map
    with `box_expired: True`. No tracking math was touched. `save_case_json`
    records `box_expired` per leg.
