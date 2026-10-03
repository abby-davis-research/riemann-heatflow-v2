# Preregistration v2.0: Mapping the high-|z| threshold event structure of de Bruijn–Newman heat flows

**Date:** 2026-10-03 · **Author:** Abby Davis (independent researcher, Tucson AZ; ORCID 0009-0002-6758-2263), with Melody/Muse as instrument
**Status:** WRITTEN BEFORE ANY v2 t≠0 COMPUTATION. No flow results for this preregistration have been observed. v1.0's data and code exist and are used ONLY for calibration design, never as v2 evidence.
**Relationship to v1.0:** This is a NEW preregistration, not a continuation or rescue of the 2026-10-02 heat-flow milestone (verdict: METHOD INVALID). It is written to be executable by an independent analyst.

## 1. What v1.0 established (and what it forbids)

v1.0 (dir: `riemann-heatflow/`, verdict METHOD INVALID per its own §6) established:

- (a) The preregistered premise "t decreasing ⇒ lowest real zero pair approaches and collides" is FALSE for the verified implementation H_t(z) = Σ c_j e^{tu_j²}cos(u_j z) (PDE residual ≤ 3.5e−22 — the code is correct; the premise was wrong). Low zeros SPREAD as t decreases.
- (b) A post-hoc wide scan found ~20+ nonreal zeros at HIGH |z|≈190–230 for t<0, consistent with the rigorous bound Λ ≥ 0 (Rodgers–Tao 2020). The threshold dynamics live at high |z|, not the lowest pair.
- (c) Verification gates V1/V2/V4/V5 were never run or logged in v1.0; synthetic controls saved no trajectories; the Davenport–Heilbronn off-axis pair was never located.
- (d) One high-z complex zero appeared t-stationary ("frozen") — possibly a quadrature artifact; unresolved.

v1.0 FORBIDS: re-pointing at high-z zeros inside the v1.0 milestone (moving the goalposts); any discriminator claim from v1.0 data; assuming a collision direction anywhere in this design.

## 2. Research question and scope

**Question.** Under the de Bruijn–Newman heat flow, what is the real↔complex transition event structure in a WIDE zero window (|z| ≤ 300) — and does any of its preregistered features distinguish the ζ flow from RH-false flows when produced by a validated, artifact-quantified apparatus?

**Scope (hard ceiling).** Methods/phenomenology paper: "the high-zero threshold event structure of de Bruijn–Newman-type heat flows." NOT a proof of RH, NOT a statement about RH's truth value, NOT a measurement of Λ. A null (no discriminator found) is a complete lab-notebook result and is not published.

## 3. Conventions

Wide frame (Rodgers–Tao / Polymath): s = (1+iz)/2. ζ: H_0(z) = ξ(1/2+iz/2)/8, zeros at z = ±2γ_n. Heat flow: H_t(z) = ∫_0^∞ Φ(u)e^{tu²}cos(zu)du, satisfying ∂_t H_t = −∂_{zz}H_t. de Bruijn–Newman constant Λ: H_t is real-rooted ⟺ t ≥ Λ; proved 0 ≤ Λ ≤ 0.2 (lower: Rodgers–Tao 2020; upper: Platt–Trudgian 2021).

**Event taxonomy (fixed in advance).**

- Type-R: two adjacent real zeros collide (gap → 0); below the collision a complex-conjugate pair emerges. Detected wherever it occurs — no direction assumed.
- Type-C: a complex-conjugate pair's |Im| → 0 at a common real point (lands on ℝ).
- Type-X: collision of zeros off the real axis that are not conjugates (codimension 2; expected never; flagged, never fitted).
- ARTIFACT: a detector trigger failing the genuineness protocol (§6). Counted and reported, never mapped.

## 4. Cases

**(a) ζ (reference).** H_t via fixed Gauss–Legendre quadrature, J=256 nodes on [0,3], mpmath dps=50, Φ(u) = Σ_{n=1}^{8} (2π²n⁴e^{9u} − 3πn²e^{5u}) e^{−πn²e^{4u}}. (This corrects v1.0's umax=1.0 deviation: the interval is [0,3] and is gated.)

**(b) Synth-δ1 ("gross" RH-false, ground truth by construction).** H_0^syn(z) = H_0(z)·R(z) with R(z) = [(1−z²/w²)(1−z²/w̄²)]/(1−z²/z_1²), w = z_1 + i·1.0. Surgically replaces the real pair ±z_1 with the quadruplet ±w, ±w̄; R is real on ℝ, R(0)=1, R→1 at ∞, and H_0·R is entire (denominator pole cancels against H_0's simple zero). N_off(0) = 4 by construction. H_t^syn = e^{−t∂_{zz}}[H_0·R] via Taylor-in-t with Leibniz rule, truncation M=48 (gated against M=96).

**(c) Synth-δ0.1 ("near-miss" RH-false, sensitivity probe).** Same as (b) with δ = 0.1 (Rodgers–Tao near-miss scale).

**(d) Davenport–Heilbronn (authentic RH-false Dirichlet series). BEST-EFFORT, GATED.** Wide-frame Φ_DH(v) = 8e^{3v} Σ_{n=1}^{20} n·a_n e^{−πn²e^{4v}/5}, a = (1, κ, −κ, −1, 0) period 5, κ = (√(10−2√5)−2)/(√5−1) (expected 0.284079043840…). Gates (ALL must pass): G1 functional equation |Λ_DH(s) − Λ_DH(1−s)| < 10⁻²⁰ at 5 random s; G2 reality |Im Λ_DH(1/2+it)| < 10⁻²⁰ at 5 random t; G3 quadrature cross-check |∫Φ_DH cos − Λ_DH((1+iz)/2)| < 10⁻²⁵ (up to one fitted overall scale) at 5 random z. Time-box: 45 minutes wall-clock. On timeout or gate failure: DH is dropped, the drop is logged as a red flag, and the comparison proceeds as ζ vs (b), (c), (e).

**(e) Null-window detector check.** The ζ up-leg interval (0.2, 0.5] is RIGOROUSLY event-free (real-rootedness is monotone in t for t ≥ Λ). Running the full event detector there measures its false-positive rate on real data under a theorem-grade null — this plays the role the Poisson control played in the anti-resonance study. Any detected event → detector bug → INVALID (apparatus).

## 5. Tracking protocol (identical across cases unless noted)

- Window: every zero with |z| ≤ 300 at t=0 (≈60+ for ζ: ±2γ_1…±2γ_30; synth cases: ±w, ±w̄ plus ζ-like rest; DH: as located, or dropped per §4).
- t=0 seeding: grid-seeded complex Newton (tol |Δz| < 10⁻³⁰, max 25 iters, require |H_t(z)| < 10⁻²⁵ at convergence), deduplicated by pairwise distance.
- Legs: down-leg t ∈ [0, −1.0]; up-leg t ∈ [0, 0.5] (synthetic landing checks; one preregistered extension to 2.0 per §8(iv); DH off-axis pair up to 1.0 if gated in).
- Step control: Δt = 0.01 initial; halve if Newton iters > 12 (floor 0.001); double if iters < 4 (cap 0.02); reject any step with |Δz| > 0.5 × neighboring gap, halve Δt and retry.
- Detection: monitor all adjacent real gaps g_k(t) and all complex pairs' |Im z_k(t)|. Trigger R when g_k < 0.03 (log (t, g_k²) each step; at g_k < 10⁻⁷ or Newton failure, fit g² = C²(t−t*) on points with g ∈ [10⁻⁷, 3×10⁻²], require ≥5 points and R² ≥ 0.999). Trigger C when |Im| < 0.03 (fit im² linearly under the same quality rules). Trigger X is logged with coordinates and never fitted.
- Trajectories: FULL JSON per case — every t-step, every tracked zero (z, |H|, Newton iters, Δt). v1.0 failed here; in v2 this is mandatory. A case without saved trajectories is VOID.

## 6. Per-event genuineness protocol (answers v1.0's frozen-zero problem)

Every triggered event must pass ALL of the following, or it is an ARTIFACT:

- G-a (t-continuation): the participating zeros track smoothly (per-step |Δz| within the §5 step-rejection bounds) for ≥10 steps on each side of the event. Non-smoothness → SUSPECT: logged, excluded from the map.
- G-b (J-refinement): recompute the event at J=512; require |Δt*| < 3σ(fit) and |Δz*| < 0.05. Failure → ARTIFACT.
- G-c (PDE): |∂_t H + ∂_{zz} H| < 10⁻²⁰ at the fitted (z*, t*).

Artifact failures are logged with reasons and tallied; the artifact rate A/(A+|E|) is a reported observable. If it exceeds 0.5 for ζ, the apparatus is declared UNFIT for high-z work → INVALID (apparatus). That is a verdict about the instrument, not the theory.

## 7. Verification gates (ALL must pass AND be logged before Phase 2 begins)

- V1 (quadrature): J=256 → J=512 changes H_t(z) by < 10⁻³⁵ at 10 test (z,t), INCLUDING 4 points with |z| ∈ [180, 280]. (v1.0 never tested high z.)
- V2 (t=0 zeros): |H_0(2γ_n)| < 10⁻³⁰ for n = 1…5; 14-digit agreement with 2γ_n.
- V3 (PDE): |∂_t H + ∂_{zz} H| < 10⁻²⁰ at 10 test (z,t) per case. (v1.0 passed this unlogged; re-run and LOG it.)
- V4 (synthetic): M=48 vs M=96 agree to < 10⁻²⁵ at all tracked (z,t); H_0^syn(w) = 0 to < 10⁻²⁵ (construction check); PDE residual < 10⁻²⁰.
- V5 (seed independence): the full ζ down-leg is repeated with Newton seeds perturbed by 10⁻⁶; the VALIDATED event map must agree — same events within 3σ, none gained, none lost.
- V6 (detector null): the full event detector is run on ζ over (0.2, 0.5]; zero events required.

Any gate failure → STOP, diagnose, red-flag. Phase 2 does not begin until every gate passes and its log exists on disk.

## 8. Calibration (ζ + synth must satisfy ALL, else METHOD INVALID)

- (i) N_off(0) = 0 for ζ in |z| ≤ 300.
- (ii) V1–V6 all pass and are logged.
- (iii) Zero validated events in ζ over (0.2, 0.5] (same computation as V6; restated here as the flow-level check).
- (iv) APPARATUS VALIDATION — the rigorous direction: synth-δ1's injected pair lands via Type-C at t*_inj with 0 < t*_inj ≤ 2.0. The up-leg is tracked to 0.5 and extended ONCE to 2.0 if and only if no landing occurred by 0.5 (preregistered extension, not peeking). The landing requires an im²-fit with R² ≥ 0.999, and N_off must drop by exactly the injected count at the landing. If no landing by t = 2.0 → INVALID.

Item (iv) replaces v1.0's dead calibration (ii) with a claim the theorems actually support: the heat flow real-ifies as t increases (de Bruijn). v1.0 calibrated on the direction the flow does NOT go; v2 calibrates on the direction it provably does.

## 9. Phases, observables, decision rule

**Phase 1 (characterization, ζ only).** Produce the validated event map E_ζ and the artifact tally A_ζ. This phase is descriptive; it makes no discriminator claim. If E_ζ is EMPTY (no genuine events survive G-a/b/c), report that plainly — Phase 2 still runs, and the comparison becomes "empty vs non-empty maps."

**Phase 2 (comparison).** The identical protocol on cases (b), (c), (e); case (d) if gated in.

**Observables.** E = validated event map {(type, t* ± σ, z*)}; A = artifact count and rate; N_off(t) curves per case; t*_inj for the synthetic cases.

**Decision rule (fixed now).**

- SUCCESS: calibration (§8) passes AND at least one preregistered feature separates ζ from at least one RH-false case at >3σ. Preregistered features: |E| (event count); the distribution of |z*| (two-sample Kolmogorov–Smirnov, α = 0.01 — if either map has fewer than 5 events this feature scores as non-separating for lack of data, never as a success); t*_max (latest event time); N_off(−1). Stability of the separating feature across the V5 seed check is required.
- FAILURE: calibration passes but no feature separates (all within 3σ) — INCLUDING the sub-case of an empty E_ζ.
- INVALID: any calibration item fails, any gate fails, or the ζ artifact rate exceeds 0.5.

**Interpretation cap.** SUCCESS means: "the high-z heat-flow event structure is a validated dynamical discriminator between RH-true and RH-false flows on the tested cases." It does NOT imply anything about the truth value of RH. It is NOT a measurement of Λ (all t* values are windowed).

## 10. Compute boxes

- DH gating: 45 minutes wall-clock (as v1.0).
- ζ wide-window down-leg: 6 hours wall-clock; on timeout, report the partial map with the box noted — a preregistered scope cut, not a failure.
- Synthetic cases: 3 hours each.
- Total study: 48 hours wall-clock from the first v2 t≠0 evaluation.

## 11. Deliverables

- `code/`: ht_flow.py (quadrature core + H_t and derivatives), track.py (wide-window Newton tracker + event detector + genuineness protocol), synthetic.py, dh.py (gated), analyze.py, make_figures.py.
- `data/`: per-case trajectory JSONs (mandatory), event table CSV including genuineness columns, gate-check log (mandatory), artifact log.
- `figures/`: (1) ζ event map (t vs z, type-encoded); (2) fold/landing fits with residuals; (3) N_off(t) for all cases; (4) artifact gallery (triggers that failed G-a/b/c, with reasons).
- `results.md`: quantitative outcomes, the full red-flag log, and a soundness self-assessment.
- This preregistration, with its sha256 recorded in results.md.

## 12. Pre-declared red flags and known risks

- R1 (v1.0 R-A): never assume a collision direction again. The protocol detects R/C/X wherever they occur; the map is the datum.
- R2 (v1.0 R-B): frozen / quadrature-stationary zeros at high z. G-a/b/c exist for exactly this; artifact rate > 0.5 → INVALID (apparatus).
- R3 (v1.0 R-C): the DH off-axis pair may again evade location → gated drop, logged, comparison proceeds without it.
- R4 (v1.0 R-D): synthetic trajectories are MANDATORY; a case without saved trajectories is VOID, not "log-only."
- R5 (v1.0 R-E): umax/J/dps are as specified in §4; any deviation is a protocol deviation — logged, and the affected gates re-run.
- R6: Newton degradation near folds (H' → 0 at the fold): fit-early strategy as v1.0; never Newton into g < 10⁻⁷.
- R7: label jumping between adjacent zeros at small gaps: step rejection plus a post-hoc trajectory-continuity check; jumps trigger re-tracking or a flag.
- R8: compute cost at |z| ~ 300 (dps=50, J=256/512) may force J/dps tradeoffs; any change is a deviation — logged and re-gated per R5.
- R9 (v1.0 R5, "the boring reason"): if ONLY synth-δ1 separates, the claim is weakened to "detects gross violations" and reported as such. The δ=0.1 case is the sensitivity probe that earns a stronger claim.
- R10: the t*_inj extension (0.5 → 2.0) is preregistered SOLELY for calibration item (iv). It is not a general license to extend windows.

---

*End of preregistration v2.0. No v2 t≠0 computation may precede the sha256 timestamp recorded in results.md.*
