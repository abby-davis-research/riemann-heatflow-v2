# Preregistration v2.0 — high-|z| threshold event structure of de Bruijn–Newman heat flows

**Status:** preregistered 2026-10-03 · **DOI:** [10.5281/zenodo.23120530](https://zenodo.org/records/23120530)
(deposited BEFORE any v2 t≠0 computation was run)

A NEW preregistration, not a continuation of the 2026-10-02 heat-flow milestone
(verdict: METHOD INVALID — its directional premise proved false for the verified
implementation). This study maps the real↔complex transition event structure of
de Bruijn–Newman-type heat flows across a wide zero window (|z| ≤ 300), where
the threshold dynamics actually live.

**Scope cap (hard ceiling):** methods/phenomenology only. NOT a proof of RH,
NOT a statement about RH's truth value, NOT a measurement of Λ.

## Layout

- `preregistration.md` — the preregistration (sha256 recorded in `TIMESTAMP.md`)
- `preregistration-v2.pdf` — typeset preregistration (7 pp)
- `TIMESTAMP.md` — Zenodo DOI, hashes, deposit record
- `code/` — `ht_flow.py` (frozen verified quadrature core), `track.py`
  (wide-window Newton tracker + event detector + genuineness protocol),
  `synthetic.py`, `dh.py` (gated Davenport–Heilbronn), `analyze.py`
  (gates, calibration, discriminator statistics), `run_gates.py`,
  `run_phase1.py`, `run_phase2.py`, `make_figures.py`, `CODE_NOTES.md`
- `data/` — gate logs, seeds, event CSV, red-flag log, decision (trajectory
  JSONs excluded from git — large files; archived with the results deposit)
- `figures/` — figure data + plots (§11)

## Reproduce

```bash
cd code
python3 run_gates.py    # V1–V6 + calibration (i)–(iii) + DH gating; STOPS on failure
python3 run_phase1.py   # ζ characterization + V5 seed-independence
python3 run_phase2.py   # comparison + calibration (iv) + §9 decision
python3 make_figures.py
```

Phase 2 is barred unless every gate passes and its log exists on disk.
A null result is a complete lab-notebook result and is not published.

## Standing rules

- Preregister before computing. Red-flag every step.
- Null/invalid results stay in the lab notebook; publish only genuinely sound work.
