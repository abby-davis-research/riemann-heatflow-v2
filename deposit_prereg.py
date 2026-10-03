#!/usr/bin/env python3
"""Zenodo preregistration deposit for the heat-flow v2.0 study (Oct 3, 2026).

Deposits the v2.0 preregistration PDF as a dated public record BEFORE any
v2 t!=0 computation is run, per the preregistration's own timestamp rule and
Abby's standing independent check (Claude review 2026-10-03: no complaints).

Auth: stored custom.zenodo connector via surrogate query param. Never prints
or logs the credential surrogate.
"""
import hashlib
import json
import os
import sys
import time
import urllib.request
import urllib.error

sys.path.insert(0, "/opt/hatch/skills/skill-creator/bin")
from dynamic_credentials import (  # noqa: E402
    url_with_surrogate_query_param,
    read_json_response,
)

CRED = "custom.zenodo"
HOSTS = ["zenodo.org"]
BASE = "https://zenodo.org/api/deposit/depositions"

V2DIR = os.path.expanduser(
    "~/workspace/goals/a-scaling-theory-of-the-critical-line/hidden_files/riemann-heatflow-v2")
PDF = os.path.join(V2DIR, "preregistration-v2.pdf")
MD = os.path.join(V2DIR, "preregistration.md")

CREATOR = [{"name": "Davis, Abby",
            "affiliation": "Independent researcher, Tucson AZ",
            "orcid": "0009-0002-6758-2263"}]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def auth_url(url):
    return url_with_surrogate_query_param(
        url, CRED, entry_name="access_token", allowed_hosts=HOSTS)


def call(method, url, data=None, ctype="application/json", retries=8):
    last = None
    for a in range(retries):
        try:
            req = urllib.request.Request(auth_url(url), data=data, method=method,
                                         headers={"Content-Type": ctype})
            with urllib.request.urlopen(req, timeout=90) as r:
                return r.status, read_json_response(r)
        except urllib.error.HTTPError as e:
            eb = e.read().decode()[:400]
            print(f"  HTTP {e.code}: {eb}", flush=True)
            last = e
            if e.code in (400, 401, 403, 404):
                raise
        except Exception as e:
            print(f"  attempt {a+1} dropped ({type(e).__name__}); retrying...",
                  flush=True)
            last = e
        time.sleep(5)
    raise RuntimeError(f"failed after {retries}: {last}")


def main():
    assert os.path.isfile(PDF), f"missing {PDF}"
    assert os.path.isfile(MD), f"missing {MD}"
    pdf_hash = sha256(PDF)
    md_hash = sha256(MD)
    print(f"PDF sha256: {pdf_hash}", flush=True)
    print(f"MD  sha256: {md_hash}", flush=True)

    print("1. creating deposition...", flush=True)
    st, d = call("POST", BASE, data=b"{}")
    dep_id = d["id"]
    bucket = d["links"]["bucket"]
    print(f"  id={dep_id}", flush=True)

    print("2. uploading files...", flush=True)
    for path, name in [(PDF, "preregistration-v2.pdf"), (MD, "preregistration.md")]:
        with open(path, "rb") as f:
            raw = f.read()
        call("PUT", f"{bucket}/{name}", data=raw,
             ctype="application/octet-stream")
        print(f"  uploaded {name} ({len(raw)} bytes)", flush=True)

    metadata = {
        "title": ("Preregistration v2.0: Mapping the high-|z| threshold event "
                  "structure of de Bruijn-Newman heat flows"),
        "upload_type": "publication",
        "publication_type": "workingpaper",
        "publication_date": "2026-10-03",
        "description": (
            "PREREGISTRATION (not a result). Written 2026-10-03 BEFORE any "
            "v2 t!=0 computation was run, and deposited the same day as a "
            "public timestamp of the study design. This is a NEW "
            "preregistration, not a continuation or rescue of the author's "
            "2026-10-02 heat-flow milestone, which returned METHOD INVALID "
            "under its own decision rule (its directional premise - lowest "
            "real zero pair collides going down in t - proved empirically "
            "false for the verified implementation; PDE residual <= 3.5e-22).\n\n"
            "The v2 study maps the real<->complex transition event structure "
            "of de Bruijn-Newman-type heat flows across a WIDE zero window "
            "(|z| <= 300), where post-hoc scans suggest the threshold dynamics "
            "actually live (|z| ~ 190-230, consistent with the rigorous bound "
            "Lambda >= 0). Design: no directional assumption (characterize the "
            "zeta event map first, compare second); all six verification gates "
            "must pass AND be logged before Phase 2; per-event genuineness "
            "protocol (t-continuation, J-refinement, PDE residual) with a "
            ">50% artifact rate invalidating the apparatus; calibration on the "
            "rigorous real-ification direction via synthetic injected-pair "
            "landings; a null-window detector check on zeta over (0.2, 0.5] "
            "as the false-positive control.\n\n"
            "SCOPE CAP (hard ceiling): methods/phenomenology only. NOT a proof "
            "of RH, NOT a statement about RH's truth value, NOT a measurement "
            "of Lambda. A null result is a complete lab-notebook result.\n\n"
            f"File integrity: preregistration-v2.pdf sha256 = {pdf_hash}; "
            f"preregistration.md sha256 = {md_hash}."),
        "creators": CREATOR,
        "license": "cc-by-4.0",
        "keywords": ["preregistration", "de Bruijn-Newman", "heat flow",
                     "Riemann zeta", "research methods", "reproducibility"],
        "prereserve_doi": True,
    }

    print("3. setting metadata + reserving DOI...", flush=True)
    st, d = call("PUT", f"{BASE}/{dep_id}",
                 data=json.dumps({"metadata": metadata}).encode())
    doi = d["metadata"]["prereserve_doi"]["doi"]
    print(f"  reserved DOI: {doi}", flush=True)

    print("4. publishing...", flush=True)
    st, d = call("POST", f"{BASE}/{dep_id}/actions/publish")
    final_doi = d.get("doi", doi)
    print(f"  PUBLISHED: {final_doi}", flush=True)

    ts_path = os.path.join(V2DIR, "TIMESTAMP.md")
    with open(ts_path, "w") as f:
        f.write(
            "# v2.0 preregistration timestamp\n\n"
            f"- Deposited: 2026-10-03 (America/Phoenix)\n"
            f"- Zenodo DOI: {final_doi}\n"
            f"- Record: https://zenodo.org/records/{final_doi.split('/')[-1]}\n"
            f"- preregistration-v2.pdf sha256: {pdf_hash}\n"
            f"- preregistration.md sha256: {md_hash}\n"
            "- No v2 t!=0 computation had been run at deposit time.\n")
    print(f"wrote {ts_path}", flush=True)
    print(f"\nRESULT: {final_doi}")


if __name__ == "__main__":
    main()
