#!/usr/bin/env python3
"""Export a flat candidate CSV and maintain the result SHA-256 manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def csv_rows(payload: dict) -> tuple[list[str], list[dict]]:
    fields = [
        "rank",
        "name",
        "framework",
        "sequence",
        "cdr1",
        "cdr2",
        "cdr3",
        "parent",
        "max_positive_identity",
        "calibrated_prior_score",
        "weight_sensitivity_score_p10",
        "weight_sensitivity_score_p90",
        "weight_sensitivity_rank_p10",
        "weight_sensitivity_rank_p90",
        "theoretical_pI",
        "gravy",
        "evidence_level",
    ]
    rows = []
    for item in payload["candidates"]:
        rows.append({
            "rank": item["rank"],
            "name": item["name"],
            "framework": item["framework"],
            "sequence": item["sequence"],
            "cdr1": item["cdr1"],
            "cdr2": item["cdr2"],
            "cdr3": item["cdr3"],
            "parent": item["parent"]["name"],
            "max_positive_identity": item["max_positive_identity"],
            "calibrated_prior_score": item["calibrated_prior_score"],
            "weight_sensitivity_score_p10": item["weight_sensitivity_score_p10"],
            "weight_sensitivity_score_p90": item["weight_sensitivity_score_p90"],
            "weight_sensitivity_rank_p10": item["weight_sensitivity_rank_p10"],
            "weight_sensitivity_rank_p90": item["weight_sensitivity_rank_p90"],
            "theoretical_pI": item["full_sequence_metrics"]["theoretical_pI"],
            "gravy": item["full_sequence_metrics"]["gravy"],
            "evidence_level": item["evidence_level"],
        })
    return fields, rows


def expected_csv(payload: dict) -> str:
    from io import StringIO

    fields, rows = csv_rows(payload)
    handle = StringIO(newline="")
    writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return handle.getvalue()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    source = RESULTS / "final_candidates.json"
    target = RESULTS / "final_candidates.csv"
    payload = json.loads(source.read_text(encoding="utf-8"))
    rendered = expected_csv(payload)
    if args.check:
        if not target.exists() or target.read_text(encoding="utf-8") != rendered:
            raise SystemExit("results/final_candidates.csv is stale")
    else:
        target.write_text(rendered, encoding="utf-8")

    manifest_paths = [
        RESULTS / "final_candidates.json",
        RESULTS / "final_candidates.csv",
        RESULTS / "final_audit.json",
        RESULTS / "validator_report.json",
        ROOT / "snapshots/stage1_candidates.json",
        ROOT / "models/pvrig_portfolio_rank_v2.json",
    ]
    lines = [f"{sha256(path)}  {path.relative_to(ROOT)}" for path in manifest_paths]
    manifest = "\n".join(lines) + "\n"
    manifest_path = RESULTS / "manifest.sha256"
    if args.check:
        if not manifest_path.exists() or manifest_path.read_text(encoding="utf-8") != manifest:
            raise SystemExit("results/manifest.sha256 is stale")
    else:
        manifest_path.write_text(manifest, encoding="utf-8")
    print("result export verified" if args.check else "result export updated")


if __name__ == "__main__":
    main()
