#!/usr/bin/env python3
"""Run exact replay or the complete PVRIG candidate-generation workflow."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "models/pvrig_portfolio_rank_v2.json"
CONFIG = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def run(command: list[str], *, env: dict[str, str] | None = None) -> None:
    print("+", " ".join(command))
    subprocess.run(command, cwd=ROOT, env=env, check=True)


def require_inputs() -> tuple[Path, Path, Path]:
    pdb_8x6b = ROOT / "inputs/8X6B.pdb"
    pdb_9e6y = ROOT / "inputs/9E6Y.pdb"
    positive = (
        ROOT
        / "vendor/ab-data-validator/src/ab_data_validator/data/positive.csv"
    )
    missing = [path for path in (pdb_8x6b, pdb_9e6y, positive) if not path.exists()]
    if missing:
        run([sys.executable, "scripts/fetch_inputs.py"])
    return pdb_8x6b, pdb_9e6y, positive


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("replay", "full"), default="replay")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    pdb_8x6b, pdb_9e6y, positive = require_inputs()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir
        else ROOT / "output" / args.mode
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    stage1_json = output_dir / "stage1_candidates.json"
    stage1_csv = output_dir / "stage1_candidates.csv"

    if args.mode == "replay":
        shutil.copyfile(ROOT / "snapshots/stage1_candidates.json", stage1_json)
    else:
        stage1 = CONFIG["stage1"]
        run([
            sys.executable,
            "-m",
            "pvrig_portfolio_rank.design",
            "--positive-csv",
            str(positive),
            "--pdb-8x6b",
            str(pdb_8x6b),
            "--pdb-9e6y",
            str(pdb_9e6y),
            "--model",
            stage1["esm_model_id"],
            "--seed",
            str(CONFIG["random_seed"]),
            "--pool-size",
            str(stage1["pool_size"]),
            "--output-json",
            str(stage1_json),
            "--output-csv",
            str(stage1_csv),
        ])

    final_json = output_dir / "final_candidates.json"
    audit_json = output_dir / "final_audit.json"
    run([
        sys.executable,
        "-m",
        "pvrig_portfolio_rank.finalize",
        "--input-json",
        str(stage1_json),
        "--positive-csv",
        str(positive),
        "--model-config",
        str(CONFIG_PATH),
        "--output-json",
        str(final_json),
        "--audit-json",
        str(audit_json),
    ])

    validator_src = ROOT / "vendor/ab-data-validator/src"
    validation_env = dict(os.environ)
    current_path = validation_env.get("PYTHONPATH", "")
    validation_env["PYTHONPATH"] = (
        str(validator_src) if not current_path else f"{validator_src}{os.pathsep}{current_path}"
    )
    run([
        sys.executable,
        "-m",
        "pvrig_portfolio_rank.validate",
        "--candidates-json",
        str(final_json),
        "--positive-csv",
        str(positive),
        "--threshold",
        str(CONFIG["stage2"]["validation_identity_threshold"]),
        "--output-json",
        str(output_dir / "validator_report.json"),
        "--output-csv",
        str(output_dir / "validator_failures.csv"),
    ], env=validation_env)
    print(f"reproduction complete: {output_dir}")


if __name__ == "__main__":
    main()
