#!/usr/bin/env python3
"""Fetch and verify the external inputs required for reproduction."""

from __future__ import annotations

import hashlib
import json
import subprocess
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONFIG = json.loads(
    (ROOT / "models/pvrig_portfolio_rank_v2.json").read_text(encoding="utf-8")
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch(url: str, output: Path, expected: str) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and sha256(output) == expected:
        print(f"verified {output.relative_to(ROOT)}")
        return
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "pvrig-portfolio-rank/2.0"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        output.write_bytes(response.read())
    actual = sha256(output)
    if actual != expected:
        raise RuntimeError(
            f"checksum mismatch for {output}: expected {expected}, got {actual}"
        )
    print(f"downloaded and verified {output.relative_to(ROOT)}")


def run(*command: str, cwd: Path | None = None) -> str:
    result = subprocess.run(
        command,
        cwd=cwd,
        check=True,
        text=True,
        capture_output=True,
    )
    return result.stdout.strip()


def fetch_validator() -> None:
    inputs = CONFIG["inputs"]
    destination = ROOT / "vendor/ab-data-validator"
    expected = inputs["organizer_validator_commit"]
    if not destination.exists():
        destination.parent.mkdir(parents=True, exist_ok=True)
        run(
            "git",
            "clone",
            "--filter=blob:none",
            inputs["organizer_validator_repository"],
            str(destination),
        )
    actual = run("git", "rev-parse", "HEAD", cwd=destination)
    if actual != expected:
        run("git", "fetch", "origin", expected, cwd=destination)
        run("git", "checkout", "--detach", expected, cwd=destination)
        actual = run("git", "rev-parse", "HEAD", cwd=destination)
    if actual != expected:
        raise RuntimeError(f"validator commit mismatch: expected {expected}, got {actual}")
    print(f"verified validator commit {actual}")


def main() -> None:
    inputs = CONFIG["inputs"]
    fetch(
        "https://files.rcsb.org/download/8X6B.pdb",
        ROOT / "inputs/8X6B.pdb",
        inputs["pdb_8x6b_sha256"],
    )
    fetch(
        "https://files.rcsb.org/download/9E6Y.pdb",
        ROOT / "inputs/9E6Y.pdb",
        inputs["pdb_9e6y_sha256"],
    )
    fetch_validator()


if __name__ == "__main__":
    main()
