#!/usr/bin/env python3
"""Run the organizer's validator logic with task-local ANARCII2/Bio.Align adapters."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from anarcii import Anarcii
from Bio.Align import PairwiseAligner

from ab_data_validator.models import AntibodyRow
from ab_data_validator.numbering import NumberedResidue
from ab_data_validator.positive_library import load_positive_library
from ab_data_validator.validation import Validator


class CachedAnarciiNumberer:
    def __init__(self, sequences: list[str]) -> None:
        engine = Anarcii()
        unique = list(dict.fromkeys(sequences))
        results = engine.number(unique)
        self.cache: dict[str, list[NumberedResidue]] = {}
        if len(results) != len(unique):
            raise RuntimeError(f"ANARCII returned {len(results)} results for {len(unique)} sequences")
        for sequence, result in zip(unique, results.values(), strict=True):
            if result.get("error") or not result.get("numbering"):
                raise RuntimeError(f"ANARCII failed: {result.get('error')}")
            self.cache[sequence] = [
                NumberedResidue(
                    position=int(position),
                    insertion=str(insertion).strip(),
                    residue=residue,
                )
                for (position, insertion), residue in result["numbering"]
                if residue not in {"-", "."}
            ]

    def number(self, sequence_id: str, sequence: str, chain: str):
        del sequence_id, chain
        return self.cache[sequence]


class GlobalAligner:
    def __init__(self) -> None:
        self.engine = PairwiseAligner(
            mode="global",
            match_score=2.0,
            mismatch_score=-1.0,
            open_gap_score=-1.0,
            extend_gap_score=-1.0,
        )

    def align(self, cdr_name: str, candidate_cdr: str, positive_cdr: str):
        del cdr_name
        alignment = self.engine.align(candidate_cdr, positive_cdr)[0]
        return str(alignment[0]), str(alignment[1])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates-json", type=Path, required=True)
    parser.add_argument("--positive-csv", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=0.75)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.candidates_json.read_text(encoding="utf-8"))
    candidates = [
        AntibodyRow(name=item["name"], vh=item["sequence"], vl=None)
        for item in payload["candidates"]
    ]
    positives = load_positive_library(args.positive_csv)
    sequences = [row.vh for row in positives + candidates]
    sequences.extend(row.vl for row in positives if row.vl)
    numberer = CachedAnarciiNumberer(sequences)

    progress = []
    validator = Validator(
        numberer=numberer,
        aligner=GlobalAligner(),
        identity_threshold=args.threshold,
        max_workers=1,
        progress_logger=progress.append,
    )
    failures = validator.validate(candidates, positives)

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    result = {
        "validator": "organizer ab-data-validator validation.Validator",
        "numbering_adapter": "ANARCII 2.0.8, IMGT",
        "alignment_adapter": "Bio.Align PairwiseAligner global; match=2, mismatch=-1, gap=-1",
        "identity_threshold": args.threshold,
        "candidate_count": len(candidates),
        "positive_reference_count": len(positives),
        "failure_count": len(failures),
        "passed_count": len(candidates) - len({failure.name for failure in failures}),
        "progress": progress,
        "failures": [
            {
                "name": failure.name,
                "reason_type": failure.reason_type,
                "chain": failure.chain,
                "cdr": failure.cdr,
                "positive_name": failure.positive_name,
                "identity": failure.identity,
                "threshold": failure.threshold,
                "details": failure.details,
            }
            for failure in failures
        ],
    }
    args.output_json.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    with args.output_csv.open("w", newline="", encoding="utf-8-sig") as handle:
        fieldnames = [
            "name", "reason_type", "chain", "cdr", "positive_name",
            "identity", "threshold", "details",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for item in result["failures"]:
            writer.writerow(item)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    raise SystemExit(0 if not failures else 1)


if __name__ == "__main__":
    main()
