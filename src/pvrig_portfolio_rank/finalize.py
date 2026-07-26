#!/usr/bin/env python3
"""Create the disclosure-complete, portfolio-ranked SICBC final candidate set.

The first-pass pool is preserved so the organizer-compatible sequence checks
remain valid.  This script fixes the submission semantics:

* candidates are correctly described as optimisation of disclosed VHH parents;
* the top-ten portfolio is balanced across four parent frameworks and H3 lengths;
* score labels are calibrated as sequence/structure-composition priors rather
  than affinity predictions; and
* CDR-only versus full-length developability checks are reported separately.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis


PARENT_IDS = {
    "F1": "20",
    "F2": "30",
    "F3": "38",
    "F4": "39",
    "F5": "151",
}

# Multi-objective top-ten portfolio selected from the first-pass ranked pool.
# Quotas: F1×2, F3×4, F4×2, F5×2; only one H3 longer than 18 residues.
TOP10_OLD_RANKS = [1, 2, 3, 6, 8, 10, 11, 14, 26, 29]

CDR_RISK_PATTERNS = {
    "N糖基化": re.compile(r"N[^P][ST]"),
    "脱酰胺/异构化": re.compile(r"NS|NG|DG|DS|DP"),
    "连续疏水": re.compile(r"[AILMFWVY]{4,}"),
    "三连重复": re.compile(r"(.)\1\1"),
}

FULL_SEQUENCE_PATTERNS = {
    "N糖基化": re.compile(r"N[^P][ST]"),
    "Asn热点": re.compile(r"NS|NG"),
    "Asp热点": re.compile(r"DG|DS|DP"),
    "氧化性Met": re.compile(r"M"),
    "游离Cys": re.compile(r"C"),
    "连续疏水": re.compile(r"[AILMFWVY]{4,}"),
}


def percentile(values: list[float]) -> list[float]:
    """Average-rank percentile in [0, 100], with higher values preferred."""
    array = np.asarray(values, dtype=float)
    result = []
    for value in array:
        less = float(np.sum(array < value))
        equal = float(np.sum(array == value))
        result.append(100.0 * (less + 0.5 * equal) / len(array))
    return result


def load_parents(path: Path, parent_ids: dict[str, str]) -> dict[str, dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = {str(row["抗体名称"]): row for row in reader}
    parents = {}
    for framework, parent_id in parent_ids.items():
        row = rows[parent_id]
        parents[framework] = {
            "id": parent_id,
            "name": f"SHR-2002 VHH-{parent_id}",
            "source": (
                f"{row['来自专利']} / {row['开发公司']} / {row['相关药物管线']}"
            ),
            "sequence": row["抗体重链氨基酸"].strip().upper(),
            "antibody_type": row["类型(IgG/VHH)"],
        }
    return parents


def motif_counts(sequence: str, patterns: dict[str, re.Pattern]) -> dict[str, int]:
    return {label: len(list(pattern.finditer(sequence))) for label, pattern in patterns.items()}


def sensitivity_intervals(
    candidates: list[dict],
    *,
    weights: list[float] | None = None,
    draws_count: int = 2000,
    concentration: float = 80.0,
    seed: int = 20260726,
) -> None:
    """Attach weight-perturbation score and rank intervals.

    These are not statistical confidence intervals for affinity.  They quantify
    only how much the sequence-prior rank changes when component weights vary.
    """
    fields = [
        "esm_target_score",
        "interface_score",
        "cdr_prior_score",
        "anarcii_score",
        "developability_score",
    ]
    matrix = np.column_stack(
        [percentile([float(item[field]) for item in candidates]) for field in fields]
    )
    novelty = np.asarray(
        [100.0 * max(0.0, min(1.0, (0.75 - item["max_positive_identity"]) / 0.25))
         for item in candidates]
    )
    h3_length = np.asarray([len(item["cdr3"]) for item in candidates], dtype=float)
    length_score = 100.0 * np.exp(-((h3_length - 15.0) / 7.0) ** 2)
    matrix = np.column_stack([matrix, novelty, length_score])

    base_weights = np.asarray(
        weights or [0.24, 0.18, 0.16, 0.10, 0.14, 0.10, 0.08]
    )
    rng = np.random.default_rng(seed)
    draws = rng.dirichlet(base_weights * concentration, size=draws_count)
    scores = matrix @ draws.T
    ranks = np.empty_like(scores)
    for column in range(scores.shape[1]):
        order = np.argsort(-scores[:, column], kind="stable")
        ranks[order, column] = np.arange(1, len(candidates) + 1)

    central = matrix @ base_weights
    for index, item in enumerate(candidates):
        item["calibrated_prior_score"] = round(float(central[index]), 2)
        item["weight_sensitivity_score_p10"] = round(float(np.percentile(scores[index], 10)), 2)
        item["weight_sensitivity_score_p90"] = round(float(np.percentile(scores[index], 90)), 2)
        item["weight_sensitivity_rank_p10"] = int(math.floor(np.percentile(ranks[index], 10)))
        item["weight_sensitivity_rank_p90"] = int(math.ceil(np.percentile(ranks[index], 90)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-json", type=Path, required=True)
    parser.add_argument("--positive-csv", type=Path, required=True)
    parser.add_argument("--model-config", type=Path)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--audit-json", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.input_json.read_text(encoding="utf-8"))
    candidates = payload["candidates"]
    if args.model_config:
        model_config = json.loads(args.model_config.read_text(encoding="utf-8"))
        stage2 = model_config["stage2"]
        parent_ids = stage2["parent_ids"]
        top10_old_ranks = stage2["top10_stage1_ranks"]
        quotas = stage2["top10_framework_quotas"]
        max_long_h3 = stage2["top10_max_h3_over_18"]
        weight_map = stage2["calibrated_percentile_weights"]
        weights = list(weight_map.values())
        sensitivity_intervals(
            candidates,
            weights=weights,
            draws_count=stage2["weight_sensitivity_draws"],
            concentration=stage2["weight_dirichlet_concentration"],
            seed=model_config["random_seed"],
        )
    else:
        parent_ids = PARENT_IDS
        top10_old_ranks = TOP10_OLD_RANKS
        quotas = {"F1": 2, "F3": 4, "F4": 2, "F5": 2}
        max_long_h3 = 1
        sensitivity_intervals(candidates)
    parents = load_parents(args.positive_csv, parent_ids)

    selected = [next(item for item in candidates if item["rank"] == rank)
                for rank in top10_old_ranks]
    selected_ids = {id(item) for item in selected}
    remaining = [item for item in candidates if id(item) not in selected_ids]
    remaining.sort(key=lambda item: (-item["composite_score"], item["name"]))
    reordered = selected + remaining

    for rank, item in enumerate(reordered, start=1):
        old_rank = int(item["rank"])
        parent = parents[item["framework"]]
        digest = hashlib.sha256(item["sequence"].encode()).hexdigest()[:5].upper()
        item["old_rank"] = old_rank
        item["rank"] = rank
        item["name"] = f"PVRIG-OPT{rank:02d}-{digest}"
        item["design_type"] = "优化改造"
        item["parent"] = parent
        item["cdr_lengths"] = [len(item["cdr1"]), len(item["cdr2"]), len(item["cdr3"])]
        cdr_counts = Counter()
        for cdr in (item["cdr1"], item["cdr2"], item["cdr3"]):
            cdr_counts.update(motif_counts(cdr, CDR_RISK_PATTERNS))
        item["cdr_risk_motifs"] = {
            label: int(cdr_counts[label]) for label in CDR_RISK_PATTERNS
        }
        item["full_sequence_risk_motifs"] = motif_counts(
            item["sequence"], FULL_SEQUENCE_PATTERNS
        )
        analysis = ProteinAnalysis(item["sequence"])
        item["full_sequence_metrics"] = {
            "length": len(item["sequence"]),
            "molecular_weight_da": round(float(analysis.molecular_weight()), 1),
            "theoretical_pI": round(float(analysis.isoelectric_point()), 2),
            "gravy": round(float(analysis.gravy()), 3),
        }
        item["evidence_level"] = (
            "B：序列/接触组成先验；未进行候选-PVRIG复合物几何判别，须以表达和阻断实验确认"
        )
        item["portfolio_reason"] = (
            f"{item['framework']}/{parent['name']}；H3={len(item['cdr3'])} aa；"
            + ("前10家族-长度配额入选" if rank <= 10 else "单体先验顺序")
        )

    top10 = reordered[:10]
    framework_counts = Counter(item["framework"] for item in top10)
    h3_lengths = [len(item["cdr3"]) for item in top10]
    assert framework_counts == Counter(quotas)
    assert sum(length > 18 for length in h3_lengths) <= max_long_h3
    assert len({item["sequence"] for item in reordered}) == 50
    assert all(item["max_positive_identity"] < 0.75 for item in reordered)
    assert all(sum(item["cdr_risk_motifs"].values()) == 0 for item in reordered)

    output = {
        **{key: value for key, value in payload.items() if key != "candidates"},
        "model_name": "PVRIG-PortfolioRank v2.0",
        "design_semantics": "已披露阳性VHH母本框架上的三CDR优化改造",
        "structure_evidence": {
            "8X6B": "PVRIG–Nectin-2 complex, X-ray, 2.0 Å",
            "9E6Y": "CD112 domain 1–CD112R complex, X-ray, 2.2 Å",
            "scope": "仅用于PVRIG接触壳层与界面组成先验，不是候选抗体复合物结构",
        },
        "score_disclaimer": (
            "校准分及权重敏感性区间只用于序列先验排序；不代表Kd、IC50、阻断率、"
            "实验成功概率或候选-PVRIG复合物几何可信度。"
        ),
        "portfolio_rule": (
            "前10在单体先验基础上施加家族与H3长度配额："
            + "、".join(f"{framework}×{count}" for framework, count in quotas.items())
            + f"；H3>18 aa不超过{max_long_h3}条，以降低共同失效和长环表达风险。"
        ),
        "parents": parents,
        "candidates": reordered,
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    audit = {
        "candidate_count": len(reordered),
        "rank_contiguous": [item["rank"] for item in reordered] == list(range(1, 51)),
        "unique_names": len({item["name"] for item in reordered}),
        "unique_sequences": len({item["sequence"] for item in reordered}),
        "design_type_values": sorted({item["design_type"] for item in reordered}),
        "parent_framework_mapping": parent_ids,
        "top10_old_ranks": top10_old_ranks,
        "top10_framework_counts": dict(sorted(framework_counts.items())),
        "top10_h3_lengths": h3_lengths,
        "top10_long_h3_count": sum(length > 18 for length in h3_lengths),
        "max_positive_cdr_identity": max(item["max_positive_identity"] for item in reordered),
        "all_cdr_risk_motif_counts_zero": all(
            sum(item["cdr_risk_motifs"].values()) == 0 for item in reordered
        ),
        "full_sequence_risk_motifs_are_disclosed": True,
        "claim_scope": output["score_disclaimer"],
    }
    args.audit_json.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
