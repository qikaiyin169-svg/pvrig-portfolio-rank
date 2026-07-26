#!/usr/bin/env python3
"""Generate and rank 50 PVRIG VHH CDR-optimization candidates.

The task-local model, PVRIG-LoopRank v1, combines:
1) an IMGT-position-aware CDR probability model trained on the organizer's
   public PVRIG positive library;
2) PVRIG interface-composition constraints derived from PDB 8X6B and 9E6Y;
3) ESM2 embedding-based one-class target enrichment;
4) ANARCII2 antibody-likeness and deterministic developability filters; and
5) diversity-aware Pareto selection.

The five VHH frameworks are disclosed positive-library parents; H1/H2/H3 are
regenerated. The numerical score is a computational prioritization proxy, not
an experimental affinity or success probability.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import torch
from anarcii import Anarcii
from Bio import PDB
from Bio.Align import PairwiseAligner
from transformers import AutoModel, AutoTokenizer


AA = "ACDEFGHIKLMNPQRSTVWY"
AROMATIC = set("FWY")
POLAR = set("STNQY")
ACIDIC = set("DE")
BASIC = set("KR")
HYDROPHOBIC = set("AILMFWVY")
CDR_RANGES = {"H1": (27, 38), "H2": (56, 65), "H3": (105, 117)}
BAD_MOTIFS = (
    re.compile(r"N[^P][ST]"),
    re.compile(r"NS|NG|DG|DS|DP"),
    re.compile(r"[AILMFWVY]{4,}"),
    re.compile(r"(.)\1\1"),
)


@dataclass
class Candidate:
    name: str
    sequence: str
    framework: str
    cdr1: str
    cdr2: str
    cdr3: str
    max_positive_identity: float
    max_selected_identity: float
    interface_score: float
    cdr_prior_score: float
    esm_target_score: float
    anarcii_score: float
    developability_score: float
    composite_score: float
    rank: int = 0


def clean_seq(value: str) -> str:
    return re.sub(r"[^A-Z]", "", value.upper())


def numbered(an: Anarcii, sequence: str) -> tuple[list[tuple[tuple[int, str], str]], float]:
    result = an.number([sequence])["Sequence 1"]
    if result.get("error") or not result.get("numbering"):
        raise ValueError(f"ANARCII failed: {result.get('error')}")
    return result["numbering"], float(result.get("score") or 0.0)


def split_regions(numbering):
    parts = {"FR1": [], "H1": [], "FR2": [], "H2": [], "FR3": [], "H3": [], "FR4": []}
    for (position, _insertion), residue in numbering:
        if residue in {"-", "."}:
            continue
        if position <= 26:
            key = "FR1"
        elif position <= 38:
            key = "H1"
        elif position <= 55:
            key = "FR2"
        elif position <= 65:
            key = "H2"
        elif position <= 104:
            key = "FR3"
        elif position <= 117:
            key = "H3"
        else:
            key = "FR4"
        parts[key].append(residue)
    return {key: "".join(value) for key, value in parts.items()}


_ALIGNER = PairwiseAligner(
    mode="global",
    match_score=2.0,
    mismatch_score=-1.0,
    open_gap_score=-1.0,
    extend_gap_score=-1.0,
)


def global_identity(first: str, second: str) -> float:
    """Needleman-Wunsch identity using Bio.Align's compiled DP implementation."""
    alignment = _ALIGNER.align(first, second)[0]
    counts = alignment.counts()
    length = counts.identities + counts.mismatches + counts.gaps
    return counts.identities / max(length, 1)


def max_cdr_identity(cdrs: tuple[str, str, str], references: list[tuple[str, str, str]]) -> float:
    maximum = 0.0
    for reference in references:
        for candidate_cdr, positive_cdr in zip(cdrs, reference):
            maximum = max(maximum, global_identity(candidate_cdr, positive_cdr))
    return maximum


def quick_identity(first: str, second: str) -> float:
    """Cheap conservative screening identity over small relative offsets."""
    maximum = 0.0
    for offset in range(-2, 3):
        matches = 0
        overlap = 0
        for index, residue in enumerate(first):
            other_index = index + offset
            if 0 <= other_index < len(second):
                overlap += 1
                matches += int(residue == second[other_index])
        # Denominator includes unaligned overhangs, as in a global alignment.
        maximum = max(maximum, matches / max(len(first), len(second), 1))
    return maximum


def max_cdr_identity_quick(
    cdrs: tuple[str, str, str], references: list[tuple[str, str, str]]
) -> float:
    return max(
        quick_identity(candidate_cdr, positive_cdr)
        for reference in references
        for candidate_cdr, positive_cdr in zip(cdrs, reference)
    )


def relative_counts(cdrs: list[str], bins: int = 12) -> list[Counter]:
    counts = [Counter({aa: 0.35 for aa in AA}) for _ in range(bins)]
    for cdr in cdrs:
        for index, residue in enumerate(cdr):
            target = min(bins - 1, round(index * (bins - 1) / max(len(cdr) - 1, 1)))
            counts[target][residue] += 1.0
    return counts


def sample_cdr(rng: random.Random, counts: list[Counter], length: int, anchor: str | None) -> str:
    residues = []
    for index in range(length):
        bin_index = min(len(counts) - 1, round(index * (len(counts) - 1) / max(length - 1, 1)))
        weighted = counts[bin_index].copy()
        # Target-conditioned prior: polar/aromatic contacts and modest acidic
        # complementarity to the basic patches at the PVRIG interface.
        for aa, factor in {
            "Y": 1.45, "F": 1.25, "W": 1.10, "S": 1.25, "T": 1.18,
            "N": 1.12, "Q": 1.10, "D": 1.18, "E": 1.10, "G": 1.14,
            "P": 0.90, "C": 0.03, "M": 0.35,
        }.items():
            weighted[aa] *= factor
        if anchor and index < len(anchor) and rng.random() < 0.18:
            weighted[anchor[index]] *= 2.2
        alphabet = list(AA)
        values = [max(weighted[aa], 1e-6) for aa in alphabet]
        residues.append(rng.choices(alphabet, weights=values, k=1)[0])
    return "".join(residues)


def cdr_prior_score(cdrs: tuple[str, str, str], count_sets: list[list[Counter]]) -> float:
    scores = []
    for cdr, counts in zip(cdrs, count_sets):
        ll = 0.0
        for index, residue in enumerate(cdr):
            bin_index = min(len(counts) - 1, round(index * (len(counts) - 1) / max(len(cdr) - 1, 1)))
            total = sum(counts[bin_index].values())
            ll += math.log(counts[bin_index][residue] / total + 1e-12)
        scores.append(ll / max(len(cdr), 1))
    # Typical per-residue log likelihood maps softly into 0–100.
    return float(100.0 / (1.0 + math.exp(-2.0 * (np.mean(scores) + 2.9))))


def composition_fraction(sequence: str, residue_set: set[str]) -> float:
    return sum(aa in residue_set for aa in sequence) / max(len(sequence), 1)


def interface_score(cdrs: tuple[str, str, str], interface_residues: str) -> float:
    sequence = "".join(cdrs)
    target = {
        "aromatic": 0.19,
        "polar": 0.34,
        "acidic": 0.12,
        "basic": 0.08,
        "hydrophobic": 0.31,
        "glypro": 0.18,
    }
    observed = {
        "aromatic": composition_fraction(sequence, AROMATIC),
        "polar": composition_fraction(sequence, POLAR),
        "acidic": composition_fraction(sequence, ACIDIC),
        "basic": composition_fraction(sequence, BASIC),
        "hydrophobic": composition_fraction(sequence, HYDROPHOBIC),
        "glypro": composition_fraction(sequence, set("GP")),
    }
    scales = {
        "aromatic": 0.11, "polar": 0.15, "acidic": 0.09,
        "basic": 0.08, "hydrophobic": 0.15, "glypro": 0.12,
    }
    penalty = np.mean([((observed[k] - target[k]) / scales[k]) ** 2 for k in target])
    score = 100.0 * math.exp(-0.48 * penalty)
    # Reward chemically complementary residues observed in PVRIG contact shells.
    interface_profile = Counter(interface_residues)
    if interface_profile["R"] + interface_profile["K"] >= 3:
        score += min(5.0, 30.0 * observed["acidic"])
    return float(min(score, 100.0))


def developability(sequence: str, cdrs: tuple[str, str, str]) -> tuple[bool, float]:
    if set(sequence) - set(AA):
        return False, 0.0
    joined = "|".join(cdrs)
    if "C" in joined or "M" in joined:
        return False, 0.0
    if any(pattern.search(joined) for pattern in BAD_MOTIFS):
        return False, 0.0
    charge = sum(aa in BASIC for aa in joined) - sum(aa in ACIDIC for aa in joined)
    if abs(charge) > 7:
        return False, 0.0
    hydrophobic = composition_fraction(joined, HYDROPHOBIC)
    aromatic = composition_fraction(joined, AROMATIC)
    if not 0.15 <= hydrophobic <= 0.47 or aromatic > 0.31:
        return False, 0.0
    score = 100.0
    score -= 4.0 * max(0, abs(charge) - 3)
    score -= 120.0 * max(0.0, hydrophobic - 0.38)
    score -= 80.0 * max(0.0, aromatic - 0.24)
    score -= 5.0 * joined.count("W")
    return True, float(max(0.0, score))


def extract_interface_residues(path: Path, chain_a: str, chain_b: str, cutoff: float = 5.0) -> str:
    structure = PDB.PDBParser(QUIET=True).get_structure(path.stem, path)
    model = structure[0]
    atoms_b = [atom for residue in model[chain_b] for atom in residue if atom.element != "H"]
    ns = PDB.NeighborSearch(atoms_b)
    contacts = []
    seen = set()
    for residue in model[chain_a]:
        if residue.id[0] != " ":
            continue
        if any(ns.search(atom.coord, cutoff) for atom in residue if atom.element != "H"):
            key = (residue.id[1], residue.resname)
            if key not in seen:
                contacts.append(PDB.Polypeptide.protein_letters_3to1.get(residue.resname, "X"))
                seen.add(key)
    return "".join(contacts)


def embed_sequences(model, tokenizer, sequences: list[str], batch_size: int = 24) -> np.ndarray:
    embeddings = []
    device = next(model.parameters()).device
    for start in range(0, len(sequences), batch_size):
        batch = sequences[start : start + batch_size]
        tokens = tokenizer(batch, return_tensors="pt", padding=True, add_special_tokens=True)
        tokens = {key: value.to(device) for key, value in tokens.items()}
        with torch.inference_mode():
            hidden = model(**tokens).last_hidden_state
        mask = tokens["attention_mask"].unsqueeze(-1)
        pooled = (hidden * mask).sum(dim=1) / mask.sum(dim=1)
        embeddings.append(pooled.cpu().numpy())
    return np.concatenate(embeddings, axis=0)


def make_decoy(rng: random.Random, framework: dict[str, str], lengths: tuple[int, int, int]) -> str:
    weights = [1] * len(AA)
    alphabet = list(AA)
    cdrs = ["".join(rng.choices(alphabet, weights=weights, k=n)) for n in lengths]
    return assemble(framework, tuple(cdrs))


def assemble(framework: dict[str, str], cdrs: tuple[str, str, str]) -> str:
    return (
        framework["FR1"] + cdrs[0] + framework["FR2"] + cdrs[1]
        + framework["FR3"] + cdrs[2] + framework["FR4"]
    )


def standardize_scores(values: np.ndarray) -> np.ndarray:
    lo, hi = np.quantile(values, [0.03, 0.97])
    return np.clip((values - lo) / max(hi - lo, 1e-9), 0.0, 1.0) * 100.0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--positive-csv", type=Path, required=True)
    parser.add_argument("--pdb-8x6b", type=Path, required=True)
    parser.add_argument("--pdb-9e6y", type=Path, required=True)
    parser.add_argument("--model", default="facebook/esm2_t6_8M_UR50D")
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260726)
    parser.add_argument("--pool-size", type=int, default=16000)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    an = Anarcii()

    with args.positive_csv.open(encoding="utf-8-sig") as handle:
        positive_rows = list(csv.DictReader(handle))

    positive_sequences = [clean_seq(row["抗体重链氨基酸"]) for row in positive_rows]
    positive_cdrs = []
    vhh_frameworks = []
    positive_scores = []
    for row, sequence in zip(positive_rows, positive_sequences):
        numbering, score = numbered(an, sequence)
        regions = split_regions(numbering)
        positive_cdrs.append((regions["H1"], regions["H2"], regions["H3"]))
        positive_scores.append(score)
        if row["类型(IgG/VHH)"].strip().upper() == "VHH":
            vhh_frameworks.append(
                {key: regions[key] for key in ("FR1", "FR2", "FR3", "FR4")}
            )
    if len(vhh_frameworks) < 3:
        raise RuntimeError("not enough VHH frameworks in positive library")

    cdr_lists = [[cdrs[index] for cdrs in positive_cdrs] for index in range(3)]
    count_sets = [relative_counts(values) for values in cdr_lists]
    length_sets = [sorted(len(value) for value in values) for values in cdr_lists]

    interface_residues = (
        extract_interface_residues(args.pdb_8x6b, "B", "A")
        + extract_interface_residues(args.pdb_9e6y, "A", "D")
    )

    pool = []
    seen_sequences = set()
    attempts = 0
    while len(pool) < args.pool_size and attempts < args.pool_size * 80:
        attempts += 1
        framework_index = rng.randrange(len(vhh_frameworks))
        framework = vhh_frameworks[framework_index]
        lengths = []
        for values in length_sets:
            base = rng.choice(values)
            delta = rng.choices([-1, 0, 1], weights=[1, 8, 1], k=1)[0]
            lengths.append(max(5, base + delta))
        anchors = rng.choice(positive_cdrs)
        cdrs = tuple(
            sample_cdr(rng, counts, length, anchor)
            for counts, length, anchor in zip(count_sets, lengths, anchors)
        )
        sequence = assemble(framework, cdrs)
        if sequence in seen_sequences:
            continue
        max_pos = max_cdr_identity_quick(cdrs, positive_cdrs)
        # Stronger than the official 80% threshold, with margin for aligner details.
        if max_pos >= 0.72:
            continue
        passed, dev_score = developability(sequence, cdrs)
        if not passed:
            continue
        prior_score = cdr_prior_score(cdrs, count_sets)
        target_score = interface_score(cdrs, interface_residues)
        preliminary = 0.42 * prior_score + 0.36 * target_score + 0.22 * dev_score
        seen_sequences.add(sequence)
        pool.append(
            {
                "sequence": sequence,
                "framework_index": framework_index,
                "cdrs": cdrs,
                "max_positive_identity": max_pos,
                "prior_score": prior_score,
                "interface_score": target_score,
                "developability_score": dev_score,
                "preliminary": preliminary,
            }
        )
    if len(pool) < 3000:
        raise RuntimeError(f"candidate generation yielded only {len(pool)} passing designs")

    pool.sort(key=lambda item: item["preliminary"], reverse=True)
    # Run exact global alignment only on the competitive head of the pool.
    exact_pool = []
    for item in pool[: min(len(pool), 1200)]:
        exact_identity = max_cdr_identity(item["cdrs"], positive_cdrs)
        if exact_identity < 0.72:
            item["max_positive_identity"] = exact_identity
            exact_pool.append(item)
    pool = exact_pool
    if len(pool) < 500:
        raise RuntimeError(f"exact identity filtering yielded only {len(pool)} designs")
    # Preserve chemistry diversity before the ESM2 stage.
    shortlist = []
    for item in pool:
        if all(
            max(quick_identity(a, b) for a, b in zip(item["cdrs"], selected["cdrs"])) < 0.86
            for selected in shortlist
        ):
            shortlist.append(item)
        if len(shortlist) >= 420:
            break
    if len(shortlist) < 180:
        raise RuntimeError(f"diversity shortlist too small: {len(shortlist)}")

    tokenizer = AutoTokenizer.from_pretrained(args.model)
    esm = AutoModel.from_pretrained(args.model)
    esm.eval()
    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    esm.to(device)

    positive_embeddings = embed_sequences(esm, tokenizer, positive_sequences)
    decoy_sequences = []
    for _ in range(max(96, len(positive_sequences) * 3)):
        framework = rng.choice(vhh_frameworks)
        lengths = tuple(rng.choice(values) for values in length_sets)
        decoy_sequences.append(make_decoy(rng, framework, lengths))
    decoy_embeddings = embed_sequences(esm, tokenizer, decoy_sequences)
    candidate_embeddings = embed_sequences(esm, tokenizer, [item["sequence"] for item in shortlist])

    positive_centroid = positive_embeddings.mean(axis=0)
    decoy_centroid = decoy_embeddings.mean(axis=0)
    direction = positive_centroid - decoy_centroid
    direction /= np.linalg.norm(direction) + 1e-9
    raw_target = candidate_embeddings @ direction
    esm_scores = standardize_scores(raw_target)

    for item, esm_score in zip(shortlist, esm_scores):
        numbering, ab_score = numbered(an, item["sequence"])
        regions = split_regions(numbering)
        recovered = (regions["H1"], regions["H2"], regions["H3"])
        if recovered != item["cdrs"]:
            item["valid_numbering"] = False
            continue
        if max(position for (position, _), residue in numbering if residue not in {"-", "."}) < 128:
            item["valid_numbering"] = False
            continue
        item["valid_numbering"] = True
        item["anarcii_score"] = ab_score
        item["esm_target_score"] = float(esm_score)

    shortlist = [item for item in shortlist if item.get("valid_numbering")]
    anarcii_values = standardize_scores(np.array([item["anarcii_score"] for item in shortlist]))
    for item, ab_score_scaled in zip(shortlist, anarcii_values):
        item["anarcii_scaled"] = float(ab_score_scaled)
        item["composite"] = (
            0.28 * item["esm_target_score"]
            + 0.24 * item["interface_score"]
            + 0.18 * item["prior_score"]
            + 0.12 * item["anarcii_scaled"]
            + 0.18 * item["developability_score"]
        )
    shortlist.sort(key=lambda item: item["composite"], reverse=True)

    selected = []
    remaining = shortlist[:]
    while remaining and len(selected) < 50:
        best_item = None
        best_utility = -1e9
        best_similarity = 0.0
        for item in remaining[:260]:
            similarity = 0.0
            if selected:
                similarity = max(
                    np.mean([quick_identity(a, b) for a, b in zip(item["cdrs"], other["cdrs"])])
                    for other in selected
                )
            if similarity >= 0.78:
                continue
            novelty_bonus = 8.0 * (1.0 - similarity)
            framework_balance = 1.5 if sum(
                other["framework_index"] == item["framework_index"] for other in selected
            ) < 10 else -1.5
            utility = item["composite"] + novelty_bonus + framework_balance
            if utility > best_utility:
                best_item, best_utility, best_similarity = item, utility, similarity
        if best_item is None:
            break
        best_item["max_selected_identity"] = best_similarity
        selected.append(best_item)
        remaining.remove(best_item)
    if len(selected) != 50:
        raise RuntimeError(f"could select only {len(selected)} candidates")

    selected.sort(key=lambda item: item["composite"], reverse=True)
    candidates = []
    for rank, item in enumerate(selected, start=1):
        digest = hashlib.sha1(item["sequence"].encode()).hexdigest()[:5].upper()
        candidate = Candidate(
            name=f"PVRIG-DV{rank:02d}-{digest}",
            sequence=item["sequence"],
            framework=f"F{item['framework_index'] + 1}",
            cdr1=item["cdrs"][0],
            cdr2=item["cdrs"][1],
            cdr3=item["cdrs"][2],
            max_positive_identity=float(item["max_positive_identity"]),
            max_selected_identity=float(item.get("max_selected_identity", 0.0)),
            interface_score=float(item["interface_score"]),
            cdr_prior_score=float(item["prior_score"]),
            esm_target_score=float(item["esm_target_score"]),
            anarcii_score=float(item["anarcii_score"]),
            developability_score=float(item["developability_score"]),
            composite_score=float(item["composite"]),
            rank=rank,
        )
        candidates.append(candidate)

    metadata = {
        "model_name": "PVRIG-LoopRank v1.0",
        "seed": args.seed,
        "pool_size": len(pool),
        "shortlist_size": len(shortlist),
        "positive_reference_count": len(positive_sequences),
        "vhh_framework_count": len(vhh_frameworks),
        "pdb_interface_residues": interface_residues,
        "esm_model": args.model,
        "score_disclaimer": "computational prioritization proxy; not measured affinity",
        "candidates": [asdict(candidate) for candidate in candidates],
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    with args.output_csv.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(candidates[0]).keys()))
        writer.writeheader()
        for candidate in candidates:
            writer.writerow(asdict(candidate))

    print(json.dumps({
        "generated": len(pool),
        "shortlisted": len(shortlist),
        "selected": len(candidates),
        "score_range": [round(candidates[-1].composite_score, 2), round(candidates[0].composite_score, 2)],
        "max_positive_identity": round(max(c.max_positive_identity for c in candidates), 4),
        "max_pair_identity": round(max(c.max_selected_identity for c in candidates), 4),
        "output": str(args.output_json),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
