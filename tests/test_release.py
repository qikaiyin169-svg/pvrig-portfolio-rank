from __future__ import annotations

import json
import unittest
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class ReleaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = json.loads(
            (ROOT / "models/pvrig_portfolio_rank_v2.json").read_text(encoding="utf-8")
        )
        cls.result = json.loads(
            (ROOT / "results/final_candidates.json").read_text(encoding="utf-8")
        )
        cls.audit = json.loads(
            (ROOT / "results/final_audit.json").read_text(encoding="utf-8")
        )
        cls.validator = json.loads(
            (ROOT / "results/validator_report.json").read_text(encoding="utf-8")
        )

    def test_model_weights_sum_to_one(self) -> None:
        self.assertAlmostEqual(
            sum(self.model["stage1"]["preliminary_weights"].values()), 1.0
        )
        self.assertAlmostEqual(
            sum(self.model["stage1"]["composite_weights"].values()), 1.0
        )
        self.assertAlmostEqual(
            sum(self.model["stage2"]["calibrated_percentile_weights"].values()), 1.0
        )

    def test_candidate_count_rank_and_uniqueness(self) -> None:
        candidates = self.result["candidates"]
        self.assertEqual(len(candidates), 50)
        self.assertEqual([item["rank"] for item in candidates], list(range(1, 51)))
        self.assertEqual(len({item["name"] for item in candidates}), 50)
        self.assertEqual(len({item["sequence"] for item in candidates}), 50)

    def test_sequence_alphabet_and_cdr_threshold(self) -> None:
        amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
        for item in self.result["candidates"]:
            self.assertFalse(set(item["sequence"]) - amino_acids)
            self.assertLess(
                item["max_positive_identity"],
                self.model["stage2"]["validation_identity_threshold"],
            )

    def test_top10_portfolio_constraints(self) -> None:
        top10 = self.result["candidates"][:10]
        counts = Counter(item["framework"] for item in top10)
        self.assertEqual(
            dict(sorted(counts.items())),
            self.model["stage2"]["top10_framework_quotas"],
        )
        self.assertLessEqual(
            sum(len(item["cdr3"]) > 18 for item in top10),
            self.model["stage2"]["top10_max_h3_over_18"],
        )

    def test_validator_release_result(self) -> None:
        self.assertEqual(self.validator["candidate_count"], 50)
        self.assertEqual(self.validator["passed_count"], 50)
        self.assertEqual(self.validator["failure_count"], 0)
        self.assertEqual(self.validator["failures"], [])
        self.assertAlmostEqual(
            self.audit["max_positive_cdr_identity"], 0.7142857142857143
        )


if __name__ == "__main__":
    unittest.main()
