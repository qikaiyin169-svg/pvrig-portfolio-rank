from __future__ import annotations

import ast
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class SourceTests(unittest.TestCase):
    def test_all_python_files_parse(self) -> None:
        paths = list((ROOT / "src").rglob("*.py"))
        paths += list((ROOT / "scripts").glob("*.py"))
        paths += list((ROOT / "tests").glob("*.py"))
        self.assertGreaterEqual(len(paths), 8)
        for path in paths:
            with self.subTest(path=path.relative_to(ROOT)):
                ast.parse(path.read_text(encoding="utf-8"), filename=str(path))

    def test_model_artifact_has_pinned_inputs(self) -> None:
        model = json.loads(
            (ROOT / "models/pvrig_portfolio_rank_v2.json").read_text(encoding="utf-8")
        )
        self.assertEqual(len(model["inputs"]["organizer_validator_commit"]), 40)
        self.assertEqual(len(model["inputs"]["pdb_8x6b_sha256"]), 64)
        self.assertEqual(len(model["inputs"]["pdb_9e6y_sha256"]), 64)


if __name__ == "__main__":
    unittest.main()
