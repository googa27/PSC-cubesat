from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


def test_portfolio_architecture_contract() -> None:
    result = subprocess.run(
        [sys.executable, "scripts/check_portfolio_architecture.py"],
        check=False,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_origin_and_upstream_remote_claims_are_truthful() -> None:
    root = Path(__file__).resolve().parents[2]
    contract = json.loads((root / "docs" / "ARCHITECTURE.yaml").read_text(encoding="utf-8"))
    assert contract["repository"]["origin"] == "https://github.com/googa27/PSC-cubesat.git"
    assert contract["repository"]["upstream_remote"] == "none configured locally"

    readme = (root / "README.md").read_text(encoding="utf-8")
    assert "Configured canonical origin: https://github.com/googa27/PSC-cubesat.git" in readme
    assert "No separate `upstream` remote is configured locally" in readme
