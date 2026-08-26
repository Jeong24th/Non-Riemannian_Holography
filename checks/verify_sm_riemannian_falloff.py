#!/usr/bin/env python3
"""Compatibility regression for the Riemannian fixed-frame falloff.

The former saddle-adapted Riemannian system has been replaced by the
single fixed-frame fluctuation h_{p qbar} used for both exact saddles.
This entry point preserves the historical command name while checking
the new LaTeX contract and running the exact symbolic projection test.
"""

from __future__ import annotations

import re
import runpy
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEX_PATH = ROOT / "NR_Holography.tex"
EXACT_CHECK = ROOT / "checks" / "verify_exact_projected_fluctuations.py"


def check_tex_contract() -> bool:
    if not TEX_PATH.is_file():
        print("LaTeX contract: SKIP (NR_Holography.tex is not in the archive)")
        return False
    compact = re.sub(r"\s+", "", TEX_PATH.read_text(encoding="utf-8"))
    required = {
        "exact R projection": r"\label{SMexactprojectionR}",
        "R fixed type-changing coefficient": (
            r"h_{\ominus\bar\oplus}^{(2)}=2"
        ),
        "R plus chiral coefficient": (
            r"h_{\oplus\bar\oplus}^{(2)}=2L_{+}"
        ),
        "R minus chiral coefficient": (
            r"h_{\ominus\bar\ominus}^{(2)}=2L_{-}"
        ),
        "R composite upper-right coefficient": (
            r"h_{\oplus\bar\ominus}^{(2)}=2L_{+}L_{-}"
        ),
        "common linearized system": r"\label{SMNRlin}",
    }
    missing = [name for name, fragment in required.items() if fragment not in compact]
    if missing:
        raise AssertionError(
            "NR_Holography.tex is out of sync with the fixed-frame R falloff: "
            + ", ".join(missing)
        )

    forbidden = {
        "obsolete saddle-adapted system": r"\label{SMlinEDFE}",
        "obsolete mixed falloff": r"\label{SMfalloff}",
        "obsolete h-plus/minus shorthand": r"h_{\pm}:=",
    }
    present = [name for name, fragment in forbidden.items() if fragment in compact]
    if present:
        raise AssertionError(
            "NR_Holography.tex still contains obsolete R falloff material: "
            + ", ".join(present)
        )
    return True


def main() -> None:
    if check_tex_contract():
        print("Fixed-frame R LaTeX contract: PASS")
    runpy.run_path(str(EXACT_CHECK), run_name="__main__")
    print("R fixed-frame falloff regression: PASS")


if __name__ == "__main__":
    main()
