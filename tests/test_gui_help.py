"""The contextual help on the GUI fields.

Streamlit's ``help=`` renders a "?" next to a field and shows the text on hover. The texts live in
``poa/gui/help_texts.py`` so they can be reviewed as documentation; these tests keep that file and
the forms from drifting apart — an unused text is wording nobody sees, and a field that lost its
tooltip is guidance silently removed.

The list of fields below is deliberately explicit rather than "every widget": trivial fields are
meant to have no tooltip, so a count would either force noise onto them or pass vacuously.
"""
import ast
import re
from pathlib import Path

import pytest

from poa.gui.help_texts import HELP

APP = Path(__file__).resolve().parents[1] / "poa" / "gui" / "app.py"
SOURCE = APP.read_text(encoding="utf-8")

#: fields whose meaning is not evident from the label, so losing the tooltip is a regression
FIELDS_THAT_MUST_EXPLAIN_THEMSELVES = [
    "proteins_fasta", "swap_headers", "len_min", "len_max", "mhla", "mic",
    "mhcii_alleles", "mhcii_length", "mhcii_upload", "b2_json", "map_specie", "map_protein",
    "cons_threshold", "cons_operator", "cons_mode", "comparison_set", "cons_csv_upload",
    "cons_declared_threshold",
    "poa2_objective", "poa2_rf", "poa2_imin", "poa2_imax", "poa2_m", "poa2_idf",
    "poa2_declared_t", "poa2_declared_op",
    "viz_pdb", "viz_style", "viz_surface", "viz_chain",
]


@pytest.mark.parametrize("field", FIELDS_THAT_MUST_EXPLAIN_THEMSELVES)
def test_field_has_help_text_and_it_is_wired_up(field):
    """Each listed field needs a text, and the app has to actually pass it."""
    assert field in HELP, f"falta o texto de ajuda de '{field}'"
    assert f'HELP["{field}"]' in SOURCE, f"'{field}' tem texto mas não é usado em app.py"


def test_no_help_text_is_orphaned():
    """
    Every entry must reach a field.

    ``_manual_upload`` looks its text up dynamically as ``HELP.get(f"upload_{key}")``, so the
    per-method keys are matched against the method argument instead of a literal subscript.
    """
    methods = set(re.findall(r'_manual_upload\(\s*"[^"]*",\s*"(\w+)"', SOURCE))
    referenced = set(re.findall(r'HELP\["(\w+)"\]', SOURCE)) | {f"upload_{m}" for m in methods}

    assert not (set(HELP) - referenced)


def test_install_guides_cover_every_required_package():
    """
    The hand-written install lists must not drift from ``requirements.txt``.

    ``TESTING_WSL.md`` installs a list of packages by hand instead of using ``-r`` (it skips the
    optional ones), and ``COMO_TESTAR.md`` verifies the environment with an import line. Both are
    copies of the dependency list, and a package added to ``requirements.txt`` but not to them is
    a package the reader never learns they need — which is how ``py3Dmol`` came to look like a
    separate manual install when it never was one.
    """
    repo = Path(__file__).resolve().parents[1]
    optional = {"playwright", "markdown", "pytest"}      # documented as extras, installed on demand
    required = set()
    for line in (repo / "requirements.txt").read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        name = re.split(r"[<>=!\[]", line, maxsplit=1)[0].strip()
        if name and name not in optional:
            required.add(name)

    wsl = (repo / "TESTING_WSL.md").read_text(encoding="utf-8")
    missing = {pkg for pkg in required if pkg.lower() not in wsl.lower()}

    assert not missing, f"pacotes ausentes do guia de instalação da WSL: {sorted(missing)}"


def test_every_help_argument_comes_from_the_shared_table():
    """
    No tooltip text inlined in the forms.

    The point of the table is that the wording is reviewable in one place; a string written
    straight into a widget call escapes that, and escapes these tests with it.
    """
    tree = ast.parse(SOURCE)
    inlined = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        for keyword in node.keywords:
            if keyword.arg != "help":
                continue
            # allowed: HELP["x"] and HELP.get(...); the sidebar buttons are UI actions, not fields
            if isinstance(keyword.value, ast.Constant) and getattr(node.func, "attr", "") != "button":
                inlined.append(ast.unparse(keyword.value)[:60])

    assert not inlined, f"ajuda escrita direto no widget: {inlined}"
