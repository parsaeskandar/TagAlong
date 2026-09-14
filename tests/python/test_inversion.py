"""Coordinate translation across an INVERTED block.

Every other coordinate-translation test in this suite works on forward-strand
homology, so nothing currently exercises the case where the two haplotypes
traverse the shared sequence in opposite directions. That matters because the
orientation handling is uneven across the stack:

  * haplotype_coverage and the Table 2 routing probe try BOTH node
    orientations, so an inverted haplotype is still discovered and scored;
  * check_common_node probes only the SOURCE's orientation
    (``Node::encode(node_id, is_rev)``), so the anchor search inside the
    coordinate trace can miss an inverted target entirely -- a node traversed
    the other way is a different ``gbwt::node_type``;
  * translate() hardcodes ``ti.strand = '+'`` at both call sites, so even a
    correctly-located inverted block is reported as forward.

The fixture is coord_translation_tests/test2_inversion, whose ground truth is
computed by hand in its README. These tests assert that ground truth. The two
marked xfail encode the defects above: they are expected to fail today and to
start passing when orientation is handled, at which point the marker should be
removed rather than the test weakened.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE = REPO_ROOT / "coord_translation_tests" / "test2_inversion"
BUILDER = REPO_ROOT / "tests" / "fixtures" / "build_fixture_indexes.sh"

# Only graph.gfa is committed -- *.gbz is gitignored, so the GBZ is derived
# too. Everything below graph.gfa is rebuilt on demand.
COORD_FILES = ("rlbwt_rindex.ri", "sampled.tags", "fastlocate.ri",
               "output.t1", "output.t2")

# Ground truth from the fixture README. Nodes 2 and 3 are shared; S1 traverses
# them forward, S2 in reverse.
SRC_CONTIG = "S1#1#chr1"
SRC_START, SRC_END = 16, 44
TGT_HAP = "S2#2"
TGT_START, TGT_END = 12, 40
# Bin/anchor granularity means the reported span can be a little wider or
# narrower than the exact block; it must not be somewhere else entirely.
TOL = 4


def _tools_available() -> bool:
    vg = os.environ.get("VG") or shutil.which("vg")
    return bool(vg) and all(shutil.which(t) for t in ("gbz_extract", "grlbwt-cli"))


@pytest.fixture(scope="module")
def inversion_index(liftover):
    """Load the inversion fixture, building its coordinate indexes if needed."""
    gbz = FIXTURE / "graph.gbz"
    if not gbz.exists():
        vg = os.environ.get("VG") or shutil.which("vg")
        if not vg:
            pytest.skip("vg not on PATH; cannot build the inversion fixture GBZ")
        proc = subprocess.run(
            [vg, "gbwt", "--gbz-format", "-g", str(gbz),
             "-G", str(FIXTURE / "graph.gfa")],
            capture_output=True, text=True,
        )
        if proc.returncode != 0 or not gbz.exists():
            pytest.skip(f"could not build {gbz.name}:\n{proc.stderr[-1000:]}")

    if not all((FIXTURE / f).exists() for f in COORD_FILES):
        if not _tools_available():
            pytest.skip(
                "inversion fixture indexes are missing and the build toolchain "
                "(vg, gbz_extract, grlbwt-cli) is not on PATH; build them with "
                f"{BUILDER} {FIXTURE}/graph.gbz {FIXTURE} --coord-only"
            )
        env = dict(os.environ)
        env.setdefault("BIN", str(REPO_ROOT / "bin"))
        env.setdefault("VG", shutil.which("vg") or "")
        env.setdefault("THREADS", "2")
        proc = subprocess.run(
            ["bash", str(BUILDER), str(FIXTURE / "graph.gbz"), str(FIXTURE),
             "--coord-only"],
            env=env, capture_output=True, text=True,
        )
        if proc.returncode != 0 or not all((FIXTURE / f).exists() for f in COORD_FILES):
            pytest.skip(f"could not build the inversion fixture:\n{proc.stderr[-2000:]}")

    idx = liftover.Index()
    idx.load(
        str(FIXTURE / "graph.gbz"),
        str(FIXTURE / "rlbwt_rindex.ri"),
        str(FIXTURE / "sampled.tags"),
        str(FIXTURE / "fastlocate.ri"),
        str(FIXTURE / "output.t1"),
        str(FIXTURE / "output.t2"),
    )
    return idx


def _describe(recs):
    return [(r.haplotype, r.start, r.end, r.strand) for r in recs]


def test_inversion_fixture_has_both_haplotypes(inversion_index):
    """Guards the fixture itself: if the graph or its metadata changes, the
    orientation tests below would skip or pass for the wrong reason."""
    names = inversion_index.get_haplotype_names()
    assert any(n.startswith("S1#1") for n in names), names
    assert any(n.startswith("S2#2") for n in names), names


@pytest.mark.xfail(
    reason="check_common_node probes only the source's node orientation, so an "
           "inverted target has no common nodes to anchor on",
    strict=False,
)
def test_inversion_is_translated_at_all(inversion_index):
    """The inverted block must translate to SOMETHING on the other haplotype."""
    res = inversion_index.translate(SRC_CONTIG, SRC_START, SRC_END, TGT_HAP)
    assert res, (
        f"{SRC_CONTIG}[{SRC_START},{SRC_END}) -> {TGT_HAP} produced no intervals; "
        "the block is shared, only its orientation differs"
    )


@pytest.mark.xfail(
    reason="same single-orientation probe; and translate() hardcodes strand '+'",
    strict=False,
)
def test_inversion_lands_on_the_right_interval(inversion_index):
    """It must land on [12,40) -- the shared block -- not on the flanks."""
    res = inversion_index.translate(SRC_CONTIG, SRC_START, SRC_END, TGT_HAP)
    assert res, "no intervals (see test_inversion_is_translated_at_all)"
    assert any(
        abs(r.start - TGT_START) <= TOL and abs(r.end - TGT_END) <= TOL
        for r in res
    ), f"expected ~[{TGT_START},{TGT_END}), got {_describe(res)}"


@pytest.mark.xfail(
    reason="translate() sets ti.strand = '+' unconditionally at both call sites "
           "(pangenome_server.cpp); no reverse strand is ever reported",
    strict=False,
)
def test_inversion_reports_minus_strand(inversion_index):
    """The whole point: an inverted correspondence must be reported as '-'."""
    res = inversion_index.translate(SRC_CONTIG, SRC_START, SRC_END, TGT_HAP)
    assert res, "no intervals (see test_inversion_is_translated_at_all)"
    assert any(r.strand == "-" for r in res), (
        f"inverted block reported as forward: {_describe(res)}"
    )


def test_inverted_haplotype_is_still_discoverable(inversion_index):
    """Routing DOES try both orientations, so the inverted haplotype must at
    least be listed as reachable even while the trace above cannot place it.

    This is the asymmetry worth pinning: a haplotype that is advertised as
    reachable but yields no coordinates is worse than one that is absent, and
    this test fails if routing ever regresses to forward-only too.
    """
    names = inversion_index.get_haplotype_names()
    s2 = next((n for n in names if n.startswith("S2#2")), None)
    assert s2, names
    try:
        reachable = inversion_index.translatable_haplotypes(
            SRC_CONTIG, SRC_START, SRC_END)
    except Exception as exc:  # pragma: no cover - extension without the call
        pytest.skip(f"translatable_haplotypes unavailable: {exc}")
    assert any(h.startswith("S2#2") for h in reachable), (
        f"inverted haplotype not discoverable at all: {reachable}"
    )
