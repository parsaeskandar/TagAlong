"""Coordinate translation across INVERTED blocks.

Every other coordinate-translation test in this suite works on forward-strand
homology. These cover the case where the two haplotypes traverse shared sequence
in opposite directions, which the stack handles in three distinct places:

  * discovery (haplotype_coverage, the Table 2 routing probe) is
    orientation-BLIND -- ``probe_tag`` folds both node orientations together via
    ``seqId(v) / 2``, and Table 2 stores no strand at all;
  * the trace is orientation-SPECIFIC -- a node traversed the other way round is
    a different ``gbwt::node_type``, so an inverted block shares nothing with the
    target's forward sequence and only lines up against its REVERSE sequence
    (gbwt seq ``2*pid+1``), where every node orientation is flipped;
  * reporting must convert reverse-sequence offsets back to the target's forward
    frame and mark the result strand '-'.

Two fixtures, because they fail differently:

  test2_inversion          the query is CONTAINED in the inverted block, so the
                           forward trace finds nothing anywhere and any
                           "try the other strand" fallback rescues it.
  test3_inversion_boundary the query SPANS the inversion boundary: collinear
                           flanks around an inverted core. Forward succeeds on
                           the flanks, so a fallback gated on "forward found
                           nothing" never fires and the core is dropped
                           silently -- a clean-looking forward translation with
                           a hole in it. This is the common query shape, since
                           callers request regions rather than exact inversion
                           extents.

NOTE on what translate() returns: despite the name, ``TranslatedInterval`` from
translate() is a per-base (source, target) PAIR -- ``start`` is the source
offset and ``end`` is the target offset (see the emit loop in
pangenome_server.cpp). Asserting on pairs is far tighter than asserting on a
span, so that is what these tests do.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
BUILDER = REPO_ROOT / "tests" / "fixtures" / "build_fixture_indexes.sh"

# Graphs are embedded rather than stored as files because .gitignore ignores both
# *.gfa and *.gbz, so a fixture kept on disk could never be committed and would
# silently be absent everywhere except the machine that created it.

# --- Fixture A: query contained in the inverted block ------------------------
#   S1#1#chr1  >1>2>3>4    nodes 1,2,3,4 forward
#   S2#2#chr1  >5<3<2>6    nodes 3,2 REVERSED, flanked by private 5 and 6
# Lengths 1=16 2=12 3=16 4=13 5=12 6=12 give:
#   S1:  1 [0,16)  2 [16,28)  3 [28,44)  4 [44,57)     total 57
#   S2:  5 [0,12)  3 [12,28)  2 [28,40)  6 [40,52)     total 52
# Nodes 1 and 4 are private to S1, so NOTHING matches on the forward strand.
GFA_CONTAINED = """H\tVN:Z:1.1
S\t1\tAAAACCCCGGGGTTTT
S\t2\tACGTACGTACGT
S\t3\tTTTTGGGGCCCCAAAA
S\t4\tGATTACAGATTAC
S\t5\tCCCCCCCCCCCC
S\t6\tGGGGGGGGGGGG
W\tS1\t1\tchr1\t0\t57\t>1>2>3>4
W\tS2\t2\tchr1\t0\t52\t>5<3<2>6
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t5\t+\t3\t-\t0M
L\t3\t-\t2\t-\t0M
L\t2\t-\t6\t+\t0M
"""

# --- Fixture B: query spans the inversion boundary ---------------------------
#   S1#1#chr1  >1>2>3>4    nodes 1,2,3,4 forward
#   S2#2#chr1  >1<3<2>4    SAME flanks 1 and 4, core 2,3 REVERSED
# Lengths 1=16 2=12 3=16 4=13 give (both haplotypes 57bp):
#   S1:  1 [0,16)  2 [16,28)  3 [28,44)  4 [44,57)
#   S2:  1 [0,16)  3 [16,32)  2 [32,44)  4 [44,57)
# So a query over the whole contig must produce THREE blocks: '+' flank,
# '-' core, '+' flank.
#
# Sequences are deliberately NOT reverse-complement palindromes here (several in
# fixture A are), so a spurious orientation match cannot pass for a real one.
GFA_BOUNDARY = """H\tVN:Z:1.1
S\t1\tACGTACGTTTGGCCAA
S\t2\tGGATTCCAGTAC
S\t3\tTTACCGGATCAGGTTA
S\t4\tGATTACAGATTAC
W\tS1\t1\tchr1\t0\t57\t>1>2>3>4
W\tS2\t2\tchr1\t0\t57\t>1<3<2>4
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t1\t+\t3\t-\t0M
L\t3\t-\t2\t-\t0M
L\t2\t-\t4\t+\t0M
"""


# --- Control fixture: NO rearrangement at all --------------------------------
#   S1#1#chr1  >1>2>3>4
#   S2#2#chr1  >1>2>3>4    the SAME walk, as a second haplotype
# Every base is shared in the same orientation, so translation must be the
# identity: src == tgt everywhere, strand '+'. This exists to tell two very
# different failures apart. If the inversion tests fail but this passes, the
# orientation handling is at fault. If this fails too, the problem is below
# orientation entirely -- the trace cannot resolve even plain collinear homology
# on a graph this small -- and the inversion tests are not evidence of anything.
GFA_COLLINEAR = """H\tVN:Z:1.1
S\t1\tACGTACGTTTGGCCAA
S\t2\tGGATTCCAGTAC
S\t3\tTTACCGGATCAGGTTA
S\t4\tGATTACAGATTAC
W\tS1\t1\tchr1\t0\t57\t>1>2>3>4
W\tS2\t2\tchr1\t0\t57\t>1>2>3>4
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
"""

COORD_FILES = ("rlbwt_rindex.ri", "sampled.tags", "fastlocate.ri",
               "output.t1", "output.t2")

SRC_CONTIG = "S1#1#chr1"
TGT_HAP = "S2#2"

# Ground truth. On an inverted block the source and target offsets sum to a
# constant: target = (target_path_length - 1) - reverse_sequence_offset, and the
# two shared nodes are adjacent in both haplotypes. Derived, not guessed:
#   fixture A: S2 is 52bp, shared core is S1[16,44) -> src + tgt == 55
#   fixture B: S2 is 57bp, shared core is S1[16,44) -> src + tgt == 59
CORE_SRC_START, CORE_SRC_END = 16, 44
CONTAINED_SUM = 55
BOUNDARY_SUM = 59
BOUNDARY_LEN = 57


def _nonempty(p: Path) -> bool:
    """Exists AND has content -- the fixture builder tests with `[[ -s ]]`, so a
    zero-byte leftover must not count as present."""
    return p.is_file() and p.stat().st_size > 0


def _tools_available() -> bool:
    vg = os.environ.get("VG") or shutil.which("vg")
    return bool(vg) and all(shutil.which(t) for t in ("gbz_extract", "grlbwt-cli"))


def _build_index(liftover, fixture: Path, gfa_text: str):
    """Build the coordinate indexes for one embedded graph, then load them."""
    if not all(_nonempty(fixture / f) for f in COORD_FILES):
        if not _tools_available():
            pytest.skip(
                "inversion fixture indexes are missing and the build toolchain "
                "(vg, gbz_extract, grlbwt-cli) is not on PATH"
            )
        vg = os.environ.get("VG") or shutil.which("vg")

        # The builder stages its input with `cp -f "$GBZ_IN" "$OUT/graph.gbz"`,
        # so the source GBZ must live OUTSIDE the fixture directory -- handing it
        # the fixture's own graph.gbz makes cp fail with "are the same file".
        with tempfile.TemporaryDirectory() as td:
            src = Path(td)
            gfa = src / "graph.gfa"
            gfa.write_text(gfa_text)
            src_gbz = src / "graph.gbz"
            proc = subprocess.run(
                [vg, "gbwt", "--gbz-format", "-g", str(src_gbz), "-G", str(gfa)],
                capture_output=True, text=True,
            )
            if proc.returncode != 0 or not _nonempty(src_gbz):
                pytest.skip(
                    f"could not build the fixture GBZ (rc={proc.returncode})\n"
                    f"stdout: {proc.stdout[-600:]}\nstderr: {proc.stderr[-600:]}"
                )

            fixture.mkdir(parents=True, exist_ok=True)
            # Keep the graph next to the indexes for human inspection.
            (fixture / "graph.gfa").write_text(gfa_text)

            env = dict(os.environ)
            env.setdefault("BIN", str(REPO_ROOT / "bin"))
            env.setdefault("VG", vg or "")
            env.setdefault("THREADS", "2")
            proc = subprocess.run(
                ["bash", str(BUILDER), str(src_gbz), str(fixture), "--coord-only"],
                env=env, capture_output=True, text=True,
            )

        missing = [f for f in COORD_FILES if not _nonempty(fixture / f)]
        if proc.returncode != 0 or missing:
            pytest.skip(
                f"could not build fixture {fixture.name} (rc={proc.returncode}, "
                f"missing={missing})\nstdout: {proc.stdout[-800:]}\n"
                f"stderr: {proc.stderr[-1200:]}"
            )

    idx = liftover.Index()
    idx.load(
        str(fixture / "graph.gbz"),
        str(fixture / "rlbwt_rindex.ri"),
        str(fixture / "sampled.tags"),
        str(fixture / "fastlocate.ri"),
        str(fixture / "output.t1"),
        str(fixture / "output.t2"),
    )
    return idx



@pytest.fixture(scope="module")
def collinear_index(liftover):
    return _build_index(
        liftover,
        REPO_ROOT / "coord_translation_tests" / "test4_collinear_control",
        GFA_COLLINEAR)


@pytest.fixture(scope="module")
def contained_index(liftover):
    return _build_index(
        liftover, REPO_ROOT / "coord_translation_tests" / "test2_inversion",
        GFA_CONTAINED)


@pytest.fixture(scope="module")
def boundary_index(liftover):
    return _build_index(
        liftover,
        REPO_ROOT / "coord_translation_tests" / "test3_inversion_boundary",
        GFA_BOUNDARY)


def _pairs(recs, strand=None):
    """{source_offset: target_offset} for the records on one strand."""
    return {r.start: r.end for r in recs if strand is None or r.strand == strand}


def _describe(recs):
    return sorted((r.start, r.end, r.strand) for r in recs)



# ── Control: plain collinear homology, no rearrangement ─────────────────────

def test_control_collinear_translates_at_all(collinear_index):
    """Baseline. Two haplotypes over the identical node walk MUST translate.

    If this fails, nothing below it is interpretable: the trace cannot resolve
    collinear homology on a graph this size, so an inversion test failing says
    nothing about orientation handling.
    """
    res = collinear_index.translate(SRC_CONTIG, 0, BOUNDARY_LEN, TGT_HAP)
    assert res, (
        "two haplotypes over the SAME node walk produced no translation at all; "
        "this is a failure below orientation handling"
    )


def test_control_collinear_is_the_identity(collinear_index):
    """Identical walks means every mapped base satisfies src == tgt, '+'."""
    res = collinear_index.translate(SRC_CONTIG, 0, BOUNDARY_LEN, TGT_HAP)
    assert res, "no intervals (see test_control_collinear_translates_at_all)"
    assert all(r.strand == "+" for r in res), (
        f"collinear graph reported a reverse strand: {_describe(res)[:10]}"
    )
    wrong = {s: t for s, t in _pairs(res).items() if t != s}
    assert not wrong, f"expected src == tgt everywhere, got {sorted(wrong.items())[:10]}"
    assert len(_pairs(res)) >= 40, (
        f"only {len(_pairs(res))} of 57 bases mapped: {_describe(res)[:10]}"
    )


# ── Fixture A: inversion fully contained in the query ────────────────────────

def test_contained_fixture_has_both_haplotypes(contained_index):
    """Guards the fixture itself: if the graph or its metadata changes, the
    orientation tests below would skip or pass for the wrong reason."""
    names = contained_index.get_haplotype_names()
    assert any(n.startswith("S1#1") for n in names), names
    assert any(n.startswith("S2#2") for n in names), names


def test_contained_inversion_is_translated_at_all(contained_index):
    """The inverted block must translate to SOMETHING on the other haplotype."""
    res = contained_index.translate(SRC_CONTIG, CORE_SRC_START, CORE_SRC_END,
                                    TGT_HAP)
    assert res, (
        f"{SRC_CONTIG}[{CORE_SRC_START},{CORE_SRC_END}) -> {TGT_HAP} produced no "
        "intervals; the block is shared, only its orientation differs"
    )


def test_contained_inversion_maps_to_the_right_bases(contained_index):
    """Every mapped base must satisfy src + tgt == 55, the exact relation for
    this fixture -- not merely land near the right span."""
    res = contained_index.translate(SRC_CONTIG, CORE_SRC_START, CORE_SRC_END,
                                    TGT_HAP)
    assert res, "no intervals (see test_contained_inversion_is_translated_at_all)"
    pairs = _pairs(res)
    wrong = {s: t for s, t in pairs.items() if s + t != CONTAINED_SUM}
    assert not wrong, (
        f"expected src + tgt == {CONTAINED_SUM} for every mapped base, "
        f"got {sorted(wrong.items())[:10]} (all: {_describe(res)[:10]})"
    )
    assert len(pairs) >= 20, f"only {len(pairs)} bases mapped: {_describe(res)}"


def test_contained_inversion_reports_minus_strand(contained_index):
    """The whole point: an inverted correspondence must be reported as '-'."""
    res = contained_index.translate(SRC_CONTIG, CORE_SRC_START, CORE_SRC_END,
                                    TGT_HAP)
    assert res, "no intervals (see test_contained_inversion_is_translated_at_all)"
    assert all(r.strand == "-" for r in res), (
        f"inverted block reported as forward: {_describe(res)}"
    )


def test_inverted_haplotype_is_still_discoverable(contained_index):
    """Routing tries both orientations, so the inverted haplotype must be listed
    as reachable.

    This is the asymmetry worth pinning: a haplotype advertised as reachable but
    yielding no coordinates is worse than one that is absent. It also fails if
    routing ever regresses to forward-only.
    """
    names = contained_index.get_haplotype_names()
    assert next((n for n in names if n.startswith("S2#2")), None), names
    try:
        reachable = contained_index.translatable_haplotypes(
            SRC_CONTIG, CORE_SRC_START, CORE_SRC_END)
    except Exception as exc:  # pragma: no cover - extension without the call
        pytest.skip(f"translatable_haplotypes unavailable: {exc}")
    assert any(h.startswith("S2#2") for h in reachable), (
        f"inverted haplotype not discoverable at all: {reachable}"
    )


# ── Fixture B: query spans the inversion boundary ────────────────────────────

def test_boundary_fixture_has_both_haplotypes(boundary_index):
    names = boundary_index.get_haplotype_names()
    assert any(n.startswith("S1#1") for n in names), names
    assert any(n.startswith("S2#2") for n in names), names


def test_boundary_query_reports_both_strands(boundary_index):
    """The regression this fixture exists for.

    A fallback gated on "the forward trace found nothing" cannot fire here: the
    collinear flanks make the forward trace succeed, so the inverted core is
    dropped and the answer looks like a clean forward translation. Both strands
    must be present.
    """
    res = boundary_index.translate(SRC_CONTIG, 0, BOUNDARY_LEN, TGT_HAP)
    assert res, "whole-contig query produced no intervals at all"
    strands = {r.strand for r in res}
    assert strands == {"+", "-"}, (
        f"expected collinear flanks ('+') AND an inverted core ('-'), got "
        f"strands={sorted(strands)}: {_describe(res)}"
    )


def test_boundary_flanks_map_collinearly(boundary_index):
    """The flanks are identical in both haplotypes, so src == tgt there."""
    res = boundary_index.translate(SRC_CONTIG, 0, BOUNDARY_LEN, TGT_HAP)
    assert res, "no intervals (see test_boundary_query_reports_both_strands)"
    fwd = _pairs(res, "+")
    assert fwd, f"no forward-strand bases at all: {_describe(res)}"
    wrong = {s: t for s, t in fwd.items() if t != s}
    assert not wrong, (
        f"flanks are identical in both haplotypes so src == tgt; "
        f"got {sorted(wrong.items())[:10]}"
    )
    # Both flanks, not just one. (Source offset 0 is discarded by the tracer's
    # target_offset==0 "unmapped" sentinel, hence the >= 1 lower bound.)
    assert any(1 <= s < CORE_SRC_START for s in fwd), (
        f"left flank missing: {sorted(fwd)[:10]}")
    assert any(CORE_SRC_END <= s < BOUNDARY_LEN for s in fwd), (
        f"right flank missing: {sorted(fwd)[-10:]}")


def test_boundary_core_is_inverted_and_exact(boundary_index):
    """The core must be reported as '-' and map to exactly the right bases."""
    res = boundary_index.translate(SRC_CONTIG, 0, BOUNDARY_LEN, TGT_HAP)
    assert res, "no intervals (see test_boundary_query_reports_both_strands)"
    rev = _pairs(res, "-")
    assert rev, (
        f"the inverted core was dropped entirely -- this is the silent-hole "
        f"failure: {_describe(res)}"
    )
    off_core = {s: t for s, t in rev.items()
                if not (CORE_SRC_START <= s < CORE_SRC_END)}
    assert not off_core, (
        f"'-' strand reported outside the inverted core [16,44): "
        f"{sorted(off_core.items())[:10]}"
    )
    wrong = {s: t for s, t in rev.items() if s + t != BOUNDARY_SUM}
    assert not wrong, (
        f"expected src + tgt == {BOUNDARY_SUM} across the inverted core, "
        f"got {sorted(wrong.items())[:10]}"
    )
    assert len(rev) >= 20, (
        f"only {len(rev)} of 28 core bases mapped: {sorted(rev.items())}"
    )


def test_boundary_core_not_silently_claimed_as_forward(boundary_index):
    """The failure mode that is worse than returning nothing: reporting the
    inverted core as a forward match."""
    res = boundary_index.translate(SRC_CONTIG, 0, BOUNDARY_LEN, TGT_HAP)
    assert res, "no intervals (see test_boundary_query_reports_both_strands)"
    bogus = {s: t for s, t in _pairs(res, "+").items()
             if CORE_SRC_START <= s < CORE_SRC_END}
    assert not bogus, (
        f"inverted core reported on the forward strand: "
        f"{sorted(bogus.items())[:10]}"
    )
