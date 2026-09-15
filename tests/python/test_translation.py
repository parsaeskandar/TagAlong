"""Coordinate-translation tests (liftover_ext.Index).

Two tiers:

* **contract** tests run on any loadable fixture (including the degenerate
  in-repo default). They pin the API shape, input validation, and determinism —
  the things that must hold regardless of graph content.
* **oracle** tests (``@pytest.mark.oracle``) need a real multi-haplotype fixture
  with actual homology. They discover a translatable interval at runtime and
  self-skip if the fixture yields none, so they are harmless on the degenerate
  default and meaningful on a chrM/chrY fixture.

Run:  pytest tests/python -q
Oracle:  PANGENOME_TEST_FIXTURE=/path/to/chrM_fixture pytest tests/python -q
"""
from __future__ import annotations

import pytest


def _tuples(records):
    """Normalize TranslatedInterval records to comparable tuples."""
    return [(r.haplotype, r.start, r.end, r.strand) for r in records]


# IMPORTANT, and the source of three separate test bugs so far: despite the name,
# a ``TranslatedInterval`` from translate() is NOT an interval. It is a per-base
# correspondence -- ``start`` is an offset on the SOURCE haplotype and ``end`` is
# the offset it maps to on the TARGET. One record per base, not one per region.
#
# So ``(r.start, r.end)`` is a (source, target) PAIR. Feeding it back into
# translate() as if it were [start, end) mixes two coordinate systems and cannot
# work; comparing ``r.end`` against a source interval's end is likewise
# meaningless. Both mistakes are easy to miss because on a collinear region the
# numbers look plausible.
def _pairs(records, haplotype=None):
    """{source_offset: target_offset} for records, optionally on one haplotype."""
    return {r.start: r.end for r in records
            if haplotype is None or r.haplotype == haplotype}


# --------------------------------------------------------------------------- #
# Contract: error / validation paths
# --------------------------------------------------------------------------- #

def test_translate_before_load_raises(liftover):
    """Using the index before load() is a usage error, not a crash."""
    idx = liftover.Index()
    with pytest.raises(RuntimeError):
        idx.translate("_gbwt_ref#4294967295", 0, 10, "_gbwt_ref#4294967295")


def test_get_haplotype_names_is_nonempty_str_list(haplotype_names):
    assert isinstance(haplotype_names, list)
    assert haplotype_names
    assert all(isinstance(n, str) and n for n in haplotype_names)


def test_get_haplotype_names_is_deterministic(coord_index):
    assert coord_index.get_haplotype_names() == coord_index.get_haplotype_names()


def test_translate_returns_list_of_typed_intervals(coord_index, haplotype_names):
    src = haplotype_names[0]
    result = coord_index.translate(src, 0, 50, src)
    assert isinstance(result, list)
    for rec in result:
        assert isinstance(rec.haplotype, str) and rec.haplotype
        assert isinstance(rec.start, int)
        assert isinstance(rec.end, int)
        assert rec.start <= rec.end
        # NOTE: this cannot currently fail -- translate() hardcodes '+' at both
        # call sites (pangenome_server.cpp). It is kept as a type guard, and
        # test_inversion.py carries the test that actually exercises '-'.
        assert rec.strand in ("+", "-")


@pytest.mark.parametrize("start,end", [(100, 99), (50, 0)])
def test_translate_rejects_inverted_interval(coord_index, haplotype_names, start, end):
    with pytest.raises(ValueError):
        coord_index.translate(haplotype_names[0], start, end, haplotype_names[0])


def test_translate_rejects_overlong_interval(coord_index, haplotype_names):
    # The binding rejects intervals longer than MAX_INTERVAL_LENGTH (1e7).
    with pytest.raises(ValueError):
        coord_index.translate(haplotype_names[0], 0, 10_000_001, haplotype_names[0])


def test_translate_rejects_missing_source(coord_index, haplotype_names):
    with pytest.raises(ValueError):
        coord_index.translate("__no_such_haplotype__", 0, 10, haplotype_names[0])


def test_translate_is_deterministic(coord_index, haplotype_names):
    src = haplotype_names[0]
    first = coord_index.translate(src, 0, 100, src)
    second = coord_index.translate(src, 0, 100, src)
    assert _tuples(first) == _tuples(second)


# --------------------------------------------------------------------------- #
# Oracle: real homology required (self-skips on the degenerate fixture)
# --------------------------------------------------------------------------- #

def _first_translatable(idx, src, tgt, span=200, step=500, max_offset=200_000):
    """Scan offsets on ``src`` for the first [off, off+span) that yields >=1
    translation to ``tgt``. Returns (start, end, records) or None."""
    off = 0
    while off <= max_offset:
        try:
            res = idx.translate(src, off, off + span, tgt)
        except ValueError:
            res = []
        if res:
            return off, off + span, res
        off += step
    return None


@pytest.mark.oracle
def test_identity_translation_preserves_interval(coord_index, haplotype_names):
    """translate(A, i, j, A) must return the same interval on A, forward strand.

    Self-oracle: the source haplotype's own sequence is ground truth, so an
    identity lift is exact by construction.
    """
    src = haplotype_names[0]
    found = _first_translatable(coord_index, src, src)
    if not found:
        # Not a degenerate fixture -- Table 2 deliberately stores no route from a
        # haplotype to itself (build_translation_tables.cpp: `if (!allow_same_hap
        # && h == path_hap[sp]) return;`, and --allow-same-haplotype defaults
        # off). So the T2 serving path cannot answer an identity query at all.
        # The table-free path can, since it discovers candidates by probing.
        pytest.skip(
            f"no self-translatable interval on {src!r}: Table 2 excludes "
            "same-haplotype pairs by design (build with --allow-same-haplotype, "
            "or run this under PANGENOME_TRANSLATE_NO_T2=1)"
        )
    start, end, res = found
    same = [r for r in res if r.haplotype == src]
    assert same, f"identity produced no record on {src!r}: {_tuples(res)}"

    # Records are (source, target) pairs, so an identity lift means every base
    # maps to itself: src == tgt, forward strand.
    pairs = _pairs(same)
    wrong = {a: b for a, b in pairs.items() if a != b}
    assert not wrong, (
        f"identity did not map bases to themselves on {src!r}: "
        f"{sorted(wrong.items())[:10]}"
    )
    assert all(r.strand == "+" for r in same), (
        f"identity reported a reverse strand on {src!r}: "
        f"{[(r.start, r.end, r.strand) for r in same if r.strand != '+'][:10]}"
    )
    assert pairs, f"identity produced no usable pairs on {src!r}"


@pytest.mark.oracle
def test_translation_round_trips(coord_index, haplotype_names):
    """translate(A->B) then translate(B->A) must recover the original interval
    (exact for SNP-only regions, within an indel-sized tolerance elsewhere)."""
    if len(haplotype_names) < 2:
        pytest.skip("round-trip needs >=2 haplotypes")

    src = haplotype_names[0]
    found = None
    for tgt in haplotype_names[1:]:
        found = _first_translatable(coord_index, src, tgt)
        if found:
            break
    if not found:
        pytest.skip(f"no cross-haplotype homology found from {src!r}")

    start, end, forward = found
    span = end - start
    tolerance = 50  # allow for indels between the two haplotypes
    recovered = False
    strand_mismatch = None
    # Each `fwd` is a (source_offset -> target_offset) pair. To go back, ask for
    # an interval on the TARGET that begins at the target offset, then look for a
    # returned pair whose own target offset lands back on the source offset we
    # started from. Passing (fwd.start, fwd.end) as an interval -- as this test
    # used to -- spans two different haplotypes' coordinate systems.
    for fwd in forward:
        src_off, tgt_off = fwd.start, fwd.end
        back = coord_index.translate(fwd.haplotype, tgt_off, tgt_off + span, src)
        hits = [
            b for b in back
            if b.haplotype == src
            and abs(b.start - tgt_off) <= tolerance
            and abs(b.end - src_off) <= tolerance
        ]
        if hits:
            # Strand must be symmetric: if A->B is inverted then B->A is too.
            # Checked only on recovered records, so a fixture with no inverted
            # homology simply sees '+' both ways and passes.
            if not any(b.strand == fwd.strand for b in hits):
                strand_mismatch = (fwd.strand, [b.strand for b in hits])
            recovered = True
            break
    assert recovered, (
        f"round-trip {src} -> {forward[0].haplotype} -> {src} did not recover any "
        f"of the {len(forward)} base correspondences from [{start},{end}) within "
        f"{tolerance}bp; first pair was source {forward[0].start} -> target "
        f"{forward[0].end}"
    )
    assert strand_mismatch is None, (
        f"round-trip strand is not symmetric: forward reported "
        f"{strand_mismatch[0]!r} but the return leg reported "
        f"{strand_mismatch[1]!r}"
    )
