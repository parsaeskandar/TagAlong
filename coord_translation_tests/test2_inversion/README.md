# test2_inversion — coordinate translation across an inverted block

Two haplotypes over one contig. The shared block (nodes 2 and 3) is traversed
in the OPPOSITE order by the two, i.e. an inversion:

    S1#1#chr1   >1>2>3>4     nodes 1,2,3,4 forward
    S2#2#chr1   >5<3<2>6     nodes 3,2 reversed, flanked by private 5 and 6

Node lengths: 1=16  2=12  3=16  4=13  5=12  6=12

    S1:  1 [0,16)   2 [16,28)   3 [28,44)   4 [44,57)      total 57
    S2:  5 [0,12)   3 [12,28)   2 [28,40)   6 [40,52)      total 52

Shared nodes are 2 and 3. On S1 they run 2 then 3 (forward); on S2 they run
3 then 2 (reversed). So the ground truth is:

    translate(S1#1#chr1, 16, 44, S2#2)  ->  [12, 40)  strand '-'

and symmetrically

    translate(S2#2#chr1, 12, 40, S1#1)  ->  [16, 44)  strand '-'

Nodes 1 and 4 are private to S1, nodes 5 and 6 private to S2, so nothing
outside the inverted block can be confused for it.

## Building the coordinate indexes

`graph.gfa` and `graph.gbz` are committed; the coordinate indexes are derived
and are built on demand by the test (see tests/python/test_inversion.py), or
by hand:

    BIN=$PWD/bin VG=$(command -v vg) \
      tests/fixtures/build_fixture_indexes.sh \
        coord_translation_tests/test2_inversion/graph.gbz \
        coord_translation_tests/test2_inversion --coord-only

That needs `gbz_extract` and `grlbwt-cli` on PATH.
