# Lifting a region to every haplotype — browser integration

`POST /api/v1/liftover` can translate a source region to **all 464 haplotypes in
one call**. The server makes a single pass over the source region for all of them,
so asking for every haplotype costs about as much as asking for a few.

Authentication is the same as the other endpoints: `X-Pangenome-Token: <token>` or
`Authorization: Bearer <token>`.

---

## Request

```json
{
  "src": "CHM13#0#chr17",
  "start": 43000000,
  "end": 43100000,
  "tgt": "*",
  "summary": true
}
```

| field | type | default | meaning |
|---|---|---|---|
| `src` | string | — | full source contig path, `sample#hap#contig` |
| `start`, `end` | int | — | 0-based, half-open `[start, end)`; at most 10,000,000 bp |
| `tgt` | string or array | — | `"*"` = every haplotype; or one name / a list (max 1024). A name may be a haplotype (`HG00097#1`) or a contig (`GRCh38#0#chr17`). |
| `summary` | bool | `false` | one compact row per target (see below). **Use this for `"*"`.** |
| `blocks` | bool | `false` | full colinear block alignment instead of intervals. Cannot be combined with `summary` (400). |
| `timeout_ms` | number | `5000` | time budget **per target**, max 60000. For many targets the budget is shared: `timeout_ms × number of targets`. |

---

## Response: `summary: true`

One row per requested target, in request order, **including targets that the
region does not reach** (those have `pieces: []`, `mapped_bp: 0`). For `"*"`, the
source's own haplotype is in the list with `mapped_bp: 0`.

```json
{
  "summary": [
    {
      "haplotype": "GRCh38#0",
      "pieces": [
        {"contig": "chr17", "start": 42143496, "end": 42243512, "strand": "+",
         "source_start": 43000000, "source_end": 43100000}
      ],
      "mapped_bp": 99960,
      "source_covered": [43000000, 43100000]
    }
  ],
  "warnings": []
}
```

(real output)

| field | meaning |
|---|---|
| `haplotype` | the target name as requested |
| `pieces` | where the region lands: `contig` is the short contig name, `start`/`end` are target coordinates, `source_start`/`source_end` are the part of the source that piece covers. More than one piece means the region is split (contig break, inversion, duplication). |
| `mapped_bp` | number of source bases that have a position on this haplotype. Use `mapped_bp / (end - start)` as the "how much of the region is there" fraction. |
| `source_covered` | `[first, last)` source base that maps, or `null` if nothing maps |

Pieces from the same haplotype may sit on different contigs. A haplotype can also
have pieces on both strands. When part of the region is duplicated on that
haplotype, each copy is its own piece, so pieces can overlap in source.

---

## Response: without `summary`

The old forms still work, now with many targets in one call:

- default → `"intervals": [{haplotype, start, end, strand}, ...]`
- `"blocks": true` → `"blocks": [{haplotype, source_start, source_end, target_start, target_end, strand}, ...]`

Here `haplotype` is the **full target contig path** (`GRCh38#0#chr17`), not the
two-field name. Rows are sorted by that path, then by position. For `"*"` these
lists get large (about 118,000 blocks for 1 Mb), so prefer `summary` for an
overview and request `blocks` for the one or two haplotypes the user opens.

---

## Warnings and errors

- A timeout does not fail the request. It returns 200 with what was done and a
  warning `{haplotype, reason: "timeout", elapsed_ms, message}`. In summary mode
  there is one warning with `haplotype: "*"`, meaning the results may be
  incomplete.
- 400: bad interval, span over 10 Mb, unknown `src`, bad `tgt`, or `blocks` and
  `summary` both set.
- 503: server still loading its indexes. Retry later.

---

## Timing

All 464 haplotypes, measured on the v2 graph:

| region size | time |
|---|---|
| 1 kb | ~10 ms |
| 10 kb | ~50 ms |
| 100 kb | ~0.4 s |
| 1 Mb | ~4.5 s |
| 10 Mb | ~50 s |

Concurrent requests run in parallel (4 × 100 kb all-haplotype summaries took
0.48 s together). The time grows with region size, so set your HTTP client
timeout to at least ~6 s per Mb, or limit the size of region you offer for
`"*"`. A 100 kb summary for all haplotypes is about 107 kB of JSON.

---

## Suggested flow

```
1. POST /api/v1/liftover {src, start, end, tgt: "*", summary: true}
   -> table of all haplotypes: mapped fraction, contig(s), strand
2. user opens a haplotype
   -> POST /api/v1/liftover {src, start, end, tgt: "<haplotype>", blocks: true}
      for the detailed alignment (one target, 1 Mb ≈ 0.1 s)
```
