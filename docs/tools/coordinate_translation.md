# coordinate_translation

Translate coordinates from a source haplotype interval to a target haplotype using:

- r-index (`.ri`)
- GBWT FastLocate (`.ri`)
- sampled tags (`.tags`)
- GBZ graph (`.gbz`)
- translation tables (Table 1 + Table 2)

---

## What you need before running

You need all required artifacts in place. Build/order is:

1. **Graph and RL-BWT**
  - `graph.gbz`
  - `graph_info.rl_bwt`
  - See: [Quick Start Pipeline](../getting-started/quick-start.md)
2. **R-index**
  - `whole_compact2.ri` (or equivalent)
  - Built with `build_rindex`
  - See: [build_rindex](build_rindex.md)
3. **Sampled tags**
  - `sampled.tags` (the sampled tag array used by coordinate translation)
  - If your upstream pipeline creates whole-genome/compressed tags first, make sure the sampled tags used here are from the same graph/index build.
  - Upstream references:
    - [build_tags](build_tags.md)
    - [convert_tags](convert_tags.md)
    - [Whole-Genome Pipeline](../pipelines/whole-genome-pipeline.md)
4. **GBWT FastLocate index**
  - `fastlocate.ri` (GBWT-side fast locate index used by translation – build using [`vg`](https://github.com/vgteam/vg) from vgteam)
5. **Translation tables**
  - Table 1 file (`output.t1`)
  - Table 2 file (`output.t2`)
  - **How to create:**  
    Both tables are built together by one positional command:
    ```bash
    ./bin/build_translation_tables \
      <graph.gbz> \
      <fastlocate.ri> \
      <output.t1> \
      <output.t2> \
      [--threads N] [--only-table1]
    ```
    The GBWT FastLocate r-index is required: Table 2 records which target paths a
    source range can reach, and resolves them with `decompressSA`. The sampled
    tags are NOT an input here — they are consumed at query time, not at build
    time. Build every artifact from the *same* graph. Prefer the non-filtered
    graph: Table 2 scales with GBWT path count, and frequency filtering multiplies
    that by ~2,400 without adding information. For whole-genome workflows see the
    [Whole-Genome Pipeline](../pipelines/whole-genome-pipeline.md).

---

## Single translation command

```bash
/usr/bin/time -v ./bin/coordinate_translation \
  <rlbwt_rindex.ri> \
  <gbwt_fastlocate.ri> \
  <sampled.tags> \
  --gbz <graph.gbz> \
  --interval <START>..<END> \
  --table1 <output.t1> \
  --table2 <output.t2> \
  --source-haplotype-name <SRC_HAP> \
  --source-contig <SRC_CONTIG> \
  --target-haplotype-name <TGT_HAP> \
  --target-contig <TGT_CONTIG>
```

Example:

```bash
/usr/bin/time -v ./bin/coordinate_translation \
  ../whole_genome_rindex/whole_compact2.ri \
  fastlocate.ri \
  ../merge_tags_Feb-11/sampled_FEB-2026_correct.tags \
  --gbz /path/to/hprc-v2.0-mc-chm13.gbz \
  --interval 38368638..43368638 \
  --table1 ../translation_tables/output.t1 \
  --table2 ../translation_tables/output.t2 \
  --source-haplotype-name CHM13#0 \
  --source-contig chr10 \
  --target-haplotype-name GRCh38#0 \
  --target-contig chr10
```

---

## Table-driven algorithm used by this tool

1. **Resolve source interval with Table 1**
  - Source haplotype + interval are resolved to one or more `(src_path_id, local_start, local_end)` sub-intervals.
2. **Get candidate targets with Table 2**
  - For each source sub-interval, lookup `(src_path_id, target_haplotype, local interval)` to get candidate target path IDs.
  - Apply optional target contig filter.
  - Deduplicate target path IDs.
3. **Restrict translation range using Table 2 extent**
  - For each distinct target path, compute min/max overlap from Table 2 segments.
  - Run `find_tags_in_interval` and tracing only inside that extent.
4. **Run translation**
  - `find_tags_in_interval`
  - `find_first_and_last_common_nodes_gbwt`
  - `trace_coordinates_gbwt`
5. **Convert to haplotype coordinates**
  - Path-local offsets are converted back using Table 1 subpath starts.
  - Output is sorted by source coordinate in single-run mode.

---

## Benchmark mode (table-driven)

```bash
./bin/coordinate_translation \
  <rlbwt_rindex.ri> <gbwt_fastlocate.ri> <sampled.tags> \
  --gbz <graph.gbz> \
  --benchmark \
  --benchmark-intervals 5000000 \
  --table1 <output.t1> \
  --table2 <output.t2> \
  --source-haplotype-name <SRC_HAP> \
  --source-contig <SRC_CONTIG> \
  --target-haplotype-name <TGT_HAP> \
  --target-contig <TGT_CONTIG>
```

This mode prints:

- how many source intervals/path IDs came from Table 1
- how many target paths came from Table 2 (before/after filtering)
- per-run timing breakdown

