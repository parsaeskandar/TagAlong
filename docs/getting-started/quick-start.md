# Quick Start Pipeline

This is the shortest path to build core assets and run query tools.

## 1) Prepare graph and RL-BWT

```bash
gbz_extract -t 8 -b your_graph.gbz > graph_info
grlbwt-cli -t 8 graph_info
# produces graph_info.rl_bwt
```

## 2) Build tags and r-index

```bash
./bin/build_tags your_graph.gbz graph_info.rl_bwt output.tags
./bin/build_rindex graph_info.rl_bwt > output.ri
```

## 3) Convert tags (required by query/find_mems)

```bash
./bin/convert_tags output.tags output_compressed.tags --gbz graph.gbz
```

## 4) Run downstream tools

- MEM search:
```bash
./bin/find_mems output.ri output_compressed.tags reads.txt 10 1
```

- Tag query:
```bash
./bin/query_tags output.ri output_compressed.tags reads.txt
```

For coordinate translation, follow the full prerequisite flow:
- [coordinate_translation Guide](../tools/coordinate_translation.md)

