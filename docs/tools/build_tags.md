# build_tags

Constructs tag arrays from a GBZ graph and RL-BWT.

## Usage

```bash
./bin/build_tags <graph.gbz> <graph_info.rl_bwt> <output.tags>
```

## Inputs
- `graph.gbz`
- `graph_info.rl_bwt` (from `grlbwt-cli`)

## Output
- `output.tags` (algorithm-format tags)

## Next step
If using query/find_mems, convert tags:

```bash
./bin/convert_tags output.tags output_compressed.tags --gbz graph.gbz
```

