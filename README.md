# Overview

Modified TREX pipeline, with only `trex run10x` kept.

## Example of run

```
trex run10x \
  --chromosome <name_of_added_sequence> \
  --start <start_position> \
  --end <end_position> \
  --delete \
  --output <output_folder> \
  --amplicon <path_to_Cell_Ranger_amp_alignment> \
  <path_to_Cell_Ranger_alignment>
```

Additionally to the original pipeline, .h5ad-file with the result of barcode calling per each cell will be generated.
