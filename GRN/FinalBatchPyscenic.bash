#!/bin/bash
# Run in background example:
# nohup bash batchPyscenic.bash > output.txt 2> error.txt &

output_dir="/home/chenliqun/FinalGene5000_tf1800"
loom_dir="$output_dir/loomfiles"  # Directory containing .loom files
result_dir="$output_dir/pyscenicResults"  # Directory for output results

# Ensure output directory exists
mkdir -p "$result_dir"

# Database file paths
dir="/home/chenliqun/pyscenic/cisTarget_databases"
tfs="$dir/hs_hgnc_tfs.txt"
feather="$dir/hg19-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather"
tbl="$dir/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

# Iterate through all .loom files
for input_loom in "$loom_dir"/*.loom; do
  # Get file base name (without path and extension)
  base_name=$(basename "$input_loom" .loom)

  # Define output file paths
  adj_output="$result_dir/adj.${base_name}.tsv"
  reg_output="$result_dir/reg.${base_name}.csv"
  auc_output="$result_dir/out_SCENIC.${base_name}.loom"

  # Check required files
  echo "Processing $input_loom"
  ls "$tfs" "$feather" "$tbl" "$input_loom"

  # 2.1 GRN inference
  echo "Running GRN inference for $input_loom"
  pyscenic grn \
    --num_workers 20 \
    --output "$adj_output" \
    --method grnboost2 \
    "$input_loom" \
    "$tfs"

  # 2.2 cisTarget pruning and validation
  echo "Running cisTarget pruning for $input_loom"
  pyscenic ctx \
    "$adj_output" \
    "$feather" \
    --annotations_fname "$tbl" \
    --expression_mtx_fname "$input_loom" \
    --mode "dask_multiprocessing" \
    --output "$reg_output" \
    --num_workers 20 \
    --mask_dropouts

  # 2.3 AUCell calculation of regulon activity
  echo "Running AUCell for $input_loom"
  pyscenic aucell \
    "$input_loom" \
    "$reg_output" \
    --output "$auc_output" \
    --num_workers 20
done

echo "All files processed!"