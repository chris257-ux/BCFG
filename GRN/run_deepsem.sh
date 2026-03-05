INPUT_CSV="/home/chenliqun/FinalGene5000_tf1800/AnaplasticThyroidCancer_malignant_aneuploid_cells_expression_transposed_5000TF_transpose.csv"
OUTPUT_DIR="/home/chenliqun/FinalGene5000_tf1800/deepsemResults/AnaplasticThyroidCancer_result"
DEEPSEM_PATH="/home/chenliqun/deepSEM/DeepSEM"

# DeepSEM parameters
TASK="celltype_GRN"
SETTING="test"
ALPHA=0.1
BETA=0.01
EPOCHS=150

# Create output directory
mkdir -p "${OUTPUT_DIR}" || { echo "Failed to create output directory"; exit 1; }

# Enter DeepSEM directory
cd "${DEEPSEM_PATH}" || { echo "Failed to change to DeepSEM directory"; exit 1; }

# Run command
echo "Processing: $(basename ${INPUT_CSV})"
echo "Output directory: ${OUTPUT_DIR}"

python main.py \
--task "${TASK}" \
--data_file "${INPUT_CSV}" \
--setting "${SETTING}" \
--alpha "${ALPHA}" \
--beta "${BETA}" \
--n_epochs "${EPOCHS}" \
--save_name "${OUTPUT_DIR}/GRN_output" 2>&1 | tee "${OUTPUT_DIR}/run.log"

# Remove trailing slash if present
OUTPUT_DIR=$(echo "${OUTPUT_DIR}" | sed 's:/*$::')

if [ $? -eq 0 ] && [ -f "${OUTPUT_DIR}/GRN_output/GRN_inference_result.tsv" ]; then
echo -e "\n=== Processing succeeded ==="
echo "Generated result files:"
ls -lh "${OUTPUT_DIR}/GRN_output"*
else
echo -e "\n=== Processing failed ==="
echo "Last 10 lines of error log:"
tail -n 10 "${OUTPUT_DIR}/run.log"
exit 1
fi