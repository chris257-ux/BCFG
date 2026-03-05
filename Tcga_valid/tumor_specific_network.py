import pandas as pd
import os

BASE_DIR="/home/chenliqun/newFinalResults"
WEIGHT_TAG="SD053047"

# Output directory
OUTPUT_DIR=os.path.join(BASE_DIR,"tumor_specific")
os.makedirs(OUTPUT_DIR,exist_ok=True)

tumor_normal_pairs={
"BladderUrothelialCarcinoma":"Bladder",
"CervicalSquamousCellCarcinoma":"Uterus",
"ChromophobeRenalCellCarcinoma":"Kidney",
"ClearCellRenalCellCarcinoma":"Kidney",
"ColorectalCancer":"Colon",
"CutaneousSquamousCellCarcinoma":"Skin",
"EndometrialCarcinoma":"Uterus",
"ERPositiveBreastCancer":"Breast"
}

all_results=[]

for tumor_type,normal_tissue in tumor_normal_pairs.items():
    print(f"Processing {tumor_type} vs {normal_tissue}")

    tumor_file=os.path.join(BASE_DIR,f"cancer_{tumor_type}_consensus_{WEIGHT_TAG}.tsv")
    normal_file=os.path.join(BASE_DIR,f"normal_{normal_tissue}_consensus_{WEIGHT_TAG}.tsv")

    if not os.path.exists(tumor_file):
        print(f"  Warning: tumor file {tumor_file} not found, skipping")
        continue
    if not os.path.exists(normal_file):
        print(f"  Warning: normal file {normal_file} not found, skipping")
        continue

    df_tumor=pd.read_csv(tumor_file,sep="\t",usecols=["TF","target","borda_score"])
    df_normal=pd.read_csv(normal_file,sep="\t",usecols=["TF","target","borda_score"])
    df_tumor["source"]="tumor"
    df_normal["source"]="normal"

    df_merged=pd.merge(
        df_tumor[["TF","target","borda_score"]],
        df_normal[["TF","target","borda_score"]],
        on=["TF","target"],how="outer",suffixes=("_tumor","_normal")
    )

    tumor_specific=df_merged[df_merged["borda_score_normal"].isna()].copy()
    tumor_specific["edge_type"]="Tumor_Specific"
    tumor_specific["borda_score"]=tumor_specific["borda_score_tumor"]

    normal_specific=df_merged[df_merged["borda_score_tumor"].isna()].copy()
    normal_specific["edge_type"]="Normal_Specific"
    normal_specific["borda_score"]=normal_specific["borda_score_normal"]

    shared=df_merged[df_merged["borda_score_tumor"].notna() & df_merged["borda_score_normal"].notna()].copy()
    shared["edge_type"]="Shared"
    shared["borda_score"]=shared["borda_score_tumor"]

    for df in [tumor_specific,normal_specific,shared]:
        df["tumor_type"]=tumor_type
        df["normal_tissue"]=normal_tissue

    def prepare_final_df(df,score_col="borda_score"):
        base_cols=["TF","target",score_col,"edge_type","tumor_type","normal_tissue"]
        if "borda_score_tumor" in df.columns and "borda_score_normal" in df.columns:
            base_cols.extend(["borda_score_tumor","borda_score_normal"])
        return df[base_cols].rename(columns={score_col:"borda_score"})

    tumor_specific_final=prepare_final_df(tumor_specific)
    normal_specific_final=prepare_final_df(normal_specific)
    shared_final=prepare_final_df(shared)

    current_result=pd.concat([tumor_specific_final,normal_specific_final,shared_final],ignore_index=True)
    all_results.append(current_result)

    current_result.to_csv(os.path.join(OUTPUT_DIR,f"{tumor_type}_complete_network_analysis.tsv"),sep="\t",index=False)
    tumor_specific_final.to_csv(os.path.join(OUTPUT_DIR,f"{tumor_type}_tumor_specific_network.tsv"),sep="\t",index=False)
    normal_specific_final.to_csv(os.path.join(OUTPUT_DIR,f"{tumor_type}_normal_specific_network.tsv"),sep="\t",index=False)
    shared_final.to_csv(os.path.join(OUTPUT_DIR,f"{tumor_type}_shared_network.tsv"),sep="\t",index=False)

    print(f"  Done {tumor_type}: Tumor_Specific={len(tumor_specific_final)}, Normal_Specific={len(normal_specific_final)}, Shared={len(shared_final)}")

if all_results:
    all_combined=pd.concat(all_results,ignore_index=True)
    all_combined.to_csv(os.path.join(OUTPUT_DIR,"ALL_TUMORS_complete_network_analysis.tsv"),sep="\t",index=False)

    all_tumor_specific=all_combined[all_combined["edge_type"]=="Tumor_Specific"]
    all_normal_specific=all_combined[all_combined["edge_type"]=="Normal_Specific"]
    all_shared=all_combined[all_combined["edge_type"]=="Shared"]

    all_tumor_specific.to_csv(os.path.join(OUTPUT_DIR,"ALL_TUMORS_tumor_specific_networks.tsv"),sep="\t",index=False)
    all_normal_specific.to_csv(os.path.join(OUTPUT_DIR,"ALL_TUMORS_normal_specific_networks.tsv"),sep="\t",index=False)
    all_shared.to_csv(os.path.join(OUTPUT_DIR,"ALL_TUMORS_shared_networks.tsv"),sep="\t",index=False)

    print("\n=== Summary ===")
    print(f"Total edges: {len(all_combined)}")
    print(f"Tumor-Specific edges: {len(all_tumor_specific)}")
    print(f"Normal-Specific edges: {len(all_normal_specific)}")
    print(f"Shared edges: {len(all_shared)}")

    print("\n=== Per-tumor statistics ===")
    stats=all_combined.groupby(["tumor_type","edge_type"]).size().unstack(fill_value=0)
    print(stats)
else:
    print("No tumor types were successfully processed")