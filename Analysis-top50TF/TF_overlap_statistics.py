import pandas as pd
import glob
import os
from itertools import combinations
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
hits_dir="/home/chenliqun/newFinalResults/HITS/HITS_results"
output_dir="./TF_overlap_results"
os.makedirs(output_dir,exist_ok=True)
files=glob.glob(os.path.join(hits_dir,"*HITS_results.tsv"))
print("Detected files:")
for f in files:
    print(f)
tf_lists={}
top50_lists={}
for f in files:
    cancer=os.path.basename(f).replace("_HITS_results.tsv","")
    df=pd.read_csv(f,sep="\t")
    tf_df=df[df["final_rank"]=="TF"].copy()
    tf_lists[cancer]=set(tf_df["gene"])
    tf_df=tf_df.sort_values("hub_score",ascending=False)
    top50=tf_df.head(50)["gene"].tolist()
    top50_lists[cancer]=set(top50)
    print(f"{cancer}: TF={len(tf_lists[cancer])}, Top50 extracted")
tf_universe=set()
for tf_set in tf_lists.values():
    tf_universe|=tf_set
N=len(tf_universe)
print("\nTotal TF universe:",N)
results=[]
for c1,c2 in combinations(top50_lists.keys(),2):
    set1=top50_lists[c1]
    set2=top50_lists[c2]
    overlap=set1&set2
    k=len(overlap)
    K=len(set1)
    M=len(set2)
    p_value=hypergeom.sf(k-1,N,K,M)
    results.append({"Cancer1":c1,"Cancer2":c2,"Top50_A":K,"Top50_B":M,"Overlap":k,"TF_universe":N,"p_value":p_value,"Shared_TFs":";".join(sorted(overlap))})
results_df=pd.DataFrame(results)
results_df["FDR"]=multipletests(results_df["p_value"],method="fdr_bh")[1]
results_df=results_df.sort_values("FDR")
supp_table=os.path.join(output_dir,"Supplementary_Table_TF_overlap_statistics.csv")
results_df.to_csv(supp_table,index=False)
print("\nSaved Supplementary Table:")
print(supp_table)
cancers=list(top50_lists.keys())
heatmap_matrix=pd.DataFrame(np.zeros((len(cancers),len(cancers))),index=cancers,columns=cancers)
for _,row in results_df.iterrows():
    c1=row["Cancer1"]
    c2=row["Cancer2"]
    value=-np.log10(row["FDR"]+1e-10)
    heatmap_matrix.loc[c1,c2]=value
    heatmap_matrix.loc[c2,c1]=value
np.fill_diagonal(heatmap_matrix.values,0)
plt.figure(figsize=(10,8))
sns.heatmap(heatmap_matrix,cmap="Reds",annot=True,fmt=".2f",linewidths=0.5,cbar_kws={"label":"-log10(FDR)"})
plt.title("TF Overlap Enrichment Between Cancer Types")
heatmap_file=os.path.join(output_dir,"TF_overlap_heatmap.pdf")
plt.tight_layout()
plt.savefig(heatmap_file)
plt.close()
print("\nSaved heatmap:")
print(heatmap_file)