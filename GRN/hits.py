import pandas as pd
import numpy as np
import os
from scipy.sparse import csr_matrix
import glob
import time
from datetime import datetime

def setup_logging(output_dir):
    """Set up log file"""
    os.makedirs(output_dir, exist_ok=True)
    log_file=os.path.join(output_dir,f"analysis_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt")
    return log_file

def log_message(log_file,message,print_to_console=True):
    """Write log message"""
    timestamp=datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    full_message=f"[{timestamp}] {message}"
    with open(log_file,'a',encoding='utf-8') as f:
        f.write(full_message+'\n')
    if print_to_console:
        print(full_message)

def build_network(edges_df):
    """Build weighted directed network matrix from TF-target-borda_score"""
    all_tfs=edges_df['TF'].unique()
    all_targets=edges_df['target'].unique()
    all_genes=np.unique(np.concatenate([all_tfs,all_targets]))
    n=len(all_genes)
    gene_to_idx={gene:idx for idx,gene in enumerate(all_genes)}
    idx_to_gene={idx:gene for gene,idx in gene_to_idx.items()}
    rows=[gene_to_idx[tf] for tf in edges_df['TF']]
    cols=[gene_to_idx[tg] for tg in edges_df['target']]
    weights=edges_df['borda_score'].values
    A=csr_matrix((weights,(rows,cols)),shape=(n,n))
    weighted_out_degree=A.sum(axis=1).A1
    weighted_in_degree=A.sum(axis=0).A1
    return A,all_genes,gene_to_idx,idx_to_gene,weighted_out_degree,weighted_in_degree

def hits_algorithm(edges_df,max_iter=100,tol=1e-6):
    """
    HITS algorithm: compute Hub and Authority scores for each gene
    Returns:
      - results: {gene:{hub_score,auth_score,...}}
      - convergence_data: convergence information for each iteration
    """
    A,all_genes,gene_to_idx,idx_to_gene,weighted_out_degree,weighted_in_degree=build_network(edges_df)
    n=len(all_genes)
    hub_scores=np.ones(n)
    auth_scores=np.ones(n)
    convergence_data=[]
    for iteration in range(max_iter):
        auth_new=A.T@hub_scores
        auth_norm=np.linalg.norm(auth_new)
        auth_new=auth_new/(auth_norm+1e-8)
        hub_new=A@auth_new
        hub_norm=np.linalg.norm(hub_new)
        hub_new=hub_new/(hub_norm+1e-8)
        hub_diff=np.linalg.norm(hub_new-hub_scores)
        auth_diff=np.linalg.norm(auth_new-auth_scores)
        hub_scores,auth_scores=hub_new,auth_new
        top_10_indices=np.argsort(hub_scores)[-10:][::-1]
        tf_count=sum(1 for idx in top_10_indices if weighted_out_degree[idx]>0)
        convergence_data.append({
            'iteration':iteration,
            'hub_diff':hub_diff,
            'auth_diff':auth_diff,
            'top10_tf_percent':tf_count/10*100,
            'top10_tf_count':tf_count
        })
        if hub_diff<tol and auth_diff<tol:
            break
    results={}
    for idx in range(n):
        gene=idx_to_gene[idx]
        results[gene]={
            'hub_score':hub_scores[idx],
            'auth_score':auth_scores[idx],
            'weighted_out_degree':weighted_out_degree[idx],
            'weighted_in_degree':weighted_in_degree[idx],
            'final_rank':'TF' if weighted_out_degree[idx]>0 else 'TG'
        }
    return results,convergence_data

def analyze_and_save_results(results,convergence_data,tumor_name,output_dir,log_file):
    """Analyze HITS results and save outputs"""
    os.makedirs(output_dir,exist_ok=True)
    results_df=pd.DataFrame.from_dict(results,orient='index')
    results_df=results_df.sort_values('hub_score',ascending=False)
    results_df['gene']=results_df.index
    results_df=results_df.reset_index(drop=True)
    top_10=results_df.head(10)
    top_20=results_df.head(20)
    top_50=results_df.head(50)
    tf_count_10=sum(1 for _,row in top_10.iterrows() if row['final_rank']=='TF')
    tf_count_20=sum(1 for _,row in top_20.iterrows() if row['final_rank']=='TF')
    tf_count_50=sum(1 for _,row in top_50.iterrows() if row['final_rank']=='TF')
    log_message(log_file,f"\n{'='*60}")
    log_message(log_file,f"HITS method - tumor: {tumor_name}")
    log_message(log_file,f"Total genes: {len(results_df)}")
    log_message(log_file,f"TF count: {(results_df['final_rank']=='TF').sum()}")
    log_message(log_file,f"TG count: {(results_df['final_rank']=='TG').sum()}")
    log_message(log_file,f"Iterations: {len(convergence_data)}")
    log_message(log_file,f"TF ratio in Top10: {tf_count_10}/10 = {tf_count_10/10*100:.1f}%")
    log_message(log_file,f"TF ratio in Top20: {tf_count_20}/20 = {tf_count_20/20*100:.1f}%")
    log_message(log_file,f"TF ratio in Top50: {tf_count_50}/50 = {tf_count_50/50*100:.1f}%")
    log_message(log_file,"Top 10 important genes (sorted by hub_score):")
    for i,(_,row) in enumerate(top_10.iterrows(),1):
        log_message(
            log_file,
            f"  {i:2d}. {row['gene']}: hub={row['hub_score']:.6f}, auth={row['auth_score']:.6f}, out_degree={row['weighted_out_degree']:.3f}, type={row['final_rank']}"
        )
    convergence_df=pd.DataFrame(convergence_data)
    convergence_file=os.path.join(output_dir,f"{tumor_name}_convergence.tsv")
    convergence_df.to_csv(convergence_file,sep='\t',index=False)
    results_file=os.path.join(output_dir,f"{tumor_name}_HITS_results.tsv")
    results_df.to_csv(results_file,sep='\t',index=False)
    return {
        'tumor_name':tumor_name,
        'method':'HITS',
        'total_genes':len(results_df),
        'tf_count':(results_df['final_rank']=='TF').sum(),
        'tg_count':(results_df['final_rank']=='TG').sum(),
        'iterations':len(convergence_data),
        'top10_tf_percent':tf_count_10/10*100,
        'top20_tf_percent':tf_count_20/20*100,
        'top50_tf_percent':tf_count_50/50*100
    }

def run_hits_only(input_dir,output_base_dir):
    """Run HITS method only"""
    hits_output_dir=os.path.join(output_base_dir,"HITS_results")
    os.makedirs(hits_output_dir,exist_ok=True)
    main_log_file=setup_logging(output_base_dir)
    pattern=os.path.join(input_dir,"*_tumor_specific_network.tsv")
    tumor_files=glob.glob(pattern)
    log_message(main_log_file,f"Found {len(tumor_files)} tumor network files")
    log_message(main_log_file,f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    summary_results=[]
    for file_path in tumor_files:
        base_name=os.path.basename(file_path)
        tumor_name=base_name.replace('_tumor_specific_network.tsv','')
        log_message(main_log_file,f"\n{'#'*80}")
        log_message(main_log_file,f"Processing tumor: {tumor_name}")
        log_message(main_log_file,f"File: {file_path}")
        try:
            df=pd.read_csv(file_path,sep='\t')
            required_cols=['TF','target','borda_score']
            if not all(col in df.columns for col in required_cols):
                log_message(main_log_file,f"Skipping {tumor_name}: missing required columns {required_cols}")
                continue
            log_message(main_log_file,f"\n--- Running HITS method ---")
            start_time=time.time()
            results,convergence_data=hits_algorithm(df)
            end_time=time.time()
            method_log_file=os.path.join(hits_output_dir,f"{tumor_name}_HITS_log.txt")
            analysis_result=analyze_and_save_results(
                results,convergence_data,tumor_name,hits_output_dir,method_log_file
            )
            analysis_result['computation_time']=end_time-start_time
            summary_results.append(analysis_result)
            log_message(main_log_file,f"HITS completed, runtime: {end_time-start_time:.2f} seconds")
        except Exception as e:
            log_message(main_log_file,f"Error processing {tumor_name}: {str(e)}")
            import traceback
            traceback.print_exc()
    summary_df=pd.DataFrame(summary_results)
    summary_file=os.path.join(output_base_dir,"HITS_summary.tsv")
    summary_df.to_csv(summary_file,sep='\t',index=False)
    log_message(main_log_file,f"\nAll HITS processing finished!")
    log_message(main_log_file,f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log_message(main_log_file,f"Summary saved to: {summary_file}")
    return summary_df

if __name__=="__main__":
    input_directory="/home/chenliqun/newFinalResults/tumor_specific"
    output_base_directory="/home/chenliqun/newFinalResults/HITS"
    summary_results=run_hits_only(input_directory,output_base_directory)
    print("\n"+"="*80)
    print("HITS summary:")
    print("="*80)
    if not summary_results.empty:
        avg_tf_percent=summary_results['top10_tf_percent'].mean()
        print(f"HITS - average TF ratio in Top10: {avg_tf_percent:.1f}%")
        avg_tf_percent2=summary_results['top50_tf_percent'].mean()
        print(f"HITS - average TF ratio in Top50: {avg_tf_percent2:.1f}%")
    else:
        print("No results (possibly no valid input files).")