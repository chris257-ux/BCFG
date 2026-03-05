import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Tuple,Dict,Optional,Set
import re
import sys
from datetime import datetime

def setup_logging(output_dir:Path):
    """Setup logging to file and console"""
    log_file=output_dir/f"borda_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    class DualLogger:
        def __init__(self,terminal,log_file):
            self.terminal=terminal
            self.log_file=log_file
        def write(self,message):
            self.terminal.write(message)
            with open(self.log_file,'a') as f:
                f.write(message)
        def flush(self):
            self.terminal.flush()
    sys.stdout=DualLogger(sys.stdout,log_file)
    print("="*60)
    print("Borda Consensus Analysis Log")
    print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Log file: {log_file}")
    print("="*60+"\n")
    return log_file

def extract_tumor_name(filename:str)->str:
    """Extract tumor type from filename"""
    return filename.split('_')[0]

def get_all_tumor_types(dir_path:Path)->Set[str]:
    """Get all tumor types from directory"""
    return {extract_tumor_name(f.name) for f in dir_path.glob("*.tsv") if f.is_file()}

def find_matching_files(tumor:str,scenic_dir:Path,deepsem_dir:Path)->Tuple[Path,Path]:
    """Match pySCENIC and deepSEM files for a tumor"""
    scenic_files=list(scenic_dir.glob(f"{tumor}_*.tsv"))
    deepsem_files=list(deepsem_dir.glob(f"processed_{tumor}_*.tsv"))
    if not scenic_files:
        raise FileNotFoundError(f"pySCENIC file not found: {tumor}")
    if not deepsem_files:
        raise FileNotFoundError(f"deepSEM file not found: processed_{tumor}")
    return scenic_files[0],deepsem_files[0]

def load_network_file(file_path:Path,algorithm:str)->pd.DataFrame:
    """Load network file"""
    try:
        df=pd.read_csv(file_path,sep='\t',header=0,usecols=['TF','Target','EdgeWeight'],dtype={'TF':str,'Target':str,'EdgeWeight':float})
    except Exception as e:
        with open(file_path) as f:
            sample=f.readlines()[:3]
        raise ValueError(f"Failed to read {algorithm} file: {file_path}\nError: {str(e)}\nFile preview:\n{''.join(sample)}")
    df=df.rename(columns={'TF':'TF','Target':'target','EdgeWeight':'weight'})
    if df['weight'].isnull().any():
        raise ValueError(f"{algorithm} file contains invalid weight values (NaN)")
    return df[['TF','target','weight']]

def calculate_borda_consensus(scenic_df:pd.DataFrame,deepsem_df:pd.DataFrame,weights:Tuple[float,float]=(0.5,0.5),output_ratio:float=1.0,deepsem_ratio:float=0.002,gold_standard:Optional[Path]=None)->Tuple[pd.DataFrame,Dict[str,int]]:
    """
    Modified Borda consensus calculation
    - pySCENIC: use all edges
    - deepSEM: use top deepsem_ratio edges
    - return result dataframe and statistics
    """
    scenic_df=scenic_df.sort_values('weight',ascending=False).reset_index(drop=True)
    scenic_df['rank']=scenic_df.index+1
    scenic_df['edge']=scenic_df['TF']+'|'+scenic_df['target']
    scenic_edges=set(scenic_df['edge'])
    scenic_rank=dict(zip(scenic_df['edge'],scenic_df['rank']))
    deepsem_top_k=max(int(len(deepsem_df)*deepsem_ratio),500)
    deepsem_df=deepsem_df.sort_values('weight',ascending=False).head(deepsem_top_k)
    deepsem_df['rank']=deepsem_df.index+1
    deepsem_df['edge']=deepsem_df['TF']+'|'+deepsem_df['target']
    deepsem_edges=set(deepsem_df['edge'])
    deepsem_rank=dict(zip(deepsem_df['edge'],deepsem_df['rank']))
    stats={'scenic_total':len(scenic_df),'deepsem_total':len(deepsem_df),'deepsem_used':deepsem_top_k,'total_candidate':len(scenic_edges.union(deepsem_edges)),'both_edges':len(scenic_edges.intersection(deepsem_edges))}
    all_edges=scenic_edges.union(deepsem_edges)
    N=len(all_edges)
    both_edges=scenic_edges.intersection(deepsem_edges)
    results=[]
    for edge in all_edges:
        r_scenic=scenic_rank.get(edge,N+1)
        r_deepsem=deepsem_rank.get(edge,N+1)
        score=weights[0]*(N-r_scenic)+weights[1]*(N-r_deepsem)
        results.append({'edge':edge,'TF':edge.split('|')[0],'target':edge.split('|')[1],'scenic_rank':r_scenic if r_scenic<=N else np.nan,'deepsem_rank':r_deepsem if r_deepsem<=N else np.nan,'borda_score':score,'consensus_type':'Both' if edge in both_edges else ('SCENIC' if edge in scenic_edges else 'deepSEM')})
    borda_df=pd.DataFrame(results).sort_values('borda_score',ascending=False)
    if output_ratio<1.0:
        keep_num=max(int(len(borda_df)*output_ratio),1)
        borda_df=borda_df.head(keep_num)
    if gold_standard:
        try:
            gold_edges=set(pd.read_csv(gold_standard,sep='\t').iloc[:,:2].apply(lambda x:f"{x[0]}|{x[1]}",axis=1))
            borda_df['in_gold_standard']=borda_df['edge'].isin(gold_edges)
            stats['gold_hits']=borda_df['in_gold_standard'].sum()
            stats['gold_total']=len(gold_edges)
        except Exception as e:
            print(f"Gold standard validation failed: {str(e)}")
    return borda_df,stats

def process_single_tumor(tumor:str,args)->Optional[str]:
    """Process a single tumor type"""
    try:
        print(f"\n=== Processing {tumor} ===")
        scenic_file,deepsem_file=find_matching_files(tumor,args.pyscenic_dir,args.deepsem_dir)
        print(f"  pySCENIC file: {scenic_file}")
        print(f"  deepSEM file: {deepsem_file}")
        scenic_df=load_network_file(scenic_file,"pySCENIC")
        deepsem_df=load_network_file(deepsem_file,"deepSEM")
        borda_df,stats=calculate_borda_consensus(scenic_df,deepsem_df,weights=(args.scenic_weight,args.deepsem_weight),output_ratio=args.output_ratio,deepsem_ratio=args.deepsem_ratio,gold_standard=args.gold_standard)
        output_file=args.output_dir/f"{tumor}_consensus_network.tsv"
        borda_df.to_csv(output_file,sep='\t',index=False)
        print(f"  Result saved: {output_file}")
        print("  Input statistics:")
        print(f"    - pySCENIC edges: {stats['scenic_total']}")
        print(f"    - deepSEM top edges: {stats['deepsem_used']}")
        print("  Candidate edge statistics:")
        print(f"    - Shared edges: {stats['both_edges']}")
        print(f"    - Total candidate edges: {stats['total_candidate']}")
        final_count=len(borda_df)
        both_count=(borda_df['consensus_type']=='Both').sum()
        scenic_only=(borda_df['consensus_type']=='SCENIC').sum()
        deepsem_only=(borda_df['consensus_type']=='deepSEM').sum()
        print("  Final result statistics:")
        print(f"    - Total retained edges: {final_count} ({(final_count/stats['total_candidate']):.1%})")
        print(f"    Both supported: {both_count} ({(both_count/final_count):.1%})")
        print(f"    SCENIC only: {scenic_only}")
        print(f"    deepSEM only: {deepsem_only}")
        if 'gold_hits' in stats:
            print("  Gold standard validation:")
            print(f"    - Gold edges covered: {stats['gold_hits']} ({(stats['gold_hits']/stats['gold_total']):.1%})")
            if both_count>0:
                both_gold=borda_df[borda_df['consensus_type']=='Both']['in_gold_standard'].mean()
                print(f"    - Gold ratio in Both edges: {both_gold:.1%}")
        return f"{tumor}: Success"
    except Exception as e:
        print(f"  Processing failed: {str(e)}")
        return f"{tumor}: Failed - {str(e)}"

def run_pipeline(args):
    """Run full analysis pipeline"""
    args.output_dir.mkdir(parents=True,exist_ok=True)
    log_file=setup_logging(args.output_dir)
    tumor_types=sorted(get_all_tumor_types(args.pyscenic_dir))
    print("\n"+"="*50)
    print(f"Weight settings: SCENIC={args.scenic_weight:.2f}, deepSEM={args.deepsem_weight:.2f}")
    print(f"Detected {len(tumor_types)} tumor types:")
    print("\n".join(f"  - {t}" for t in tumor_types))
    print(f"Candidate edge retention ratio: {args.output_ratio*100:.1f}%")
    print(f"deepSEM retention ratio: {args.deepsem_ratio*100:.1f}%")
    print("="*50)
    results=[]
    for tumor in tumor_types:
        results.append(process_single_tumor(tumor,args))
    print("\n=== Summary ===")
    print("\n".join(results))

if __name__=="__main__":
    parser=argparse.ArgumentParser(description='Integrate pySCENIC and deepSEM results')
    parser.add_argument('--pyscenic_dir',type=str,required=True,help='pySCENIC result directory')
    parser.add_argument('--deepsem_dir',type=str,required=True,help='deepSEM result directory')
    parser.add_argument('--output_dir',type=str,required=True,help='output directory')
    parser.add_argument('--gold_standard',type=str,default=None,help='gold standard file path (optional)')
    parser.add_argument('--output_ratio',type=float,default=1.0,help='final retention ratio (0.2=20%%,1.0=100%%)')
    parser.add_argument('--deepsem_ratio',type=float,default=0.002,help='deepSEM usage ratio')
    parser.add_argument('--scenic_weight',type=float,default=0.5,help='pySCENIC weight')
    parser.add_argument('--deepsem_weight',type=float,default=0.5,help='deepSEM weight')
    args=parser.parse_args()
    args.pyscenic_dir=Path(args.pyscenic_dir)
    args.deepsem_dir=Path(args.deepsem_dir)
    args.output_dir=Path(args.output_dir)
    args.gold_standard=Path(args.gold_standard) if args.gold_standard else None
    run_pipeline(args)