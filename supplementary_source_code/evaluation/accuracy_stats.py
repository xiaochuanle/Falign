#!/usr/bin/env python
# coding: utf-8

import argparse
import os
import pandas as pd
import numpy as np
import datatable as dt
import gc
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

def get_fragment_information(Csim_fastq_file):
    readid_list, chrom_list, chrsite_list, endsite_list, strand_list, fragid_list = [],[],[],[],[],[]
    with open(Csim_fastq_file,'r') as cq:
        for line in cq:
            if line.startswith("@"):
                readid = line.split(chr(32))[0][1:]
                raw_info = line.split(chr(32))[2].strip()
                frag_info_list = raw_info.split(';')
                frag_num = len(frag_info_list)
                readid_list.extend([readid]*frag_num)
                fragid_list.extend(list(range(0,frag_num)))
                chrom = [x.split('-')[0] for x in frag_info_list]
                chrom_list.extend(chrom)
                chrsite = [x.split('-')[1] for x in frag_info_list]
                chrsite_list.extend(chrsite)
                endsite = [x.split('-')[2] for x in frag_info_list]
                endsite_list.extend(endsite)
                strand = [x.split('-')[3] for x in frag_info_list]
                strand_list.extend(strand)
                
    raw_info = pd.DataFrame({"read_name":readid_list,
                             "fidx":fragid_list,
                             "strand":strand_list,
                             "chrom":chrom_list,
                             "start":chrsite_list,
                             "end":endsite_list})
    raw_info.sort_values(by=["read_name","fidx"],inplace=True)
    raw_info.reset_index(drop=True,inplace=True)
    ridx = raw_info[['read_name']].drop_duplicates().reset_index(drop=True)
    ridx['ridx'] = ridx.index
    raw_info = raw_info.merge(ridx,on='read_name')
    intcol = ['fidx','start', 'end','ridx']
    raw_info[intcol] = raw_info[intcol].astype(int)
    return raw_info

def load_Snake_par(Snake_par_file):
    Snake_par = pd.read_parquet(Snake_par_file)
    Snake_par = Snake_par[Snake_par['pass_filter'] == True]
    colselect = ['align_type','chrom', 'start', 'end','strand', 'read_name', 'read_length', 'read_start', 'read_end','mapping_quality','align_score','align_base_qscore','fragment_start', 'fragment_end']
    Snake_par = Snake_par[colselect]
    Snake_par["strand"] = Snake_par["strand"].apply(lambda x: "0" if x==True else "16")
    intcol = ['start', 'end','read_length', 'read_start', 'read_end','mapping_quality','align_score','align_base_qscore','fragment_start', 'fragment_end']
    Snake_par[intcol] = Snake_par[intcol].astype(int)
    strcol = ["read_name","chrom","strand"]
    Snake_par[strcol] = Snake_par[strcol].astype(str)
    return Snake_par

def load_batch_file(Snake_par_file_list):
    Snake_par_list = []
    for file in Snake_par_file_list:
        Snake_par = load_Snake_par(file)
        Snake_par_list.append(Snake_par)
    Snake_par = pd.concat(Snake_par_list)
    Snake_par.sort_values(by=["read_name","read_start","read_end"],inplace=True)
    Snake_par.reset_index(drop=True,inplace=True)
    return Snake_par

def get_mapping_information(Cmap_paf_file):
    awk_cmd = "awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $8, $9, $12, $15, $16}' " + Cmap_paf_file
    #print(awk_cmd)
    dt_df = dt.fread(cmd = awk_cmd,header=False)
    map_info = dt_df.to_pandas()
    del dt_df
    gc.collect()
    map_info.columns=["read_name","read_length","read_start","read_end","strand","chrom","start","end","mapping_quality","fragment_start","fragment_end"]
    map_info["strand"] = map_info["strand"].apply(lambda x: "0" if x=="+" else "16")
    map_info['fragment_start']=map_info['fragment_start'].apply(lambda x: x[5:])
    map_info['fragment_end']=map_info['fragment_end'].apply(lambda x: x[5:])
    intcol = ['read_length','read_start', 'read_end','start', 'end','mapping_quality','fragment_start', 'fragment_end']
    map_info[intcol] = map_info[intcol].astype(int)
    strcol = ["read_name","chrom","strand"]
    map_info[strcol] = map_info[strcol].astype(str)
    map_info.sort_values(by=["read_name","read_start","read_end"],inplace=True)
    map_info.reset_index(drop=True,inplace=True)
    return map_info

## confident end
def cal_virtual_dist(map_info,end_dist):
    ## calculate virtual distance
    map_info['Ldist_virtual'] = abs(map_info['fragment_start'] - map_info['start'])
    map_info['Rdist_virtual'] = abs(map_info['fragment_end'] - map_info['end'])
    ## judge confident end
    Lbool = (map_info['Ldist_virtual'] <= end_dist)
    Rbool = (map_info['Rdist_virtual'] <= end_dist)
    map_info['confident_end'] = 'UN'
    map_info.loc[(Lbool),"confident_end"] = 'L'
    map_info.loc[(Rbool),"confident_end"] = 'R'
    map_info.loc[(Lbool&Rbool),"confident_end"] = 'LR'
    ## assign ridx, fidx
    ridx = map_info.groupby('read_name')['strand'].count().reset_index()
    ridx['ridx'] = ridx.index
    fidx = []
    for i in ridx['strand']:
        rfidx = [v for v in range(i)]
        fidx.extend(rfidx)
    ridx_fidx = map_info[['read_name']].merge(ridx[['read_name','ridx']])
    ridx_fidx['fidx'] = fidx
    map_info[['ridx','fidx']] = ridx_fidx[['ridx','fidx']]
    return map_info

## match end
def get_frags_stats(raw_info,map_info,end_dist):
    raw_info = raw_info.copy()
    map_info = map_info.copy()
    raw_info = raw_info.rename(columns={"start":"start_raw","end":"end_raw","ridx":"ridx_raw","fidx":"fidx_raw"})
    map_info = map_info.rename(columns={"start":"start_map","end":"end_map","ridx":"ridx_map","fidx":"fidx_map"})
    raw_col = raw_info.columns
    map_col = ["read_name","chrom","strand","start_map","end_map",'mapping_quality',"ridx_map","fidx_map"]
    #map_col = ["read_name","chrom","strand","start_map","end_map",'mapping_quality',"align_score","align_base_qscore","ridx_map","fidx_map"]

    frags_stats = pd.merge(raw_info, map_info[map_col],on=["read_name","chrom","strand"])
    ## calculate distance
    frags_stats['Ldist_true'] = abs(frags_stats.start_raw - frags_stats.start_map)
    frags_stats['Rdist_true'] = abs(frags_stats.end_raw - frags_stats.end_map)
    ## judge match end
    Lbool = (frags_stats['Ldist_true'] <= end_dist)
    Rbool = (frags_stats['Rdist_true'] <= end_dist)
    frags_stats['true_end'] = 'UN'
    frags_stats.loc[(Lbool&Rbool),"true_end"] = 'LR'
    ## filter and keep unique candidate map, others are considered as incorrect map
    frags_stats = frags_stats[frags_stats['true_end']=='LR']
    frags_stats.sort_values(['read_name','fidx_raw','true_end'],inplace=True)
    frags_stats.reset_index(drop=True,inplace=True)
    frags_stats.drop_duplicates(['read_name','fidx_raw'],inplace=True)
    frags_stats.reset_index(drop=True,inplace=True)

    ## integrate stats
    frags_stats['type'] = 'correct'
    int_col = list(frags_stats.select_dtypes(include='number').columns)
 
    unmap = raw_info.append(frags_stats[raw_col])
    unmap.drop_duplicates(subset=['read_name','ridx_raw','fidx_raw'],keep=False,inplace=True)
    unmap['type'] = 'unmap'
    
    incorrect = map_info[map_col].append(frags_stats[map_col])
    incorrect.drop_duplicates(subset=['read_name','ridx_map','fidx_map'],keep=False,inplace=True)
    incorrect['type'] = 'incorrect'
    
    frags_stats = pd.concat([frags_stats,unmap,incorrect])
    frags_stats[int_col] = frags_stats[int_col].fillna(-1)
    frags_stats[int_col] = frags_stats[int_col].astype(int)
    frags_stats['true_end'] = frags_stats['true_end'].fillna('UN')
    frags_stats.sort_values(['read_name','fidx_raw','fidx_map'],inplace=True)
    frags_stats.reset_index(drop=True,inplace=True)
    frags_stats['type'] = pd.Categorical(frags_stats['type'], ["correct","unmap","incorrect"])
    
    ## fragments length
    frags_stats['length'] = frags_stats.apply(lambda x: (x.end_map - x.start_map) if x.type=='incorrect' else (x.end_raw - x.start_raw),axis=1)
    return frags_stats

def get_reads_stats(frags_stats):
    frags_stats = frags_stats.copy()
    ## sub_read_type, detail
    reads_stats = frags_stats.groupby('read_name')['type'].apply(list).reset_index()
    reads_stats['sub_read_type'] = reads_stats['type'].apply(lambda x: '_'.join(sorted(set(x))))
    ## read_type, 3 main type, full/partial/unmatch
    reads_stats['read_type'] = 'unmatch'
    partial = ['correct_incorrect_unmap','correct_incorrect','correct_unmap']
    partial_bool = reads_stats['sub_read_type'].apply(lambda x: x in partial)
    reads_stats.loc[partial_bool,"read_type"] = 'partial'
    full_bool = (reads_stats['sub_read_type'] == 'correct')
    reads_stats.loc[full_bool,"read_type"] = 'full'
    ## raw chrom count, raw frag count
    raw_reads_stats = frags_stats[frags_stats['fidx_raw'] != -1].groupby('read_name')['chrom'].apply(list).reset_index()
    raw_reads_stats['chrom_count'] = raw_reads_stats['chrom'].apply(lambda x: len(set(x)))
    raw_reads_stats['frag_count'] = raw_reads_stats['chrom'].apply(lambda x: len(x))
    ## add count
    reads_stats = reads_stats.merge(raw_reads_stats,on='read_name')
    ## pairwise count
    reads_stats["pwise_count"]=reads_stats["frag_count"].apply(lambda x: len(list(combinations(range(x),2))))
    reads_stats['mfrag_count']=reads_stats['type'].apply(lambda x: len(x)-x.count("unmap"))
    reads_stats["mpwise_count"]=reads_stats["mfrag_count"].apply(lambda x: len(list(combinations(range(x),2))))
    reads_stats['cfrag_count']=reads_stats['type'].apply(lambda x: x.count("correct"))
    reads_stats["cpwise_count"]=reads_stats["cfrag_count"].apply(lambda x: len(list(combinations(range(x),2))))
    return reads_stats

def output_stats(raw_info,map_info,frags_stats,reads_stats):
    # frag stats, P=C/C+I ; S= C/C+U ## v7 changed.
    print("fragment stats:")
    raw_count = len(raw_info)
    map_count = len(map_info)
    correct_count = frags_stats["type"].value_counts()['correct']
    unmap_count = frags_stats["type"].value_counts()['unmap']
    incorrect_count = frags_stats["type"].value_counts()['incorrect']

    Precision = round(correct_count/(correct_count + incorrect_count)*100,2)
    Sensitivity = round(correct_count/(correct_count + unmap_count)*100,2)
    #Sensitivity = round(correct_count/(correct_count + incorrect_count + unmap_count)*100,2)
    
    print("raw_count","map_count","correct_count","unmap_count","incorrect_count","Precision(%)","Sensitivity(%)")
    print(raw_count,map_count,correct_count,unmap_count,incorrect_count,Precision,Sensitivity)
    print("\n")
    
    raw_length = sum(raw_info["end"] - raw_info["start"])
    map_length = sum(map_info["end"] - map_info["start"])

    correct_frags_stats = frags_stats[frags_stats['type']=='correct']
    correct_raw_length = sum(correct_frags_stats["end_raw"] - correct_frags_stats["start_raw"])
    correct_map_length = sum(correct_frags_stats["end_map"] - correct_frags_stats["start_map"])
    correct_raw_ratio = round(correct_raw_length/(raw_length)*100,2)
    correct_map_ratio = round(correct_map_length/(map_length)*100,2)

    print("raw","map")
    print(raw_length,map_length)
    print(correct_raw_length,correct_map_length)
    print(correct_raw_ratio,correct_map_ratio)
    print("\n")
    
    # read stats
    print("read stats:")
    full_match = reads_stats['read_type'].value_counts()['full']
    partial_match = reads_stats['read_type'].value_counts()['partial']
    try:
        unmatch = reads_stats['read_type'].value_counts()['unmatch']
    except:
        unmatch = 0
    total = len(reads_stats)

    full_ratio = round(full_match/total*100,2)
    partial_ratio = round(partial_match/total*100,2)
    unmatch_ratio = round(unmatch/total*100,2)

    print("full","partial","unmatch","total")
    print(full_match,partial_match,unmatch,total)
    print(full_ratio,partial_ratio,unmatch_ratio)
    print('\n')

    # pairwise stats
    print("pairwise stats:")
    raw_pwise,map_pwise,correct_pwise=reads_stats['pwise_count'].sum(),reads_stats['mpwise_count'].sum(),reads_stats['cpwise_count'].sum()
    print("raw","map","correct")
    print(raw_pwise,map_pwise,correct_pwise)
    print('\n')

    # ouput stats df
    pd.set_option('display.float_format', lambda x: '%.2f' % x)
    stats_cols=["raw_count","map_count","correct_count","unmap_count","incorrect_count","Precision","Sensitivity",\
    "raw_bases","correct_raw_bases","correct_raw_bases_ratio","map_bases","correct_map_bases","correct_map_bases_ratio",\
    "read_count","full_read_count","partial_read_count","unmatch_read_count","full_read_ratio","partial_read_ratio","unmatch_read_ratio",\
    "raw_pairwise_count","map_pairwise_count","correct_pairwise_count"]
    stats_vals=[raw_count,map_count,correct_count,unmap_count,incorrect_count,Precision,Sensitivity,\
    raw_length,correct_raw_length,correct_raw_ratio,map_length,correct_map_length,correct_map_ratio,\
    total,full_match,partial_match,unmatch,full_ratio,partial_ratio,unmatch_ratio,\
    raw_pwise,map_pwise,correct_pwise]
    stats_df=pd.DataFrame({0:stats_cols,1:stats_vals})
    stats_df=stats_df.T
    stats_df.columns=stats_df.loc[0,]
    stats_df.drop(0,inplace=True)
    intcol=stats_df.columns.difference(["Precision","Sensitivity","correct_raw_bases_ratio","correct_map_bases_ratio","full_read_ratio","partial_read_ratio","unmatch_read_ratio"])
    stats_df[intcol]=stats_df[intcol].astype(int)
    return stats_df

#============== main =====================#
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Evaluation of alignment tools for long noisy 3C data: Precision and Sensitivity.')

    #-> required arguments
    parser.add_argument('-cq', type=str, dest='Csim_fastq_file', required=True, help='path to simulated long noisy 3C fastq file')
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument('-dir', type=str, dest='Snake_par_dir', default=None, help='path to Pore-C-Snakemake align_table directory')
    inputs.add_argument('-paf', type=str, dest='Cmap_paf_file', default=None, help='path to Falign alignment paf file')
    #-> optional arguments
    parser.add_argument('-end', type=int, dest='end_dist', default=20, help='distance threshold of confident/match end [20].')
    parser.add_argument('-save',type=str,dest='save',default='no',help='save details of each reads for further analysis [no].')
    parser.add_argument('-outdir',type=str,dest='outdir',default=None,help='path to save results. By default, output to the directory of alignment files.')
    args = parser.parse_args()

    # threshold
    end_dist = args.end_dist
    print('Setting thresholds for Stats:')
    print('end_dist')
    print(end_dist,'\n')
    ## load file
    # raw fragments info
    print("Loading simulated long noisy 3C fastq file at,")
    Csim_fastq_file = os.path.abspath(args.Csim_fastq_file)
    print(Csim_fastq_file,'\n')
    raw_info = get_fragment_information(Csim_fastq_file)
    # map fragments info
    if args.Snake_par_dir:
        Snake_par_dir = os.path.abspath(args.Snake_par_dir)
        Snake_par_file_list = [x for x in os.listdir(Snake_par_dir) if x.endswith("pore_c.parquet")]
        print("Loading Pore-C-Snakemake alignment parquet file at,")
        print(Snake_par_dir,'\n')
        #for file in sorted(Snake_par_file_list):print(file)
        #print('\n')
        map_info = load_batch_file(Snake_par_file_list)
        outdir=Snake_par_dir
    else:
        print("Loading Falign alignment paf file at,")
        Cmap_paf_file = os.path.abspath(args.Cmap_paf_file)
        print(Cmap_paf_file,'\n')
        map_info = get_mapping_information(Cmap_paf_file)
        outdir=os.path.dirname(Cmap_paf_file)
    map_info = cal_virtual_dist(map_info,end_dist)

    # fragments stats
    frags_stats = get_frags_stats(raw_info,map_info,end_dist)
    # reads stats
    reads_stats = get_reads_stats(frags_stats)
    print('\n')
    ## output stats
    stats_df = output_stats(raw_info,map_info,frags_stats,reads_stats)
    outdir=os.path.abspath(args.outdir) if args.outdir else outdir

    print("Outputting Stats file at,")
    stats_df_file = outdir+"/"+"accuracy.stats.txt"
    stats_df.to_csv(stats_df_file,index=None)
    print(stats_df_file)
    if args.save=='yes':
        frags_file = outdir+"/"+"fragments.accuracy.details.txt"
        frags_stats.to_csv(frags_file,index=None)
        print(frags_file)
        reads_file = outdir+"/"+"reads.accuracy.details.txt"
        reads_stats.to_csv(reads_file,index=None)
        print(reads_file)
