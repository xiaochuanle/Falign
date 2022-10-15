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

def load_Snake_par(Snake_par_file):
    Snake_par = pd.read_parquet(Snake_par_file)
    Snake_par = Snake_par[Snake_par['pass_filter'] == True]
    colselect = ['align_type','chrom', 'start', 'end','strand', 'read_name', 'read_length', 'read_start', 'read_end','mapping_quality','align_score','align_base_qscore','fragment_start', 'fragment_end']
    Snake_par = Snake_par[colselect]
    Snake_par["strand"] = Snake_par["strand"].apply(lambda x: "0" if x==True else "16")
    intcol = ['start', 'end','read_length', 'read_start', 'read_end','mapping_quality','align_score','align_base_qscore','fragment_start', 'fragment_end']
    Snake_par[intcol] = Snake_par[intcol].astype(int)
    return Snake_par

def load_batch_file(Snake_par_file_list):
    Snake_par_list = []
    for file in Snake_par_file_list:
        Snake_par = load_Snake_par(file)
        Snake_par_list.append(Snake_par)
    Snake_par = pd.concat(Snake_par_list)
    Snake_par.sort_values(by=["read_name","read_start","read_end"],inplace=True)
    Snake_par.reset_index(drop=True,inplace=True)
    del Snake_par_list
    gc.collect()
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
    map_info.sort_values(by=["read_name","read_start","read_end"],inplace=True)
    map_info.reset_index(drop=True,inplace=True)
    return map_info

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

def find_gap_overlap(fragment_dist):
    gc=len([i for i in fragment_dist if i >= gap_dist])
    f1=[i for i in range(len(fragment_dist)) if fragment_dist[i] <= -overlap]
    oc=len(f1)
    f2=[f+1 for f in f1]
    ofc=len(set(f1+f2))
    return gc,oc,ofc

def find_complete(gap_count,overlap_count):
    if gap_count==0 and overlap_count==0:
        return 'complete'
    elif overlap_count==0:
        return 'incomplete'
    else:
        return 'overlapped'

def read_stats(map_info):
    read_start_list = map_info.groupby("read_name")["read_start"].apply(list)
    read_end_list = map_info.groupby("read_name")["read_end"].apply(list)
    read_name_list = list(map_info.groupby("read_name")["read_name"].first())
    read_fragcount_list = list(map_info.groupby("read_name")["read_end"].count())
    read_info = pd.DataFrame({"read_name":read_name_list,"read_start":read_start_list,"read_end":read_end_list,"fragment_count":read_fragcount_list})
    read_info = read_info.reset_index(drop=True)

    read_info['fragment_dist']=read_info.apply(lambda x: np.array(x.read_start[1:]) - np.array(x.read_end[:-1]),axis=1)
    read_info[['gap_count','overlap_count','overlap_fcount']] = read_info.apply(lambda x: find_gap_overlap(x.fragment_dist),axis=1,result_type="expand")
    read_info["maptype"] = read_info.apply(lambda x: find_complete(x.gap_count, x.overlap_count), axis=1)
    return read_info

def output_stats(map_info,read_info):
    # frag stats, confident end: Both(LR),Single(L or R),None(UN)
    print("fragment stats:")
    total_fragcount = len(map_info)
    Bconf_fragcount = map_info['confident_end'].value_counts()['LR']
    Uconf_fragcount = map_info['confident_end'].value_counts()['UN']
    Sconf_fragcount = total_fragcount-Bconf_fragcount-Uconf_fragcount
    conf_fragcount = Bconf_fragcount + Sconf_fragcount
    
    Bconf_ratio = round(Bconf_fragcount/total_fragcount*100,2)
    Sconf_ratio = round(Sconf_fragcount/total_fragcount*100,2)
    Uconf_ratio = round(Uconf_fragcount/total_fragcount*100,2)
    conf_ratio = Bconf_ratio + Sconf_ratio
    
    map_info["length"] = map_info["end"]-map_info["start"]
    total_bases = map_info["length"].sum()
    Bconf_bases = map_info.groupby("confident_end")["length"].sum()['LR']
    Uconf_bases = map_info.groupby("confident_end")["length"].sum()['UN']
    Sconf_bases = total_bases-Bconf_bases-Uconf_bases
    conf_bases = Bconf_bases + Sconf_bases

    print("map_count","map_bases","confident_count","confident_ratio","confident_bases","both_count","single_count","none_count","both_ratio","single_ratio","none_ratio","both_bases","single_bases","none_bases")
    print(total_fragcount,total_bases,conf_fragcount,conf_ratio,conf_bases,Bconf_fragcount,Sconf_fragcount,Uconf_fragcount,Bconf_ratio,Sconf_ratio,Uconf_ratio,Bconf_bases,Sconf_bases,Uconf_bases)
    print("\n")
    
    # read stats, overlap or gap in read: complete, incomplete, overlapped
    print("read stats:")
    total_count = len(read_info)
    complete_count = read_info["maptype"].value_counts()['complete']
    incomplete_count = read_info["maptype"].value_counts()['incomplete']
    try:
        overlapped_count = read_info["maptype"].value_counts()['overlapped']
    except:
        overlapped_count = 0

    complete_ratio = round(complete_count/total_count*100,2)
    incomplete_ratio = round(incomplete_count/total_count*100,2)
    overlapped_ratio = round(overlapped_count/total_count*100,2)

    overlapped_fcount = read_info['overlap_fcount'].sum()
    overlapped_fratio = round(overlapped_fcount/total_fragcount*100,2)

    print("read_count","complete_count","incomplete_count","overlapped_count","complete_ratio","incomplete_ratio","overlapped_ratio","overlapped_fragment_count","overlapped_fragment_ratio")
    print(total_count,complete_count,incomplete_count,overlapped_count,complete_ratio,incomplete_ratio,overlapped_ratio,overlapped_fcount,overlapped_fratio)
    print("\n")
    
    # contact stats, pairwise contact in read
    print("contact stats:")
    contact_count = sum(read_info["fragment_count"] >= 2)
    read_info["pairwise_contact_count"] = read_info["fragment_count"].apply(lambda x: len(list(combinations(range(x),2))))
    pairwise_contact_count = sum(read_info["pairwise_contact_count"])
    print("contact_count","pairwise_contact_count")
    print(contact_count,pairwise_contact_count)
    
    # ouput stats df
    stats_cols=["map_count","map_bases","confident_count","confident_ratio","confident_bases","both_count","single_count","none_count","both_ratio","single_ratio","none_ratio","both_bases","single_bases","none_bases",\
    "read_count","complete_count","incomplete_count","overlapped_count","complete_ratio","incomplete_ratio","overlapped_ratio","overlapped_fragment_count","overlapped_fragment_ratio","contact_count","pairwise_contact_count"]
    stats_vals=[total_fragcount,total_bases,conf_fragcount,conf_ratio,conf_bases,Bconf_fragcount,Sconf_fragcount,Uconf_fragcount,Bconf_ratio,Sconf_ratio,Uconf_ratio,Bconf_bases,Sconf_bases,Uconf_bases,\
    total_count,complete_count,incomplete_count,overlapped_count,complete_ratio,incomplete_ratio,overlapped_ratio,overlapped_fcount,overlapped_fratio,contact_count,pairwise_contact_count]
    stats_df=pd.DataFrame({0:stats_cols,1:stats_vals})
    stats_df=stats_df.T
    stats_df.columns=stats_df.loc[0,]
    stats_df.drop(0,inplace=True)
    intcol=stats_df.columns.difference(["confident_ratio","both_ratio","single_ratio","none_ratio","complete_ratio","incomplete_ratio","overlapped_ratio","overlapped_fragment_ratio"])
    stats_df[intcol]=stats_df[intcol].astype(int)
    return stats_df

#============== main =====================#
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Evaluation of alignment tools for long noisy 3C data: Completeness and Contacts.')

    #-> required arguments
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument('-dir', type=str, dest='Snake_par_dir', default=None, help='path to Pore-C-Snakemake align_table directory')
    inputs.add_argument('-paf', type=str, dest='Cmap_paf_file', default=None, help='path to Falign alignment paf file')
    #-> optional arguments
    parser.add_argument('-end', type=int, dest='end_dist', default=20, help='distance threshold of confident end [20].')
    parser.add_argument('-gap', type=int, dest='gap_dist', default=50, help='distance threshold of gap between fragments [50].')
    parser.add_argument('-overlap', type=int, dest='overlap', default=10, help='distance threshold of overlap between fragments [10].')
    parser.add_argument('-save',type=str,dest='save',default='no',help='save details of each reads for further analysis [no].')
    parser.add_argument('-outdir',type=str,dest='outdir',default=None,help='path to save results. By default, output to the directory of alignment files.')
    args = parser.parse_args()

    # threshold
    end_dist = args.end_dist
    gap_dist = args.gap_dist
    overlap = args.overlap
    print('Setting thresholds for Stats:')
    print('end_dist','gap_dist','overlap')
    print(end_dist,gap_dist,overlap,'\n')
    ## load file
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

    # reads stats
    read_info = read_stats(map_info)

    ## output stats
    stats_df = output_stats(map_info,read_info)
    outdir=os.path.abspath(args.outdir) if args.outdir else outdir
    print('\n')
    print("Outputting Stats file at,")
    stats_df_file = outdir+"/"+"completeness.stats.txt"
    stats_df.to_csv(stats_df_file,index=None)
    print(stats_df_file)
    if args.save=='yes':
        frags_file = outdir+"/"+"fragments.completeness.details.txt"
        map_info.to_csv(frags_file,index=None)
        print(frags_file)
        reads_file = outdir+"/"+"reads.completeness.details.txt"
        usecol=["read_name","fragment_count","pairwise_contact_count","gap_count","overlap_count","overlap_fcount","maptype"]
        read_info[usecol].to_csv(reads_file,index=None)
        print(reads_file)
