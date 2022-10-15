#!/usr/bin/env python
# coding: utf-8

import os
import re
import pandas as pd
import numpy as np
import random
from joblib import Parallel, delayed
import argparse

def loading_configuration(configuration):
    """
    1. multifrag to cutting steps
    2. ratio to proportion
    """
    conf=pd.read_csv(configuration,sep='\t')
    mfs=conf.iloc[:,0].to_list()
    rts=conf.iloc[:,1].to_list()
    steps,probs=[],[]
    for i in range(len(conf)):
        s,p=mfs[i],float(rts[i])
        fs,ls=int(s.split('-')[0]),int(s.split('-')[-1])
        st=list(range(fs,ls+1))
        pr=[p/len(st)]*len(st)
        steps.extend(st)
        probs.extend(pr)
    sum_probs=sum(probs)
    probs=probs if sum_probs==1 else [p/sum_probs for p in probs]
    return steps,probs

def loading_reference(ref_genome_file):
    """
    1. filter unlocalized and unplaced (sometimes called 'random') and alternate sequences
    2. non-repeating sequence is shown in upper case;
       reapeating sequence is shown in lower case
    """
    ref_genome_pre = {}
    ## one chrom multiple sequence line
    for line in open(ref_genome_file):
        if line[0] == '>' :
            chrom = line.split()[0][1:]
            ref_genome_pre[chrom] = []
        else:
            ref_genome_pre[chrom].append(line.strip())
    ## one chrom one sequence line
    ref_genome = {}
    filter_chrom=['Un','random','alt']
    for chrom,valueL in ref_genome_pre.items():
        if not any(f in chrom for f in filter_chrom):
            ref_genome[chrom] = ''.join(valueL)
            print("Loading",chrom,len(ref_genome[chrom]))
    return ref_genome

def finding_restriction_sites(ref_genome,enzyme_seq):
    """
    1. correct enzyme sequence like ^GATC, ^ denotes pointer
    2. left-cut before first base, like ^GATC
       right-cut after last base, like CATG^
       middle-cut in the sequence, like A^AGCTT
    3. find restriction sites in the reference genome, independent of case
    """
    pointer_loc=[i.start() for i in re.finditer('\^',enzyme_seq)]
    if len(pointer_loc)==1: 
        pointer_loc=pointer_loc[0] 
        find_seq=enzyme_seq.replace('^','')
        seq_bytes=len(find_seq)
        cut_type="L" if pointer_loc==0 else "R" if pointer_loc==seq_bytes else "M"
        chrom_sites={}
        for chrom,ref_seq in ref_genome.items():
            start, end, sites = [0], [len(ref_seq)], []
            for site in re.finditer(find_seq,ref_seq,re.I):
                cutsite=site.start() if cut_type=="L" else site.end() if cut_type=="R" else (site.start()+pointer_loc)
                sites.append(cutsite)
            sites=start+sites+end
            chrom_sites[chrom]=sites
        print("Correct enzyme sequence found.")
        print("find_seq","seq_bytes","pointer_loc","cut_type")
        print(find_seq,seq_bytes,pointer_loc,cut_type)
        return(chrom_sites)
    else:
        print("Incorrect enzyme sequence found. Please check your input parameters.")

def simulating_genome_digestion(chrom_sites,seed):
    """
    1. random digestion at restriction sites
    2. digestion at adjacent 2 sites produce the smallest unit of virtual fragment
    3. digestion across multiple sites produce longer virtual fragment, that is more near to true circumstance
    """
    random.seed(seed)
    #steps,probs = [1,2,3],[0.25,0.5,0.25]
    new_chrom_sites={}
    for chrom,sites in chrom_sites.items():
        count,sitescount,new_sites=0,len(sites),[0]
        while (count < sitescount-1):
            step=int(np.random.choice(steps,1,p=probs))
            count=min(count+step,sitescount-1)
            site=sites[count]
            new_sites.append(site)
        new_chrom_sites[chrom]=new_sites
    return(new_chrom_sites)

def parallel_digestion(chrom_sites,depth,thread):
    parallel_id=list(range(depth))
    parallel_result=Parallel(n_jobs=thread)(delayed(simulating_genome_digestion)(chrom_sites,seed) for seed in parallel_id)
    new_chrom_sites_list=[]
    for rs in parallel_result:
        new_chrom_sites_list.append(rs)
    new_chrom_sites_dict=dict(zip(parallel_id,new_chrom_sites_list))
    return(new_chrom_sites_dict)

def extract_fragment_sequence(chrom,start,end):
    """
    1. extract sequence from reference genome
    2. calculate reapeating sequence ratio
    3. mark sequence with odd bases, like N S K
    """
    fragseq=ref_genome[chrom][start:end]
    allbase_count=len(fragseq)
    repbase_count=len(re.findall(r'a|g|c|t',fragseq))
    nonrepbase_count=len(re.findall(r'A|G|C|T',fragseq))
    normalbase_count=repbase_count+nonrepbase_count
    
    keep='pass' if allbase_count==normalbase_count else 'fail'
    repratio=0 if repbase_count==0 else max(1,int(repbase_count/allbase_count*100))
    return(keep,repratio)

def construct_fragment_library(new_chrom_sites):
    fragments=pd.DataFrame()
    for chrom,sites in new_chrom_sites.items():
        ## integrate fragments digestion information
        ssites=sites[0:-1]
        esites=sites[1:]
        chromf=pd.DataFrame({'chrom':chrom,'start':ssites,'end':esites})
        chromf['fragment_length']=chromf['end']-chromf['start']
        ## divide fragments into contiguous, intrachromosomal, interchromosomal intervals
        chromsize=sites[-1]
        interval=int(chromsize/10)
        interval_list=list(range(0,chromsize,interval))
        interval_list[-1]=chromsize
        interval_label=list(range(1,len(interval_list)))
        chromf['inchrom_id'] = pd.cut(chromf.start,interval_list,labels=interval_label,right=False)
        ## extract fragment sequence repeating information
        chromf[['keep','repratio']]=chromf.apply(lambda x:extract_fragment_sequence(x.chrom,x.start,x.end),axis=1,result_type="expand")
        fragments=pd.concat([fragments,chromf],ignore_index=True)

    outchrom_id=fragments.groupby(['chrom','inchrom_id'])['start'].count().reset_index()[['chrom','inchrom_id']]
    outchrom_id['outchrom_id']=outchrom_id.index+1
    fragments=fragments.merge(outchrom_id)
    return(fragments)

def parallel_construct_library(new_chrom_sites_dict,output_dir,depth,thread):
    parallel_id=list(range(depth))
    parallel_result=Parallel(n_jobs=thread)(delayed(construct_fragment_library)(new_chrom_sites) for pid,new_chrom_sites in new_chrom_sites_dict.items())
    pid=0
    for rs in parallel_result:
        file=output_dir+'/'+'fragments.parquet_'+str(pid)
        rs.to_parquet(file,index=False)
        pid+=1
        print(file)

#============== main =====================#
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate restriction enzyme digestion, extract virtual fragment library.')

    #-> required arguments
    parser.add_argument('-r', type=str, dest='ref_genome_file', required=True, help='path to reference genome file(unzip).')
    parser.add_argument('-e', type=str, dest='enzyme_seq', required=True, help='restriction enzyme sequence with pointer(^).')
    parser.add_argument('-o', type=str, dest='output_dir', required=True, help='output directory to fragment library.')
    
    #-> optional arguments
    parser.add_argument('-d', type=int, dest='depth', default=20, help='simulated depth(number of digestion) [20].')
    parser.add_argument('-t', type=int, dest='thread', default=10, help='multiple threads [10].')
    parser.add_argument('-c', type=str, dest='configuration', default=None, help='configure the ratio of fragments across several sites, as seen in fragconf.tsv.')
    args = parser.parse_args()

    if args.configuration:
        steps,probs=loading_configuration(args.configuration)
    else:
        steps,probs=[1,2,3],[0.25,0.5,0.25]
    print('Loading reference genome file at:')
    ref_genome_file = os.path.abspath(args.ref_genome_file)
    print(ref_genome_file)
    ref_genome = loading_reference(ref_genome_file)
    print('\n')
    print('Finding restriction sites')
    chrom_sites=finding_restriction_sites(ref_genome,args.enzyme_seq)
    print('\n')
    print('Simulating genome digestion')
    print('multifrag','ratio')
    for i in range(len(steps)):print(steps[i],probs[i])
    new_chrom_sites_dict=parallel_digestion(chrom_sites,depth=args.depth,thread=args.thread)
    print('\n')
    print('Constructing fragment library')
    output_dir = os.path.abspath(args.output_dir)
    parallel_construct_library(new_chrom_sites_dict,output_dir=output_dir,depth=args.depth,thread=args.thread)
    print('Done.')
