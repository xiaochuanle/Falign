#!/usr/bin/env python
# coding: utf-8

import os
import re
import pandas as pd
import numpy as np
import random
import argparse

def loading_configuration(configuration):
    """
    1. fragnum per read
    2. ratio to proportion
    """
    conf=pd.read_csv(configuration,sep='\t')
    fns=conf.iloc[:,0].to_list()
    rts=conf.iloc[:,1].to_list()
    read_fragnum,probs=[],[]
    for i in range(len(conf)):
        s,p=fns[i],float(rts[i])
        fs,ls=int(s.split('-')[0]),int(s.split('-')[-1])
        st=list(range(fs,ls+1))
        pr=[p/len(st)]*len(st)
        read_fragnum.extend(st)
        probs.extend(pr)
    sum_probs=sum(probs)
    probs=probs if sum_probs==1 else [p/sum_probs for p in probs]
    return read_fragnum,probs

def preprocess_seed_fragment(fragments_file,fltlen,read_fragnum,probs):
    """
    1. filter fragments with odd bases, that can not be recognized by sequencing simulation tools
    2. randomly assign seed fragment and contact fragment number to each 3C read
    """
    fragments=pd.read_parquet(fragments_file)
    print("All virtual fragments:")
    print(fragments['fragment_length'].describe(),'\n')
    fragments=fragments[fragments['keep']=='pass']
    fragments=fragments[fragments['fragment_length']>=fltlen]
    fragments.reset_index(drop=True,inplace=True)
    fragments['id']=list(fragments.index)
    print("Passed virtual fragments:")
    print(fragments['fragment_length'].describe())
    
    read_num=len(fragments)
    fragnum_list = []
    for i in range(len(read_fragnum)):
        fn,pr=read_fragnum[i],probs[i]
        fc=int(read_num*pr)
        fragnum_list.extend([fn]*fc)
    random.shuffle(fragnum_list)
    if read_num > len(fragnum_list):
        fc=read_num-len(fragnum_list)
        supp_list=list(np.random.choice(read_fragnum,fc,p=probs))
        fragnum_list.extend(supp_list)
    fragments['fragment_number']=fragnum_list
    
    fragstrand_list=np.random.randint(0,2,read_num)
    fragstrand_list[fragstrand_list==1]=16
    fragments['strand']=fragstrand_list
    return(fragments)

def extract_id_dict(fragments,mainchrom):
    """
    1. each fragment has corresponding contiguous,intrachromosomal,interchromosomal fragments sets
    2. with linear distance increased, the interaction frequency gradually decreased
    3. extract fragments sets based on linear distance for simulation of proximity ligation
    """
    outchrom_id=list(fragments[fragments['chrom']!=mainchrom]['id'])
    mainchrom_id=set(fragments[fragments['chrom']==mainchrom]['id'])
    mainchrom_df=fragments[fragments['chrom']==mainchrom].groupby('inchrom_id')['id'].apply(set).reset_index()
    mainchrom_df.columns=['inchrom_id','cid']
    mainchrom_df['iid']=mainchrom_df['cid'].apply(lambda x: mainchrom_id.difference(x))
    contiguous_id_dict,inchrom_id_dict={},{}
    for i in mainchrom_df.index:
        inchrom_id=mainchrom_df['inchrom_id'][i]
        cid=mainchrom_df['cid'][i]
        iid=mainchrom_df['iid'][i]
        contiguous_id_dict[inchrom_id]=list(cid)
        inchrom_id_dict[inchrom_id]=list(iid)
    return(outchrom_id,contiguous_id_dict,inchrom_id_dict)

def sample_fragment_id(fragments,mainchrom,outchrom_id,inchrom_id_dict,contiguous_id_dict,outchrom,inchrom,contiguous):
    main_fragments=fragments[fragments['chrom']==mainchrom]
    #print(main_fragments['fragment_number'].describe())
    sfrag_id=main_fragments.loc[main_fragments.fragment_number==1,"id"].to_list()
    main_fragments=main_fragments[main_fragments.fragment_number>1].reset_index(drop=True)
    fragment_id=[]
    for i in sfrag_id:
        fragment_id.append([i])
    for i in main_fragments.index:
        seed_id=main_fragments['id'][i]
        fragnum=main_fragments['fragment_number'][i]
        seed_cid=main_fragments['inchrom_id'][i]
        
        cid=random.sample(contiguous_id_dict[seed_cid],min(len(contiguous_id_dict[seed_cid]),fragnum*contiguous))
        iid=list(np.random.choice(inchrom_id_dict[seed_cid],min(len(inchrom_id_dict[seed_cid]),fragnum*inchrom)))
        oid=list(np.random.choice(outchrom_id,fragnum*outchrom))
        id_pool=[seed_id]+cid+iid+oid
        id_pool=list(set(id_pool))
        while True:
            id_list=random.sample(id_pool,fragnum)
            fraglen_array=np.array(fragments.loc[id_list,'fragment_length'])
            shortfrag_idx=np.where(fraglen_array<200)
            no_shortfrag_end=len({0,fragnum-1} & set(shortfrag_idx[0].tolist())) == 0
            no_shortfrag_link=(np.diff(shortfrag_idx)==1).sum() == 0
            if no_shortfrag_end and no_shortfrag_link:
                break
        fragment_id.append(id_list)
    return(fragment_id)

def analysing_restriction_sites(enzyme_seq):
    """
    1. correct enzyme sequence like ^GATC, ^ denotes pointer
    2. left-cut before first base, like ^GATC
       right-cut after last base, like CATG^
       middle-cut in the sequence, like A^AGCTT
    """
    pointer_loc=[i.start() for i in re.finditer('\^',enzyme_seq)]
    if len(pointer_loc)==1: 
        pointer_loc=pointer_loc[0] 
        find_seq=enzyme_seq.replace('^','')
        seq_bytes=len(find_seq)
        cut_type="L" if pointer_loc==0 else "R" if pointer_loc==seq_bytes else "M"
        print("Correct enzyme sequence found.")
        print("find_seq","seq_bytes","pointer_loc","cut_type")
        print(find_seq,seq_bytes,pointer_loc,cut_type)
        return(find_seq,seq_bytes,pointer_loc,cut_type)
    else:
        print("Incorrect enzyme sequence found. Please check your input parameters.")

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

def transform_reapeating_sequence(seq):
    comp = {'a': 'A', 'g': 'G', 't': 'T', 'c': 'C'}
    seq_list = list(seq)
    comp_seq_list = [comp.get(base, base) for base in seq_list]
    comp_seq = ''.join(comp_seq_list)
    return comp_seq

def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N':'N'} # consider N base
    seq_list = list(seq)
    reverse_seq_list = reversed([comp.get(base, base) for base in seq_list])
    reverse_seq = ''.join(reverse_seq_list)
    return reverse_seq

#============== main =====================#
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate fragment proximity ligation, construct perfect 3C read library.')

    #-> required arguments
    parser.add_argument('-r', type=str, dest='ref_genome_file', required=True, help='path to reference genome file(unzip).')
    parser.add_argument('-e', type=str, dest='enzyme_seq', required=True, help='restriction enzyme sequence with pointer(^).')
    parser.add_argument('-f', type=str, dest='fragments_file', required=True, help='path to virtual fragment library.')
    
    #-> optional arguments
    parser.add_argument('-l', type=int, dest='fltlen', default=50, help='filter short fragments [50].')
    parser.add_argument('-c', type=str, dest='configuration', default=None, help='configure the ratio of reads with different number of fragments, as seen in readconf.tsv.')
    parser.add_argument('--outchrom', type=int, dest='outchrom', default=1, help='ratio of interchromosomal interaction frequency [1].')
    parser.add_argument('--inchrom', type=int, dest='inchrom', default=2, help='ratio of intrachromosomal interaction frequency [2].')
    parser.add_argument('--contiguous', type=int, dest='contiguous', default=4, help='ratio of contiguous interaction frequency [4].')
    parser.add_argument('--mainchrom', type=str, dest='mainchrom', default='', help='specify a main chromosome to generate library [Null].')
    args = parser.parse_args()

    if args.configuration:
        read_fragnum,probs=loading_configuration(args.configuration)
    else:
        read_fragnum,probs=[1,2,3,4,5],[0.2,0.2,0.2,0.2,0.2]
    print('Loading virtual fragment library file at:')
    fragments_file = os.path.abspath(args.fragments_file)
    print(fragments_file,'\n')
    fragments=preprocess_seed_fragment(fragments_file,args.fltlen,read_fragnum,probs)
    print('\n')
    print('Stratified sampling contact fragments')
    print('read_fragnum','ratio')
    for i in range(len(read_fragnum)):print(read_fragnum[i],probs[i])
    print('\n')
    chroms=[args.mainchrom] if args.mainchrom else list(sorted(set(fragments.chrom)))
    fragment_id_dict={}
    print("chrom","read_count")
    for mainchrom in chroms:
        outchrom_id,contiguous_id_dict,inchrom_id_dict=extract_id_dict(fragments,mainchrom)
        fragment_id_dict[mainchrom]=sample_fragment_id(fragments,mainchrom,outchrom_id,inchrom_id_dict,contiguous_id_dict,args.outchrom,args.inchrom,args.contiguous)
        print(mainchrom,len(fragment_id_dict[mainchrom]))
    print('\n')
    print('Extracting fragment sequence')
    find_seq,seq_bytes,pointer_loc,cut_type=analysing_restriction_sites(args.enzyme_seq)
    print('\n')
    ref_genome_file = os.path.abspath(args.ref_genome_file)
    print('Loading reference genome file at:')
    print(ref_genome_file)
    ref_genome = loading_reference(ref_genome_file)
    fragments['fragment_sequence']=fragments.apply(lambda x: ref_genome[x.chrom][x.start:x.end] if x.repratio==0 else transform_reapeating_sequence(ref_genome[x.chrom][x.start:x.end]),axis=1)
    fragments['fragment_sequence']=fragments['fragment_sequence'].apply(lambda x: x[seq_bytes-pointer_loc:] if cut_type=='L' else x[seq_bytes-pointer_loc:-pointer_loc])
    fragments['fragment_sequence']=fragments.apply(lambda x: x.fragment_sequence if x.strand==0 else reverse_complement(x.fragment_sequence),axis=1)
    print('\n')
    strcol=['chrom','start','end','strand','repratio']
    fragments[strcol]=fragments[strcol].astype(str)
    pid='p'+args.fragments_file.rsplit('_',1)[-1]
    fasta_file=fragments_file+'.fasta'
    print('Constructing perfect 3C read library')
    with open(fasta_file,'w') as fa:
        for mainchrom in chroms:
            read_num=len(fragment_id_dict[mainchrom])
            count=0
            for id_list in fragment_id_dict[mainchrom]:
                read=fragments.loc[id_list]
                fraginfo_list=read.apply(lambda x: '-'.join([x.chrom,x.start,x.end,x.strand,x.repratio]),axis=1)
                readinfo=';'.join(fraginfo_list)
                readseq=find_seq.join(read['fragment_sequence'])
                readid='>'+pid+'_'+mainchrom+'_'+str(count).zfill(len(str(read_num)))+' '+readinfo
                fa.write(readid+'\n')
                fa.write(readseq+'\n')
                count+=1
    print('Output file at:')
    print(fasta_file)
