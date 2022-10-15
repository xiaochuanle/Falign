#!/usr/bin/env python
# coding: utf-8

### Match simulated 3C information and DeepSimulator results.
### Output final simulated Pore-C fastq.

## import function
import pandas as pd
import os
import argparse

def get_Csim_information(Csim_fasta_file):
    #print("Loading simulated 3C fasta at,")
    #print(Csim_fasta_file,"\n")
    rawid_list = []
    rawinfo_list = []
    with open(Csim_fasta_file,'r') as ca:
        for line in ca:
            if line.startswith(">"):
                rawid = line.split(chr(32))[0][1:]
                rawinfo = line.split(chr(32))[1]
                rawid_list.append(rawid)
                rawinfo_list.append(rawinfo)
    raw_info = pd.DataFrame({"rawid":rawid_list,
                             "rawinfo":rawinfo_list})
    return raw_info

def get_Dsim_information(DeepSimu_paf_file):
    #print("Loading DeepSimulator mapping.paf at,")
    #print(DeepSimu_paf_file,"\n")
    newid_list = []
    rawid_list = []
    with open(DeepSimu_paf_file,'r') as dp:
        for line in dp:
            newid = line.split('\t')[0]
            rawid = line.split('\t')[5]
            newid_list.append(newid)
            rawid_list.append(rawid)
    new_info = pd.DataFrame({"rawid":rawid_list,
                             "newid":newid_list})
    print("Checking and removing duplicated reads")
    print("Before:",len(new_info))
    new_info.drop_duplicates(["newid"],keep=False,inplace=True)
    new_info.drop_duplicates(["rawid"],keep=False,inplace=True)
    print("After:",len(new_info),"\n")
    return new_info

def get_matchup_information(DeepSimu_fastq_file,Csim_fastq_file):
    #print("Matching read id and inforamation",'\n')
    merge_info = pd.merge(new_info, raw_info,on="rawid")
    merge_info['readname'] = merge_info.apply(lambda x: '@'+' '.join([x.newid,x.rawid,x.rawinfo]),axis=1)
    merge_info.newid.to_csv("tmp.list",index=None,header=None)
    os.system("seqkit grep --quiet -j 20 -f %s %s -o %s "%("tmp.list",DeepSimu_fastq_file,"tmp.fastq"))

    rname_list = list(merge_info["readname"])
    count = 0
    with open("tmp.fastq",'r') as dq, open(Csim_fastq_file,'w') as cq:
        for line in dq:
            if line.startswith("@"):
                rname = rname_list[count]
                cq.write(rname)
                count += 1
            else:
                cq.write(line)
    os.system("rm %s %s "%("tmp.list","tmp.fastq"))

#============== main =====================#
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Match simulated 3C information and DeepSimulator results. Output final simulated Pore-C fastq.')
    
    #-> required arguments
    parser.add_argument('-ca', type=str, dest='Csim_fasta_file', required=True, help='path to simulated 3C fasta file')
    parser.add_argument('-dq', type=str, dest='DeepSimu_fastq_file', required=True, help='path to DeepSimulator fastq file')
    parser.add_argument('-dp', type=str, dest='DeepSimu_paf_file', required=True, help='path to DeepSimulator paf file')
    parser.add_argument('-cq', type=str, dest='Csim_fastq_file', required=True, help='output path to simulated Pore-C fastq file')
    args = parser.parse_args()
    
    ## get required information
    raw_info = get_Csim_information(args.Csim_fasta_file)
    new_info = get_Dsim_information(args.DeepSimu_paf_file)
    ## matchup and output
    get_matchup_information(args.DeepSimu_fastq_file,args.Csim_fastq_file)
