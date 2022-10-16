# 3Csim
long noisy chromosome conformation capture (3C) reads simulator

# <a name="S-table-of-contents"></a> Table of Contents

* [Introduction](#S-introduction)
* [Installation](#S-installation)
* [Usage](#S-usage)

# <a name="S-introduction"></a> Introduction
Fragmented long noisy reads (FLNRs), such as Pore-C, are advantageous in identifying high-order 3D structures and structural variants. Pore-C reads are composed of several chromatin fragments with lengths from tens of bps to thousands of bps, and contain nanopore sequencing errors.

Here, we propose 3Csim pipeline to generate simulated Pore-C datasets for algorithms evaluation. We first simulated long 3C datasets based on reference genome, restriction enzyme, fragment length, fragment number and the depth of sequencing. Then, we modified `DeepSimulator` to directly generate errors mimicking Oxford Nanopore Technologies (ONT) sequencing based on the concatenated sequence. It is practical to generate other simulated long nosiy 3C datasets by collocating other simulators, such as `pbsim` for PacBio.

# <a name="S-installation"></a> Installation
## <a name="SS-download-and-install-dependencies"></a> Download and Install dependencies
```shell
$ chmod u+x install.sh
$ ./install.sh
$ chmod u+x DeepSimulator/*.sh
```

# <a name="S-usage"></a> Usage
### </a> 1) Generate fragments and long 3C reads
```shell
$ python fragment_library.py -r ${reference_fasta} -e ${enzyme_sequence} -o ${out_directory} 
$ python 3C_read_library.py -r ${reference_fasta} -e ${enzyme_sequence} -f ${out_directory}/fragments.parquet_0 
```

* Required arguments
```shell
-r path to reference genome file(unzip).
-e restriction enzyme sequence with pointer(^).
-o output directory.
```
* Optional arguments
```shell
## for fragments
-d simulated depth(number of digestion) [20].
-c configure the ratio of fragments across several sites,
   as seen in fragconf.tsv.
## for long 3C reads
-l filter short fragments [50].
-c configure the ratio of reads with different number of fragments,
   as seen in readconf.tsv.
--outchrom   ratio of interchromosamal interaction frequency [1].
--inchrom    ratio of intrachromosomal interaction frequency [2].
--contiguous ratio of contiguous interaction frequency [4].
--mainchrom  specify a main chromosome to generate library [Null].
```

### </a> 2) Micmik ONT errors
```shell
$ ./DeepSimulator/deep_simulator.PoreC.sh -i ${out_directory}/fragments.parquet_0.fasta -o ${out_directory}
```
* Show details on DeepSimulator arguments
```shell
$ ./DeepSimulator/deep_simulator.PoreC.sh
```
### </a> Example
Generate A. thaliana Pore-C datasets using DpnII. 
```shell
$ wget http://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
$ gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
$ python fragment_library.py -r Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -e ^GATC -o Ara_DpnII
$ python 3C_read_library.py -r Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -e ^GATC -f Ara_DpnII/fragments.parquet_0.fasta
$ ./DeepSimulator/deep_simulator.PoreC.sh -i Ara_DpnII/fragments.parquet_0.fasta -o Ara_DpnII
```


### </a> Reference
Li, Y. et al. DeepSimulator: a deep simulator for Nanopore sequencing. Bioinformatics. 34, 2899-2908 (2018).
