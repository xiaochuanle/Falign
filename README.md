# <a name="S-table-of-contents"></a> Table of Contents

* [Introduction](#S-introduction)
* [Tools](#S-tools)
* [Installation](#S-installation)
* [Usage](#S-usage)
* [Use `u4falign` for manipulating alignment results](#S-u4falign)
* [Difference between `read-SAM` and `frag-SAM`](#S-sam-diff)

# <a name="S-current-version"></a> Current Version

'0.0.2'

# <a name="S-introduction"></a> Introduction

`Falign` is a sequence alignment toolkit for long noisy chromosome conformation capture (3C) reads, such as Pore-C.
`Falign` is written in C and C++ programming language.

# <a name="S-tools"></a> Tools

Three tools are released together in this toolkit.

* `falign`. The alignment tool.
* `falign_ngf`. Another alignment tool. It is used for benchmark only.
* `u4falign`. A utility for manipulating `SAM` or `PAF` mapping results.
* The directory `supplementary_source_code` contains scripts for generating simulated Pore-C reads and analyzing mapping results.

# <a name="S-installation"></a> Installation

## <a name="SS-download"></a> Download files

```shell
$ git clone https://github.com/xiaochuanle/Falign.git
```

## <a name="SS-install-from-binary"></a> Install from executable binaries

```shell
$ cd Falign/release
$ tar xzvf Falign-0.0.2-20230213-Linux-amd64.tar.gz
$ cd Falign/Linux-amd64/bin/
$ export PATH=$PATH:$(pwd)
```
The last command `export PATH=$PATH:$(pwd)` is used for adding the path of `falign` to the system `PATH` so that you don't have to type the full path of `falign` (such as `/data3/cy/map_test/Falign/Linux-amd64/bin/falign`) every time you used `falign`.

## <a name="SS-install-from-source-codes"></a> Install from source codes

```shell
$ cd Falign/src
$ make -j
$ cd ../Linux-amd64/bin/
$ export PATH=$PATH:$(pwd)
```

# <a name="S-usage"></a> Usage

Decompress the sample data and decompress it.
```shell
$ cd Falign
$ tar xzvf sample-data.tar.gz
$ cd sample-data
```

We provide one sample reference `ara_2_4.fa` in the `./sample-data/ref` directory:
```shell
$ ls ref
ara_2_4.fa
```
and provide three sample reads in the `./sample-data/reads/` directory:
```shell
$ ls reads
ara_reads_1.fq ara_reads_2.fq ara_reads_3.fq
```

## <a name="SS-quick-start"></a> Quick Start

By default, `falign` outputs mapping results in `SAM` format:

```shell
$ falign -num_threads 48 ^GATC ref/ara_2_4.fa \
	reads/ara_reads_1.fq \
    reads/ara_reads_2.fq \
    reads/ara_reads_3.fq > map.sam
```

Users can output the mapping results in `PAF` format by using the `-outfmt paf` option:
```shell
$ falign -num_threads 48 -outfmt paf ^GATC ref/ara_2_4.fa \
	reads/ara_reads_1.fq \
    reads/ara_reads_2.fq \
    reads/ara_reads_3.fq > map.paf
```

## <a name="some-comments-on-usage"></a> Some comments on usage

### </a> 1. The restriction enzyme sequence

In the running commands above, `^GATC` is the sequence of the DpnII restriction enzyme used for generating the reads. `falign` provides the sequences for all familiar restriction enzymes (to see the following information, just type the `falign` command) so that you don't bother to lookup for them in other places:
```shell
*** Examples of familiar <enzyme_seq>:
Enzyme_Name             Enzyme_Seq
DpnII                   ^GATC
HindIII                 A^AGCTT
NcoI                    C^CATGG
NlaIII                  CATG^
```

### </a> The input format of reads

Besides the way in the examples above, `falign` also accepts reads input in the following way:

* Read list. You can list the paths of all the reads in a file:
```shell
$ cat read_list.lst
/Users/chenying/Desktop/mwj/temp/Falign/sample-data/reads/ara_reads_1.fq
/Users/chenying/Desktop/mwj/temp/Falign/sample-data/reads/ara_reads_2.fq
/Users/chenying/Desktop/mwj/temp/Falign/sample-data/reads/ara_reads_3.fq
```
And then input `read_list.lst` to `falign`:
``` shell
falign -num_threads 48 ^GATC ref/ara_2_4.fa.gz read_list.lst > map.sam
```

* Directory. If you have many `FASTQ`s of reads in a directory, just type the name of the directory:
``` shell
falign -num_threads 48 ^GATC ref/ara_2_4.fa.gz reads > map.sam
```
Note that only read files are allowed in the `reads` directory. If you put other files in the `reads`, `falign` will complain.

### </a> Output format

`falign` supports `SAM` output format:
```shell
$ falign -outfmt sam
```
and `PAF` format:
```shell
$ falign -outfmt paf
```
In each output result (note that every alignment has for offsets: read start, read end, reference start, reference end), `falign` adds the following additional fields:
* `qS:i:` the nearest restiction enzyme site to the read start position
* `qE:i:` the nearest restiction enzyme site to the read end position
* `vS:i:` the nearest restriction enzyme site to the reference start position
* `vE:i:` the nearest restriction enzyme site to the reference end position
* `pi:f:` percentage of identity of the alignment
* `gs:i:` the global chain score of the alignment's candidate
* `hm:Z:` a homologous map of the fragment

# <a name="S-u4falign"></a> Use `u4falign` for manipulating alignment results

`u4falign` is used for manipulating output results of `falign`. It supports the following commands:

 * `sam2salsa2`      Transfer SAM results to pairwise contacts for SALSA2
 * `sam23ddna`       Transfer SAM results to pairwise contacts for 3DDNA
 * `paf2salsa2`      Transfer PAF results to pairwise contacts for SALSA2
 * `paf23ddna`       Transfer PAF results ot pairwise contacts for 3DDNA
 * `sam2frag-sam`    Transfer SAM results output by falign to fragment SAM mapping results

# <a name="S-sam-diff"></a> Difference between `read-SAM` and `frag-SAM`

By convention, in `SAM` mapping results, a read may contain multiple mapping results. For saving space, the read sequence is usually only presented in the first mapping result of this read. In the second to last mapping results of this read, the `sequence field` is filled by a start `*`. `falign` outputs `SAM` mapping results in this manner. And we call `SAM` mapping results output in this manner `read-SAM`.  Since a Pore-C read always contains multiple fragments, a Pore-C read usually contains many mapping results:
```shell
0cd79600-51cf-4255-a6f7-0e9660721e85	16	2	1794202	1 ...
0cd79600-51cf-4255-a6f7-0e9660721e85	16	2	1791136	60 ...
0cd79600-51cf-4255-a6f7-0e9660721e85	0	2	1733063	60 ...
0cd79600-51cf-4255-a6f7-0e9660721e85	16	4	12727850	60 ...
0cd79600-51cf-4255-a6f7-0e9660721e85	0	2	1713207	60 ...
```
In the example above, the read `0cd79600-51cf-4255-a6f7-0e9660721e85` contains five alignment results.

Some tools such as whatshap taking `BAM` format as input will complain if there exist too many duplicate read names. In this case we suggest transfer `read-SAM` to `frag-SAM`. In `frag-SAM` every fragment is treadted as an individual read. We can use the `sam2frag-sam` command in the `u4falign` to transfer `read-SAM` to `frag-SAM`:
```shell
$ u4falign sam2frag-sam map.sam frag-map.sam
```
After transformation, we have
```shell
0cd79600-51cf-4255-a6f7-0e9660721e85_0000000001:000:0000000000:0001794201	16	2	1794202	1 ...
0cd79600-51cf-4255-a6f7-0e9660721e85_0000000001:001:0000000000:0001791135	16	2	1791136	60 ...
0cd79600-51cf-4255-a6f7-0e9660721e85_0000000001:002:0000000000:0001733062	0	2	1733063	60 ...
0cd79600-51cf-4255-a6f7-0e9660721e85_0000000001:003:0000000001:0012727849	16	4	12727850	60 ...
0cd79600-51cf-4255-a6f7-0e9660721e85_0000000001:004:0000000000:0001713206	0	2	1713207	60 ...
```
In `frag-SAM`, the `sequence field` in every alignment results is represented by the fragment sequence (not the whose read sequence). Note that we add a suffix to every read name to avoid duplicate read names in the `SAM` file. The meanings of the  fields in the suffix string is
```
read-id:fragment-id:reference-sequence-id:reference-mapping-position
```


# <a name="S-maintainers"></a>Maintainers

* Chuan-Le Xiao, xiaochuanle@126.com
* Ying Chen, chenying2016@gmail.com

# <a name="S-citation"></a>Citation

# <a name="S-license"></a>License

GPLv3
