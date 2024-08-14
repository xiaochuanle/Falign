# <a name="S-table-of-contents"></a> Table of Contents

* [Introduction](#S-introduction)
* [Tools](#S-tools)
* [Installation](#S-installation)
* [Usage](#S-usage)
* [Use `u4falign` for manipulating alignment results](#S-u4falign)
* [Difference between `read-SAM` and `frag-SAM`](#S-sam-diff)

# <a name="S-introduction"></a> Introduction

`Falign` is a sequence alignment toolkit for long noisy chromosome conformation capture (3C) reads, such as Pore-C.
`Falign` is written in C and C++ programming language.
`Falign` relies on the [`HTSLIB`](https://github.com/samtools/htslib) library. We have dowanloaded one copy of that and save in `./htslib/htslib-1.19.1.tar.bz2` to make it convenient for off-line installation.

# <a name="S-tools"></a> Tools

Two tools are released together in this toolkit.

* `falign`. The alignment tool.
* `dip3d`. A utility used internally by [`Dip3D`](https://github.com/xiaochuanle/dip3d), a diploid human 3D genome construction algorithm.

# <a name="S-installation"></a> Installation

### Step 1: Download source codes
```shell
$ git clone https://github.com/xiaochuanle/Falign.git
```

### Step 2: Install `HTSLIB`
```shell
$ cd Falign
$ ./htslib/install-htslib.sh
```

### Step 3: Export `HTSLIB` to an environment variable

In Step 2 above, we install `HTSLIB` in directory `./htslib/htslib`:
```shell
$ ls ./htslib/htslib
bin  include  lib  share
```
We have to tell `Falign` where to link `HTSLIB`. To do this, we export the location of this directory to an environment variable name `LIBHTS`:
```shell
$ export LIBHTS=$(pwd)/htslib/htslib
$ echo $LIBHTS
/share/home/chuanlex/chenying/data/Falign/htslib/htslib
```

### Step 4: Install `Falign`
***Before installing `Falign`, make sure you have executed Step 3 above.*** Otherwise the linker will complain that it cannot find `-lhts`.
```shell
$ cd src/ && make && cd ..
```

`Falign` is installed in the directory `Linux-amd64/bin/`:
```shell
$ ls Linux-amd64/bin/
dip3d  falign
```

# <a name="S-usage"></a> Usage

## Sample reads and reference

The directory `ara-sample-data/` provides sample reference file `reference.fa.gz` and sample reads file `reads.fq.gz`.

## Step 1: Build repeat regions for the reference genome
```shell
./Linux-amd64/bin/falign build-repeat ara-sample-data/reference.fa.gz ara-sample-data/reference.fa.gz.repeat.bed
```

## Step 2: Map reads to the reference
```shell
./Linux-amd64/bin/falign -repeat_bed ara-sample-data/reference.fa.gz.repeat.bed -num_threads 4 ara-sample-data/reference.fa.gz ara-sample-data/reads.fq.gz > map.paf
```

The mapping results are output to the file `map.paf` in PAF format.

## <a name="output-format"></a> Output format

`Falign` supports different output formats (specified by the `-outfmt` option):
* PAF
* SAM
* BAM
* FRAG-SAM
* FRAG-BAM

In each output result (note that every alignment has for offsets: read start, read end, reference start, reference end), `falign` adds the following additional fields:
* `qS:i:` the nearest restiction enzyme site to the read start position
* `qE:i:` the nearest restiction enzyme site to the read end position
* `vS:i:` the nearest restriction enzyme site to the reference start position
* `vE:i:` the nearest restriction enzyme site to the reference end position
* `pi:f:` percentage of identity of the alignment
* `SA:Z:` a homologous map of the fragment

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
$ u4falign sam2frag-sam map.sam > frag-map.sam
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
