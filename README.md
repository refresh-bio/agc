Assembled Genomes Compressor
=
[![GitHub downloads](https://img.shields.io/github/downloads/refresh-bio/agc/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/refresh-bio/agc/releases)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/agc.svg?style=flag&label=Bioconda%20downloads)](https://anaconda.org/bioconda/agc)

Assembled Genomes Compressor (AGC) is a tool designed to compress collections of de-novo assembled genomes.
It can be used for various types of datasets: short genomes (viruses) as well as long (humans).

The tool offers high compression ratios, especially for high-quality genomes. 
For example the 96 haplotype sequences from the [Human Pangenome Project](https://github.com/human-pangenomics/HPP_Year1_Assemblies) (47 samples), [GRCh 38 reference](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz), and [CHM13 v.1.1 assembly](https://github.com/marbl/CHM13) containing about 290Gb are squeezed to less than 1.5GB. 
The compressed samples are easily accessible as *agc* offers extraction of single samples or contigs in a few seconds.
The compression is also fast. 
On a AMD TR 3990X-based machine (32 threads used) it takes about 12 minutes to compress the HPP collection.


## Quick start

```bash
git clone https://github.com/refresh-bio/agc
cd agc && make

# Compress a collection of 3 genomes
./agc create ref.fa in1.fa in2.fa > col.agc                         # file names given in command-line
./agc create ref.fa in1.fa.gz in2.fa.gz > col.agc                   # gzipped non-reference FASTA files
./agc create -i fn.txt ref.fa > col.agc                             # fl.txt contains 2 file names in seperate lines: 
                                                                    # in1.fa in2.fa
./agc create -a -i fn.txt ref.fa > col.agc                          # adaptive mode (use for bacterial data)
./agc create -i fn.txt -o col.agc ref.fa                            # output file name is specified as a parameter
./agc create -i fn.txt -o col.agc -k 29 -l 22 -b 100 -t 16 ref.fa   # same as above, but manual selection 
                                                                    # of compression parameters
./agc create -c -o col.agc ref.fa samples.fa                        # compress samples stored in a single file
                                                                    # (reference must be given separately)

# Add new genomes to the collection
./agc append in.agc in3.fa in4.fa > out.agc                         # add 2 genomes to the compressed archive
./agc append -i fn.txt in.agc -o out.agc                            # add genomes (fn.txt contains file names)
./agc append -a -i fn.txt in.agc -o out.agc                         # add genomes (adaptive mode)

# Extract all genomes from the compressed archive
./agc getcol in.agc > out.fa                                        # extract all samples
./agc getcol -o out_path/ in.agc                                    # extract all samples and store them in separate files

# Extract a genome or genomes from the compressed archive
./agc getset in.agc in1 > out.fa                                    # extract sample in1 from the archive
./agc getset in.agc in1 in2 > out.fa                                # extract samples in1 and in2 from the archive

# Extract contigs from the compressed archive
./agc getctg in.agc ctg1 ctg2 > out.fa                              # extract contigs ctg1 and ctg2 from the archive
./agc getctg in.agc ctg1@gn1 ctg2@gn2 > out.fa                      # extract contigs ctg1 from genome gn1 and ctg2 from gn2 
                                                                    # (useful if contig names are not unique)
./agc getctg in.agc ctg1@gn1:from1-to1 ctg2@gn2:from2-to2 > out.fa  # extract parts of contigs 
./agc getctg in.agc ctg1:from1-to1 ctg2:from2-to2 > out.fa          # extract parts of contigs 

# List genome names in the archive
./agc listset in.agc > out.txt                                      # list sample names

# List contig names in the archive
./agc listctg in.agc gn1 gn2 > out.txt                              # list contig names in genomes gn1 and gn2

# Show info about the compression archive
./agc info in.agc                                                   # show some stats, parameters, command-lines 
                                                                    # used to create and extend the archive

```

## Installation and configuration
agc should be downloaded from https://github.com/refresh-bio/agc and compiled. The supported OS are:
* Windows: Visual Studio 2022 solution provided,
* Linux: make project (G++ 9.0 or newer required),
* MacOS: make project (G++ 9.0 or newer required).

The release contains a set of [precompiled binaries](https://github.com/refresh-bio/agc/releases) for Windows, Linux, and OS X. 

The software is also available on [Bioconda](https://anaconda.org/bioconda/agc):
```
conda install -c bioconda agc
```
For detailed instructions on how to set up Bioconda, please refer to the [Bioconda manual](https://bioconda.github.io/user/install.html#install-conda).


## Version history
* 3.0 (22 Dec 2022)
  * Improved compression (slightly better ratio).
  * Improved archive format &mdash; much faster queries for archives containing large number of samples.
  * Bugfixes.
* 2.1 (9 May 2022)
  * Bugfix in append mode. (In version 2.0, running append could produce improper archive.)
* 2.0 (5 Apr 2022)
  * Optional adaptive mode (especially for bacterial data).
  * New mode: decompression of whole collection.
  * New archive format (a bit more compact): AGC 1.x tool cannot read AGC 2 archives, but AGC 2.x tool can operate on AGC 1.x and AGC 2.x archives.
* 1.1 (14 Jan 2022)
  * Small bugfixes.
* 1.0 (23 Dec 2021)
  * First public release.

## Usage

`agc <command> [options]`

Command:
* `create`   - create archive from FASTA files
* `append`   - add FASTA files to existing archive
* `getcol`   - extract all samples from archive
* `getset`   - extract sample from archive
* `getctg`   - extract contig from archive
* `listset`  - list sample names in archive
* `listctg`  - list sample and contig names in archive
* `info`     - show some statistics of the compressed data

### Creating new archive

`agc create [options] <ref.fa> [<in1.fa> ...] > <out.agc>`

Options:
* `-a`             - adaptive mode (default: false)
* `-b <int>`       - batch size (default: 50; min: 1; max: 1000000000)
* `-c`             - concatenated genomes in a single file (default: false)
* `-d`             - do not store cmd-line (default: false)
* `-i <file_name>` - file with FASTA file names (alterantive to listing file names explicitely in command line)
* `-k <int>`       - k-mer length (default: 31; min: 17; max: 32)
* `-l <int>`       - min. match length (default: 20; min: 15; max: 32)
* `-o <file_name>` - output to file (default: output is sent to stdout)
* `-s <int>`       - expected segment size (default: 60000; min: 100; max: 1000000)
* `-t <int>`       - no. of threads (default: no. logical cores / 2; min: 1; max: no. logical. cores)
* `-v <int>`       - verbosity level (default: 0; min: 0; max: 2)

#### Hints
FASTA files can be optionally gzipped. It is, however, recommended (for performance reasons) to use uncompressed reference FASTA file.

If all samples are given in a single file (concatenated genomes mode) the reference must be given in a separate file.

Setting parameters allows difference compromises, usually between compressed size and decompression time. The impact of the most important options is discussed below.
* *Batch size* specifies the internal granularity of the data. If it is set to low values then in the extraction mode *agc* needs to internally decompress just a small fraction of the archive, which is fast. From the opposite side, if batch size is huge then *agc* needs to internally decompress the complete archive. This is, however, just partial decompression (to some compacted form), so the time can still be acceptable. You can experiment with this parameter as set it to your needs. Note, however, that the parameter can be set only in the `create` mode, so it cannot be changed later, when the archive is extended.
* *k-mer length* is an internal parameter which specifies the length of _k-mers_ used to split the genomes into shorter parts (segments) for compression. The parameter should not be changed without a necessity. In fact setting it to a value between 25 and 32 should not change too much in the compression ratio and (de)compression speeds. Nevertheless, setting it to a too low value (e.g., 17 for human genomes) will make the compression much harder and you should expect poor results.
* *minimal match length* is an internal parameter specifying the minimal match length when the similarities between contigs are looked for. If you really want, you can try to change it. Nevertheless, the impact on the compression ratios and (de)compression speeds should be insignificant.
* *segment size* specifies how the contigs are splitted into shorter fragments (segments) during the compresssion. This is an expected segment size and some segments can be much longer. In general, the more similar the genomes in a collection the larger the parameter can be. Nevertheless, the impact of its value on the compression ratios and (de)compression speeds is limited. If you want, you can experiment with it. Note that for short sequences, especially for virues, the segment size should be smaller, you can try 10000 or similar values.
* *no. of threads* impacts the running time. For large genomes (e.g., human) the parallelization of the compression is realatively good and you can use 30 or more threads. Setting *segment size* to larger values can improve paralelization a bit.
* *adaptive mode* allows to look for new splitters in all genomes (not only reference). It needs more memory but give significant gains in compression ratio and speed especially for highly divergent genomes, e.g., bacterial.



### Append new genomes to the existing archive

`agc append [options] <in.agc> [<in1.fa> ...] > <out.agc>`

Options:
* `-c`             - concatenated genomes in a single file (default: false)
* `-d`             - do not store cmd-line (default: false)
* `-i <file_name>` - file with FASTA file names (alterantive to listing file names explicitely in command line)
* `-o <file_name>` - output to file (default: output is sent to stdout)
* `-t <int>`       - no. of threads (default: no. logical cores / 2; min: 1; max: no. logical. cores)
* `-v <int>`       - verbosity level (default: 0; min: 0; max: 2)

#### Hints
FASTA files can be optionally gzipped.

### Decompress whole collection
`agc getcol [options] <in.agc> > <out.fa>`

Options:\n";
* `-l <int>`         - line length (default: 80; min: 40; max: 2000000000)
* `-o <output_path>` - output to files at path (default: output is sent to stdout)
* `-t <int>`         - no. of threads (default: no. logical cores / 2; min: 1; max: no. logical. cores)
* `-v <int>`         - verbosity level (default: 0; min: 0; max: 2)

#### Hints
If output path is specified then it must be an existing directory.
Each sample will be stored in a separate file (the files in the directory will be overwritten if their names are the same as sample name).

### Extract genomes from the archive

`agc getset [options] <in.agc> <sample_name1> [<sample_name2> ...] > <out.fa>`

Options:
* `-l <int>`       - line length (default: 80; min: 40; max: 2000000000)
* `-o <file_name>` - output to file (default: output is sent to stdout)
* `-t <int>`       - no. of threads (default: no. logical cores / 2; min: 1; max: no. logical. cores)
* `-p`             - disable file prefetching (useful for short genomes)
* `-v <int>`       - verbosity level (default: 0; min: 0; max: 2)
  
#### Hints
If output file name ends with `.gz' the output file will be gzipped.
  
### Extract contigs from the archive

`agc getctg [options] <in.agc> <contig1> [<contig2> ...] > <out.fa>` <br />
`agc getctg [options] <in.agc> <contig1@sample1> [<contig2@sample2> ...] > <out.fa>` 
`agc getctg [options] <in.agc> <contig1:from-to>[<contig2:from-to> ...] > <out.fa>`
`agc getctg [options] <in.agc> <contig1@sample1:from1-to1> [<contig2@sample2:from2-to2> ...] > <out.fa>`

Options:
* `-l <int>`       - line length (default: 80; min: 40; max: 2000000000)
* `-o <file_name>` - output to file (default: output is sent to stdout)
* `-t <int>`       - no. of threads (default: no. logical cores / 2; min: 1; max: no. logical. cores)
* `-p`             - disable file prefetching (useful for short queries)
* `-v <int>`       - verbosity level (default: 0; min: 0; max: 2)

#### Hints
If output file name ends with `.gz' the output file will be gzipped.

  ### List samples in the archive
`agc listset [options] <in.agc> > <out.txt>`

Options:
* `-o <file_name>` - output to file (default: output is sent to stdout)
  
### List contigs in the archive

`agc listctg [options] <in.agc> <sample1> [<sample2> ...] > <out.txt>`

Options:
* `-o <file_name>` - output to file (default: output is sent to stdout)
  
### Show some info about the archive

`agc info [options] <in.agc> > <out.txt>`

Options:
* `-o <file_name>` - output to file (default: output is sent to stdout)


## AGC decompression library
AGC files can be accessed also with C/C++ or Python library. 

### C/C++ libraries
The C and C++ APIs are provided in src/lib-cxx/agc-api.h file (in C++ you can use C or C++ API).
You can also take a look at src/examples to see both APIs in use.

### Python library
AGC files can be accessed also with Python wrapper for AGC API, which was created using pybind11, version 2.8.1. 
It is available for Linux and macOS.  

Warning: Python binding is experimental. The library used to create binding as well as public interface may change in the future.

Python module wrapping AGC API must be compiled. To compile it run:
 ```
 make py_agc_api
 ```

As a result of pybind11 *.so file is created and may be used as a python module. 
The following file is created: 
py_agc_api.`` `python3-config --extension-suffix` ``

To be able to use this file one should make it visible for Python. 
One way to do this is to extend PYTHONPATH environment variable. It can be done by running:
```
source py_agc_api/set_path.sh
```
The example usage of Python wrapper for AGC API is presented in file: `py_agc_api/py_agc_test.py`
To test it, run:
```
python3 py_agc_api/py_agc_test.py 
```
It reads the AGC file from the toy example (`toy_ex/toy_ex.agc`).

## Toy example
There are four example FASTA files (`ref.fa`, `a.fa`, `b.fa` and `c.fa`) in the `toy_ex` folder. They can be used to test AGC. The `toy_ex/toy_ex.agc` is an AGC archive created with:
```
agc create -o toy_ex/toy_ex.agc toy_ex/ref.fa toy_ex/a.fa toy_ex/b.fa toy_ex/c.fa
```
The AGC file is read in Python test script (`py_agc_api/py_agc_test.py`) and can be used in the example using C/C++ library (`src/examples`).

For more options see Usage section.

## Large datasets
Archives of 94 haplotype human assemblies <a href="https://github.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0">released by HPRC</a> in 2021 as well as 619,750 complete SARC-Cov-2 genomes <a href="https://www.ncbi.nlm.nih.gov/datasets/coronavirus/genomes/">published by NCBI</a> can be downloaded from <a href="https://zenodo.org/record/5826274">Zenodo</a>.
  
## Citing
S. Deorowicz, A. Danek, H. Li,
AGC: Compact representation of assembled genomes with fast queries and updates.
Bioinformatics, btad097 (2023)
https://doi.org/10.1093/bioinformatics/btad097
