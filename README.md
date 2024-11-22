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
git clone --recurse-submodules https://github.com/refresh-bio/agc
cd agc

# Linux compilation
make

# MacOS compilation: must specify g++ compiler
make CXX=g++-11

# Compress a collection of 3 genomes
bin/agc create ref.fa in1.fa in2.fa > col.agc                         # file names given in command-line
bin/agc create ref.fa in1.fa.gz in2.fa.gz > col.agc                   # gzipped non-reference FASTA files
bin/agc create -i fn.txt ref.fa > col.agc                             # fl.txt contains 2 file names in seperate lines: 
                                                                    # in1.fa in2.fa
bin/agc create -a -i fn.txt ref.fa > col.agc                          # adaptive mode (use for bacterial data)
bin/agc create -i fn.txt -o col.agc ref.fa                            # output file name is specified as a parameter
bin/agc create -i fn.txt -o col.agc -k 29 -l 22 -b 100 -t 16 ref.fa   # same as above, but manual selection 
                                                                    # of compression parameters
bin/agc create -c -o col.agc ref.fa samples.fa                        # compress samples stored in a single file
                                                                    # (reference must be given separately)

# Add new genomes to the collection
bin/agc append in.agc in3.fa in4.fa > out.agc                         # add 2 genomes to the compressed archive
bin/agc append -i fn.txt in.agc -o out.agc                            # add genomes (fn.txt contains file names)
bin/agc append -a -i fn.txt in.agc -o out.agc                         # add genomes (adaptive mode)

# Extract all genomes from the compressed archive
bin/agc getcol in.agc > out.fa                                        # extract all samples
bin/agc getcol -o out_path/ in.agc                                    # extract all samples and store them in separate files

# Extract a genome or genomes from the compressed archive
bin/agc getset in.agc in1 > out.fa                                    # extract sample in1 from the archive
bin/agc getset in.agc in1 in2 > out.fa                                # extract samples in1 and in2 from the archive

# Extract contigs from the compressed archive
bin/agc getctg in.agc ctg1 ctg2 > out.fa                              # extract contigs ctg1 and ctg2 from the archive
bin/agc getctg in.agc ctg1@gn1 ctg2@gn2 > out.fa                      # extract contigs ctg1 from genome gn1 and ctg2 from gn2 
                                                                    # (useful if contig names are not unique)
bin/agc getctg in.agc ctg1@gn1:from1-to1 ctg2@gn2:from2-to2 > out.fa  # extract parts of contigs 
bin/agc getctg in.agc ctg1:from1-to1 ctg2:from2-to2 > out.fa          # extract parts of contigs 

# List genome names in the archive
bin/agc listset in.agc > out.txt                                      # list sample names

# List contig names in the archive
bin/agc listctg in.agc gn1 gn2 > out.txt                              # list contig names in genomes gn1 and gn2

# Show info about the compression archive
bin/agc info in.agc                                                   # show some stats, parameters, command-lines 
                                                                    # used to create and extend the archive

```

## Installation and configuration
agc should be downloaded from https://github.com/refresh-bio/agc and compiled. The supported OS are:
* Windows: Visual Studio 2022 solution provided,
* Linux: make project (G++ 10.x or newer required),
* MacOS: make project (G++ 11.x, 12.x, or 13.x required; GNUMake 4.3 or newer required).

### Compilation options
For better performance gzipped input is readed using [isa-l](https://github.com/intel/isa-l) library for x64 CPUs. 
This, however, requires [NASM](https://github.com/netwide-assembler/nasm) compiler to be installed (you can install it from GitHub or are `nasm` package, e.g., `sudo apt install nasm`).
If NASM is not present (or at ARM-based CPUs), the [zlib-ng](https://github.com/zlib-ng/zlib-ng) is used.

Compilation with default options optimizes the tool for the native platform. 
If you want more control, you can specify the platform:
```
make PLATFORM=arm8    # compilation for ARM-based machines (turns on `-march=armv8-a`)
make PLATFORM=m1      # compilation for M1/M2/... (turns on `-march=armv8.4-a`)
make PLATFORM=sse2    # compilation for x64 CPUs with SSE2 support
make PLATFORM=avx     # compilation for x64 CPUs with AVX support
make PLATFORM=avx2    # compilation for x64 CPUs with AVX2 support
```

You can also specify the g++ compiler version (if installed):
```
make CXX=g++-11
make CXX=g++-12
```

### Prebuild releases
The release contains a set of [precompiled binaries](https://github.com/refresh-bio/agc/releases) for Windows, Linux, and OS X. 

The software is also available on [Bioconda](https://anaconda.org/bioconda/agc):
```
conda install -c bioconda agc
```
For detailed instructions on how to set up Bioconda, please refer to the [Bioconda manual](https://bioconda.github.io/user/install.html#install-conda).


## Version history
* 3.2 (21 Nov 2024)
  * Improved compression speed.
  * Optional fallback procedure to improve compression ratio.
  * Streaming mode for decompression &ndash; slower, but less memory needed.
  * Binaries (agc command-line tool and libraries) are now in the `bin` directory.
  * Small bug fixes.
* 3.1 (18 Mar 2024)
  * Improved compression speed for gzipped input.
  * Support for ARM-based CPUs (e.g., Mac M1/M2/...).
  * Reporting reference sample name.
  * Fixed truncating .fa from gzipped input.
  * Fixed Python lib GetCtgSeq().
  * Added optional gzipping in decompression modes.
  * Small bug fixes.
* 3.0 (22 Dec 2022)
  * Improved compression (slightly better ratio).
  * Improved archive format &mdash; much faster queries for archives containing a large number of samples.
  * Bugfixes.
* 2.1 (9 May 2022)
  * Bugfix in append mode. (In version 2.0, running append could produce an improper archive.)
* 2.0 (5 Apr 2022)
  * Optional adaptive mode (especially for bacterial data).
  * New mode: decompression of the whole collection.
  * New archive format (a bit more compact): AGC 1.x tool cannot read AGC 2 archives, but AGC 2.x tool can operate on AGC 1.x and AGC 2.x archives.
* 1.1 (14 Jan 2022)
  * Small bug fixes.
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
* `listref`  - list reference sample name in archive
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
* `-f <float>`     - fraction of fall-back minimizers (default: 0.000000; min: 0.000000; max: 0.050000)
* `-i <file_name>` - file with FASTA file names (alternative to listing file names explicitly in command line)
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
* *k-mer length* is an internal parameter that specifies the length of _k-mers_ used to split the genomes into shorter parts (segments) for compression. The parameter should not be changed without a necessity. In fact setting it to a value between 25 and 32 should not change too much in the compression ratio and (de)compression speeds. Nevertheless, setting it to a too-low value (e.g., 17 for human genomes) will make the compression much harder, and you should expect poor results.
* *minimal match length* is an internal parameter specifying the minimal match length when the similarities between contigs are looked for. If you really want, you can try to change it. Nevertheless, the impact on the compression ratios and (de)compression speeds should be insignificant.
* *segment size* specifies how the contigs are split into shorter fragments (segments) during the compression. This is an expected segment size and some segments can be much longer. In general, the more similar the genomes in a collection the larger the parameter can be. Nevertheless, the impact of its value on the compression ratios and (de)compression speeds is limited. If you want, you can experiment with it. Note that for short sequences, especially for virues, the segment size should be smaller, you can try 10000 or similar values.
* *no. of threads* impacts the running time. For large genomes (e.g., human) the parallelization of the compression is realatively good and you can use 30 or more threads. Setting *segment size* to larger values can improve paralelization a bit.
* *adaptive mode* allows to look for new splitters in all genomes (not only reference). It needs more memory but give significant gains in compression ratio and speed especially for highly divergent genomes, e.g., bacterial.
* *fall-back minimizers* allow to look for matching segment when it cannot be found using splitting <i>k</i>-mers. The parameter specifies what fraction of all <i>k</i>-mers will be used in the fall-back procedure. This can be useful for highly divergent genomes. For bacterial genomes, a value of 0.01 should be a reasonable choice. The improvement of compression ratio can be up to 20%. For human data, you can try using 0.001. The potential gain can be smaller like 2&ndash;3%. This slows down the compression. Use this feature with care, as sometimes it is better not to add a segment to a group if the splitters do not match and start a new group instead.


### Append new genomes to the existing archive

`agc append [options] <in.agc> [<in1.fa> ...] > <out.agc>`

Options:
* `-c`             - concatenated genomes in a single file (default: false)
* `-d`             - do not store cmd-line (default: false)
* `-f <float>`     - fraction of fall-back minimizers (default: 0.000000; min: 0.000000; max: 0.050000)
* `-i <file_name>` - file with FASTA file names (alternative to listing file names explicitly in command line)
* `-o <file_name>` - output to file (default: output is sent to stdout)
* `-t <int>`       - no. of threads (default: no. logical cores / 2; min: 1; max: no. logical. cores)
* `-v <int>`       - verbosity level (default: 0; min: 0; max: 2)

#### Hints
FASTA files can be optionally gzipped.

### Decompress whole collection
`agc getcol [options] <in.agc> > <out.fa>`

Options:\n";
* `-g <int>`         - optional gzip with given level (default: 0; min: 0; max: 9)
* `-l <int>`         - line length (default: 80; min: 40; max: 2000000000)
* `-o <output_path>` - output to files at path (default: output is sent to stdout)
* `-t <int>`         - no. of threads (default: no. logical cores / 2; min: 1; max: no. logical. cores)
* `-v <int>`         - verbosity level (default: 0; min: 0; max: 2)

#### Hints
If output path is specified then it must be an existing directory.
Each sample will be stored in a separate file (the files in the directory will be overwritten if their names are the same as sample name).
Samples can be gzipped when `-g` flag is provided.

### Extract genomes from the archive

`agc getset [options] <in.agc> <sample_name1> [<sample_name2> ...] > <out.fa>`

Options:
* `-g <int>`       - optional gzip with given level (default: 0; min: 0; max: 9)
* `-l <int>`       - line length (default: 80; min: 40; max: 2000000000)
* `-o <file_name>` - output to file (default: output is sent to stdout)
* `-s`             - enable streaming mode (slower but needs less memory)
* `-t <int>`       - no. of threads (default: no. logical cores / 2; min: 1; max: no. logical. cores)
* `-p`             - disable file prefetching (useful for short genomes)
* `-v <int>`       - verbosity level (default: 0; min: 0; max: 2)
  
#### Hints
Samples can be gzipped when `-g` flag is provided.
  
### Extract contigs from the archive

`agc getctg [options] <in.agc> <contig1> [<contig2> ...] > <out.fa>` <br />
`agc getctg [options] <in.agc> <contig1@sample1> [<contig2@sample2> ...] > <out.fa>` 
`agc getctg [options] <in.agc> <contig1:from-to>[<contig2:from-to> ...] > <out.fa>`
`agc getctg [options] <in.agc> <contig1@sample1:from1-to1> [<contig2@sample2:from2-to2> ...] > <out.fa>`

Options:
* `-g <int>`       - optional gzip with given level (default: 0; min: 0; max: 9)
* `-l <int>`       - line length (default: 80; min: 40; max: 2000000000)
* `-o <file_name>` - output to file (default: output is sent to stdout)
* `-s`             - enable streaming mode (slower but needs less memory)
* `-t <int>`       - no. of threads (default: no. logical cores / 2; min: 1; max: no. logical. cores)
* `-p`             - disable file prefetching (useful for short queries)
* `-v <int>`       - verbosity level (default: 0; min: 0; max: 2)

#### Hints
Contigs can be gzipped when `-g` flag is provided.

### List reference sample name in the archive
`agc listref [options] <in.agc> > <out.txt>`

Options:
* `-o <file_name>` - output to file (default: output is sent to stdout)
  
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
AGC files can be accessed also with Python wrapper for AGC API, which was created using pybind11, version 2.11.1. 
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

For more options, see the Usage section.

## Large datasets
Archives of 94 haplotype human assemblies <a href="https://github.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0">released by HPRC</a> in 2021 as well as 619,750 complete SARC-Cov-2 genomes <a href="https://www.ncbi.nlm.nih.gov/datasets/coronavirus/genomes/">published by NCBI</a> can be downloaded from <a href="https://zenodo.org/record/5826274">Zenodo</a>.
  
## Citing
S. Deorowicz, A. Danek, H. Li,
AGC: Compact representation of assembled genomes with fast queries and updates.
Bioinformatics, btad097 (2023)
https://doi.org/10.1093/bioinformatics/btad097
