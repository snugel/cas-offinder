# WARNING: Cas-OFFinder 3 is not production ready yet, it is known that the result can be different from that of Cas-OFFinder 2. For production use please use the latest Cas-OFFinder 2 instead.

Cas-OFFinder
==================================

Cas-OFFinder is OpenCL based, ultrafast and versatile program that searches for
potential off-target sites of CRISPR/Cas-derived RNA-guided endonucleases
(RGEN).

Cas-OFFinder is not limited by the number of mismatches and allows variations in
protospacer-adjacent motif (PAM) sequences recognized by Cas9, the essential
protein component in RGENs.

Requires an OpenCL device to run properly.

Cas-OFFinder is distributed under new BSD license (3-clauses).

Cas-OFFinder has been tested on the platforms below:
- Microsoft Windows (7 and 8)
- GNU/Linux (CentOS, OpenSUSE, Debian, Ubuntu/Elementary OS)
- Mac OS X (Mavericks)

Download
--------

Cas-OFFinder binaries are available at: https://github.com/snugel/cas-offinder/releases

For further information about installation, check out below `Installation` section.

CRISPR/Cas-derived RNA-guided endonucleases (RGEN)
--------------------------------------------------

RGENs use complementary base pairing to recognize target sites.

RGENs consist of two parts.
1. Guide RNA, as:
  - Dual RNA components comprising sequence-invariant tracrRNA and sequence-variable guide RNA termed crRNA, or,
  - Single-chain guide RNA (sgRNA) constructed by linking essential portions of tracrRNA and crRNA
2. Cas9 Protein
  - A fixed protein component that recognizes the protospacer adjacent motif (PAM) downstream of target
    DNA sequences corresponding to guide RNA.

PAM sites:
* __SpCas9__ from *Streptococcus pyogenes*: 5’-NGG-3’ (to a lesser extent, 5’-NAG-3’)
* __StCas9__ from *Streptococcus thermophilus*: 5’-NNAGAAW-3’ (W = A or T)
* __NmCas9__ from *Neisseria meningitidis*:5’-NNNNGMTT-3’ (M = A or C)
* __SaCas9__ from *Staphylococcus aureus*: 5’-NNGRRT-3’ (R = A or G)

Usage
-----

Cas-OFFinder can run with:

    cas-offinder {input_filename|-} {G|C|A}[device_id(s)] {output_filename|-}

G stands for using GPU devices, C for using CPUs, and A for using accelerators.

Optionally, you can set device IDs in addition to `G`, `C`, `A` to limit
number of devices used by Cas-OFFinder.  A range of device IDs can be
specified by using comma or colon to separate the IDs.  See examples
below.

The special filename `-` may be used in place of the input and output filename
to read and write from `stdin` and `stdout`, respectively.

A short example may be helpful!

First, download any target organism's chromosome FASTA files. You can find one
in below links:

- http://hgdownload.soe.ucsc.edu/downloads.html (UCSC genome sequences library)
- http://ensembl.org/info/data/ftp/index.html (Ensembl sequence library)

Extract all FASTA files in a directory.

For example (human chromosomes, in POSIX environment):

    $ wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
    $ mkdir -p /var/chromosome/human_hg19
    $ tar zxf chromFa.tar.gz -C /var/chromosome/human_hg19
    $ ls -al /var/chromosome/human_hg19
      drwxrwxr-x.  2 user group      4096 2013-10-18 11:49 .
      drwxrwxr-x. 16 user group      4096 2013-11-12 12:44 ..
      -rw-rw-r--.  1 user group 254235640 2009-03-21 00:58 chr1.fa
      -rw-rw-r--.  1 user group 138245449 2009-03-21 01:00 chr10.fa
      -rw-rw-r--.  1 user group 137706654 2009-03-21 01:00 chr11.fa
      -rw-rw-r--.  1 user group 136528940 2009-03-21 01:01 chr12.fa
      -rw-rw-r--.  1 user group 117473283 2009-03-21 01:01 chr13.fa
      -rw-rw-r--.  1 user group 109496538 2009-03-21 01:01 chr14.fa
      ...

Now, download Cas-OFFinder binary here,

https://sourceforge.net/projects/cas-offinder/files/Binaries

and save it to any directory you want.

And just try running it for a short help:

    $ ./cas-offinder
    Cas-OFFinder 3.0.0 beta (Jan 24 2021)

    Copyright (c) 2021 Jeongbin Park and Sangsu Bae
    Website: http://github.com/snugel/cas-offinder

    Usage: cas-offinder {input_filename|-} {C|G|A}[device_id(s)] {output_filename|-}
    (C: using CPUs, G: using GPUs, A: using accelerators)

    Example input file (DNA bulge size 2, RNA bulge size 1):
    /var/chromosomes/human_hg38
    NNNNNNNNNNNNNNNNNNNNNRG 2 1
    GGCCGACCTGTCGCTGACGCNNN 5
    CGCCAGCGTCAGCGACAGGTNNN 5
    ACGGCGCCAGCGTCAGCGACNNN 5
    GTCGCTGACGCTGGCGCCGTNNN 5

    Available device list:
    Type: GPU, ID: 0, <GeForce GTX 980> on <NVIDIA CUDA>
    Type: GPU, ID: 1, <GeForce GTX 980> on <NVIDIA CUDA>

Also it provides a list of all available OpenCL devices!

On Windows, if you encountered a missing .dll error, you may need to download and install [Visual C++
Redistributable Packages for Visual Studio 2015, 2017 and 2019](https://support.microsoft.com/en-us/topic/the-latest-supported-visual-c-downloads-2647da03-1eea-4433-9aff-95f26a218cc0).

Now you should create an input file:

- The first line is a directory containing FASTA or 2BIT files,
- The second line indicates the desired pattern including PAM site and optional DNA or RNA bulge sizes, separated by spaces,
- The remaining lines are the query sequences and maximum mismatch numbers, separated by spaces,
- Optionally, you can specify an ID for each query sequence, which will be included in the output (`Id` column).

The length of the desired pattern and the query sequences should be the same!

For the pattern and the query sequences, mixed bases are allowed to account for the degeneracy in PAM sequences.

Also, the number of mismatched bases is not limited!

Following codes are supported:

|   A   |    C   |   G   |   T   |
|:-----:|:------:|:-----:|:-----:|
|Adenine|Cytosine|Guanine|Thymine|

|   R  |   Y  |   S  |   W  |   K  |   M  |
|:----:|:----:|:----:|:----:|:----:|:----:|
|A or G|C or T|G or C|A or T|G or T|A or C|

|     B     |     D     |     H     |     V     |   N    |
|:---------:|:---------:|:---------:|:---------:|:------:|
|C or G or T|A or G or T|A or C or T|A or C or G|any base|

An example of an input file:

    /var/chromosomes/human_hg38
    NNNNNNNNNNNNNNNNNNNNNRG 2 1
    GGCCGACCTGTCGCTGACGCNNN 3 Seq1
    CGCCAGCGTCAGCGACAGGTNNN 3 Seq2
    ...

Save it as `input.txt`.

Now you can run Cas-OFFinder as following (using GPUs):

    $ ./cas-offinder input.txt G output.txt

Optionally, you can set the ID of devices to select a specific device
used by Cas-OFFinder. For example, to use the second GPU device
(0-based):

    $ ./cas-offinder input.txt G1 output.txt

You can use commas or colons to select range of devices:

    $ ./cas-offinder input.txt G0,1 output.txt

or

    $ ./cas-offinder input.txt G0:2 output.txt

Then the `TAB` separated output file will be generated with the following columns:

- The first column is the sequence id, by default line numbers (0-based), or given IDs from the input file,
- The second column is the type of bulge, the value is one of `X`, `DNA`, or `RNA`,
- The third column is the given query sequence (including gaps in case of DNA bulge),
- The fourth column is the actual sequence at the location (including gaps in case of RNA bulge; mismatched bases are noted in lowercase letters),
- The fifth column is the sequence name (if you downloaded it from UCSC or Ensembl, it is usually a chromosome name),
- The sixth column is the 0-based location of the off-target site (same convention as [Bowtie](https://github.com/BenLangmead/bowtie), not 1-based as IGV Viewer and others),
- The seventh column is the direction of sequence, either forward (+) or reverse (-) strand of the found sequence,
- The eighth column is the number of the mismatched bases,
- The last column is the size of bulge.

`output.txt` will be (shown here as a table for convenience):

|Id    |Bulge type|crRNA                      |DNA                        |Chromosome|Location|Direction|Mismatches|Bulge Size|
|:----:|:---------|:-------------------------:|:--------------------------|:--------:|:-------|:-------:|:---------|:--------:|
|Seq2  |DNA       |`CGCCAGCGTCAGCGACAGG--TNNN`|`CcCCAGtGTCAGCcACAGGGCTCAG`|chr1      |34978273|-        |3         |2         |
|Seq2  |DNA       |`CGC--CAGCGTCAGCGACAGGTNNN`|`CaCTCCAGCcTCAGCGACAGGcAAG`|chr1      |18173251|+        |3         |2         |
|Seq2  |DNA       |`CG--CCAGCGTCAGCGACAGGTNNN`|`CaCTCCAGCcTCAGCGACAGGcAAG`|chr1      |18173251|+        |3         |2         |
|Seq2  |DNA       |`C--GCCAGCGTCAGCGACAGGTNNN`|`CACtCCAGCcTCAGCGACAGGcAAG`|chr1      |18173251|+        |3         |2         |
|Seq1  |DNA       |`GGCCGACC--TGTCGCTGACGCNNN`|`GGCCcAgCTCTGTCGCTGACGgGAG`|chr1      |40979785|+        |3         |2         |
|Seq1  |DNA       |`GGCCGACCTGTCGCTGA--CGCNNN`|`GGCCGtCCTGTtGCTGAGACtCGGG`|chr1      |17408102|-        |3         |2         |
|Seq1  |DNA       |`GGCCGACCTGTCGCTG--ACGCNNN`|`GGCCGtCCTGTtGCTGAGACtCGGG`|chr1      |17408102|-        |3         |2         |
|Seq1  |DNA       |`GGCCGACCTGTCGCT--GACGCNNN`|`GGCCGtCCTGTtGCTGAGACtCGGG`|chr1      |17408102|-        |3         |2         |
|...   |
|Seq2  |X         |`CGCCAGCGTCAGCGACAGGTNNN`  |`CtCCAGCcTCAGCGACAGGcAAG`  |chr1      |18173253|+        |3         |0         |
|Seq2  |RNA       |`CGNCCCAGCGTCAGCGACAGGTNNN`|`C-CCCCAGtGTCActGACAGGTGGG`|chr1      |5502832 |-        |3         |1         |
|Seq2  |RNA       |`CGNCCCAGCGTCAGCGACAGGTNNN`|`t-CtCCAGCcTCAGCGACAGGcAAG`|chr1      |18173254|+        |3         |1         |
|Seq2  |RNA       |`CGCCGCAGCGTCAGCGACAGGTNNN`|`CG-CGCAGCGaCAGgGAgAGGTGAG`|chr1      |1273663 |-        |3         |1         |
|Seq2  |RNA       |`CGCCGCAGCGTCAGCGACAGGTNNN`|`CGC-GCAGCGaCAGgGAgAGGTGAG`|chr1      |1273663 |-        |3         |1         |
|Seq2  |RNA       |`CGCCAGCGTCGCAGCGACAGGTNNN`|`CGCCtGCG-CGgAGCtACAGGTGAG`|chr1      |18888560|+        |3         |1         |
|Seq2  |RNA       |`CGCCAGCGTCGCAGCGACAGGTNNN`|`aGCCAGCt-CtCAGCGACAGcTGAG`|chr1      |91631879|+        |3         |1         |
|Seq2  |RNA       |`CGCCAGCGTCGTAGCGACAGGTNNN`|`CGCCtGCGg-GgAGCtACAGGTGAG`|chr1      |18888560|+        |3         |1         |
|Seq2  |RNA       |`CGCCAGCGTCAGCGGCACAGGTNNN`|`CcCCAGaGTCAGC-GCACAGaTGGG`|chr1      |4006295 |-        |3         |1         |
|Seq2  |RNA       |`CGCCAGCGTCAGCGACGAAGGTNNN`|`tGCCAGCGgCAGCGA-GAAGtTTAG`|chr1      |8462729 |+        |3         |1         |
|...   |

Advanced Usage
--------------
Cas-OFFinder is mainly designed for CRISPR/Cas9 derived RGENs, however, it can
also be used for searching off-targets of other nucleases, e.g.
TALENs(Transcription activator-like effector nucleases) or ZFNs(Zinc-finger
nucleases), by specifying pattern sequence as all 'N's.

Example input file for TALENs:

    /var/chromosomes/human_hg38
    NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    TTCTGGAGGTGCCTGAGGCCNNNNNNNNNNNNGAGGCCACCTTTCCAGTCCA 5
    TGGCCAATGTGACGCTGACGNNNNNNNNNNNNCTGGAGACTCCAGACTTCCA 5
    ....

Installation
------------

* Cas-OFFinder requires an OpenCL-enabled device and corresponding runtime API
    pre-installed to run properly.

  OpenCL is supported on various platforms, including many Intel/AMD CPUs and
  NVidia/AMD graphic cards!  Before installing Cas-OFFinder, please check
  whether your device is an OpenCL-supported one.

  Khronos group provides an extensive list of supported devices here:

  http://www.khronos.org/conformance/adopters/conformant-products/#opencl

  Cas-OFFinder usually runs faster on GPUs than CPUs.
  If you want to purchase a new graphic card for fast analyzing speed,
  please check GPU benchmark results below as a reference:

  - https://compubench.com/result.jsp
  - http://www.luxrender.net/luxmark/top/top20/Room/GPU/1

  Recently, the OpenCL runtime binaries are already shipped with the device
  drivers in many cases - so you don't have to install anything to run
  Cas-OFFinder.

  But if it wasn't, you should download and install a proper OpenCL SDK to
  install runtime APIs.  In that case, download an OpenCL SDK among the links
  below.  If you know your device's vendor name, it is enough to install only
  your vendor's one.

  - AMD: http://developer.amd.com/tools-and-sdks/heterogeneous-computing/amd-accelerated-parallel-processing-app-sdk/downloads/
  - Intel: http://software.intel.com/en-us/vcsource/tools/opencl-sdk
  - NVidia: https://developer.nvidia.com/cuda-downloads

* Download Cas-OFFinder binary here,

  https://github.com/snugel/cas-offinder/releases

  or compile it from its source code (below section).

Compile
-------
  The OpenCL library is required to compile Cas-OFFinder.

  The OpenMP library is optional, for faster processing of 2BIT files.

  To support cross-platform compilation on various operating systems, the CMake
  build system is used. See http://www.cmake.org for more information.

  If you use Ubuntu linux, you can also install it via apt-get, `apt-get install cmake`.

  If you use macOS and homebrew, you can also install it via brew, `brew install cmake`.

  Otherwise, download CMake [here](http://www.cmake.org/cmake/resources/software.html).

  Checkout the source code of Cas-OFFinder with Git client,
  or download it manually on github website.


  In POSIX environment (`g++` should be pre-installed), launch terminal and type
  the following to build Cas-OFFinder:

      cmake -G "Unix Makefiles" .
      make

  On Windows (Visual Studio should be pre-installed), launch 'Visual Studio
  Command Prompt'.  (You can find it under 'Start menu' - 'Microsoft Visual
  Studio xxxx' - 'Visual Studio Tools'.)

  Assuming the CMake binary is installed in `C:\Program Files (x86)\CMake 3.19\bin`, type the following:

      "C:\Program Files (x86)\CMake 3.19\bin\cmake.exe" -G "NMake Makefiles" .
      nmake

  Then cas-offinder binary will be generated. Copy it wherever you want.

  If you have difficulty compiling cas-offinder, cmake can output system
  information to a file.  Please open a GitHub issue and attach the file.

      $ cmake --system-information cmake-sys-info.log

Module reference
----------------

* For reading/parsing FASTA files, the kseq.h library (developed by Heng Li) is used.

  The library is distributed under MIT licence.

  More information at:
  http://lh3lh3.users.sourceforge.net/parsefastq.shtml

* For supporting Dirent API on Windows environment,
  Dirent API for Microsoft Visual Studio is used.

  More information at:
  http://softagalleria.net/dirent.php

Download & Source
--------
The binaries can be downloaded from

https://sourceforge.net/projects/cas-offinder/files/Binaries

And the source code is distributed from

https://github.com/snugel/cas-offinder

Publication
-------

Cas-OFFinder is discussed in a
[paper](https://academic.oup.com/bioinformatics/article/30/10/1473/267560) published in the journal
_Bioinformatics_, Volume 30, Issue 10, 15 May 2014, Pages 1473–1475,
https://doi.org/10.1093/bioinformatics/btu048

Changelog
-------
* 3.0.0
  - Native support of DNA/RNA bulges
  - Update output format to display bulge information
  - Code cleanup / Better library handling by CMake (@richardkmichael)
* 2.4.1
  - Corrected warnings, code cleanups, and document updates (@richardkmichael)
  - Follow symlinks (@richardkmichael)
  - Allow stdin for input (@richardkmichael)
  - Tagging of input sequences (#28)
  - Github Actions based CI & automatic release
* 2.4
  - Corrected critical bug (The last 3 bases of 2bit input could be wrong)
  - Corrected bug (Segmentation fault if the match occurs at the very first location in chromosome)
  - Corrected bug (Cas-OFFinder does not follow symbolic links)
  - Now user can limit number of devices used by Cas-OFFinder.
  - Now Cas-OFFinder reports the name of platform.
  - Now user can set '-' as output file, then the output will be redirected to stdout. All other messages from Cas-OFFinder will be directed to stderr.
* 2.3
  - Removed cl.hpp due to lack of C++ binding support in the new OpenCL 2.0 standard.
  - Constant arguments are stored in constant or local memory, rather than global memory.
  - Added support for 2bit format.
  - Removed kseq.h
  - Precise running time measurment on POSIX platform.
* 2.2
  - Corrected a critical bug (when cas-offinder finds no binding sites in the given genome chunk, it crashes).
  - Now Cas-OFFinder reads whole fasta file at once, in order to achieve faster searching speed when it searches in FASTA files which contain many small scaffolds.
* 2.1
  - Using atomic operation, reduced computing load on CPU. In our benchmark, the total computation speed increased about twice as fast as before.
  - When lowercase sequences are given, convert them uppercase sequences before computation.
  - Corrected a bug (mixed bases were shown as lowercases letters, even they had been matched with normal bases).
  - Now supports 'accelerators', with 'A' option.
* 1.1
  - When Cas-OFFinder is launched without parameters, now it display available device list.
  - If the given chromosomes directory does not exist, now it returns an error message.
  - Corrected a bug (when Cas-OFFinder couldn't find any OpenCL device it would hang).
* 1.0
  - Initial release.

License
-------
Cas-OFFinder (except dirent.h) is licensed under the new BSD licence.

Copyright (c) 2021, Jeongbin Park and Sangsu Bae
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the name of the Seoul National University nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNERS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
