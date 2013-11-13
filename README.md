Cas-OFFinder
==================================

Cas-OFFinder is OpenCL based, fast and dedicated tool for searching
CRISPR/Cas-derived RNA-guided endonucleases (RGEN) off-target sites.

Cas-OFFinder has no limitation on searching patterns and numbers of mismatched bases!

Requires any OpenCL library exist to compile and run.

Cas-OFFinder is distributed under new BSD license.

Usage
-------

Cas-OFFinder can run with:
  
    cas-offinder {input_file} {G|C} {output_file}

G stands for using all available GPU devices, and C stands for using all CPUs.

A short example may be helpful!

First, download any target organism's chromosome FASTA files. You can find one in below links:

- http://hgdownload.soe.ucsc.edu/downloads.html (UCSC genome sequences library)
- http://ensembl.org/info/data/ftp/index.html (Ensembl sequence library)

Extract all FASTA files in a directory.

For example (human chromosomes, in POSIX environment):
    
    $> wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
    $> mkdir -p /var/chromosome/human_hg19
    $> tar zxf chromFa.tar.gz -C /var/chromosome/human_hg19
    $> ls -al /var/chromosome/human_hg19
      drwxrwxr-x.  2 user group      4096 2013-10-18 11:49 .
      drwxrwxr-x. 16 user group      4096 2013-11-12 12:44 ..
      -rw-rw-r--.  1 user group 254235640 2009-03-21 00:58 chr1.fa
      -rw-rw-r--.  1 user group 138245449 2009-03-21 01:00 chr10.fa
      -rw-rw-r--.  1 user group 137706654 2009-03-21 01:00 chr11.fa
      -rw-rw-r--.  1 user group 136528940 2009-03-21 01:01 chr12.fa
      -rw-rw-r--.  1 user group 117473283 2009-03-21 01:01 chr13.fa
      -rw-rw-r--.  1 user group 109496538 2009-03-21 01:01 chr14.fa
      ...

Now, download Cas-OFFinder binary to any directory you want.

And just run it for a short help:

    $> ./cas-offinder
        Cas-OFFinder v1.0 (2013-11-03)
        
        Copyright 2013 Jeongbin Park and Sangsu Bae
        Website: http://github.com/snugel/cas-offinder
        
        Usage: cas-offinder {input_file} {C|G} {output_file}
        (C: using CPU, G: using GPU)

        ...
        
Now the input file should be created.

- The first line of the input file gives directory path containing chromosomes FASTA files,
- The second line indicates the PAM site,
- ...and following lines are the query sequences and maximum mistmatch numbers, seperated by spaces.
(The length of PAM site and the query sequences should be the same!)

For all seqeucnes, following codes (including wildcards) are supported:

   A   |    C   |   G   |   T   
:-----:|:------:|:-----:|:-----:
Adenine|Cytosine|Guanine|Thymine

   R  |   Y  |   S  |   W  |   K  |   M  
:----:|:----:|:----:|:----:|:----:|:----:
A or G|C or T|G or C|A or T|G or T|A or C

     B     |     D     |     H     |     V     |   N
:---------:|:---------:|:---------:|:---------:|:------:
C or G or T|A or G or T|A or C or T|A or C or G|any base

For example, write like this:

    /var/chromosomes/human_hg19
    NNNNNNNNNNNNNNNNNNNNNRG
    GGCCGACCTGTCGCTGACGCNNN 5
    CGCCAGCGTCAGCGACAGGTNNN 5
    ACGGCGCCAGCGTCAGCGACNNN 5
    GTCGCTGACGCTGGCGCCGTNNN 5
    ...

Save it as 'input.txt'.

Now you can run Cas-OFFinder as following:
 
    $> ./cas-offinder input.txt G out.txt
    ...
 
Then the output file will be generated :
- The first column of the output file indicates the given query sequence,
- The second column is the FASTA file name (usually chromosome name),
- The third column is the position of the off-target site,
- The forth column shows the actual sequence at the position
with indicating mismatched bases in lowercase letters,
- The fifth column is the direction of the found sequence (same convention with bowtie),
- ... and the last line is the number of the mismatched bases.

out.txt:

    GGCCGACCTGTCGCTGACGCNNN chr8    49679   GGgCatCCTGTCGCaGACaCAGG +       5
    GGCCGACCTGTCGCTGACGCNNN chr8    517739  GcCCtgCaTGTgGCTGACGCAGG +       5
    GGCCGACCTGTCGCTGACGCNNN chr8    599935  tGCCGtCtTcTCcCTGACGCCAG -       5
    GGCCGACCTGTCGCTGACGCNNN chr8    5308348 GGCaGgCCTGgCttTGACGCAGG -       5
    GGCCGACCTGTCGCTGACGCNNN chr8    9525579 GGCCcAgCTGTtGCTGAtGaAAG +       5
    GGCCGACCTGTCGCTGACGCNNN chr8    12657177        GGCCcACCTGTgGCTGcCcaTAG -       5
    GGCCGACCTGTCGCTGACGCNNN chr8    12808911        GGCCGACCaGgtGCTccCGCCGG +       5
    GGCCGACCTGTCGCTGACGCNNN chr8    21351922        GGCCcACCTGaCtCTGAgGaCAG -       5
    GGCCGACCTGTCGCTGACGCNNN chr8    21965064        GGCCGtCCTGcgGCTGctGCAGG -       5
    GGCCGACCTGTCGCTGACGCNNN chr8    22409058        GcCCGACCccTCcCcGACGCCAG +       5
    ...


Installation
----------------

* Cas-OFFinder requires OpenCL device and corresponding SDK to run properly.

  OpenCL is supported on varoius platforms, including recent Intel/AMD CPUs and NVidia/AMD graphic cards!
  Before installing Cas-OFFinder, please check the vendor's website to know whether your device is supported.
  
  Now, download OpenCL SDK. If you know your device's vendor name, it is enough that to install only your vendor's one.
  - AMD: http://developer.amd.com/tools-and-sdks/heterogeneous-computing/amd-accelerated-parallel-processing-app-sdk/downloads/
  - Intel: http://software.intel.com/en-us/vcsource/tools/opencl-sdk
  - NVidia: https://developer.nvidia.com/cuda-downloads

* Download Cas-OFFinder binary, or compile it from its source code (below section).

Compilation
----------------
* OpenCL library is required to compile Cas-OFFinder.
  
  To support cross-platform compilation on various operating systems,
  CMake build system is used (more informations on http://www.cmake.org).
  
  First, download CMake here (http://www.cmake.org/cmake/resources/software.html).
  If you use Ubuntu linux, you can also install it via apt-get.
  (apt-get install cmake)
  
  Checkout the source code of Cas-OFFinder with Git client,
  or download it manually on github website.
  
  Now, type the following in termnal (on windows, press windows_key+r and type cmd) to build Cas-OFFinder.

  In POSIX environment (g++ should be pre-installed):
  
      cmake -G "Unix Makefiles"
      make
      
  On Windows (Visual Studio should be pre-installed):
  
      cmake -G "NMake Makefiles"
      nmake
  
  Then cas-offinder binary will be generated. Copy it wherever you want.

Module reference
----------------

* For reading/parsing FASTA files, the kseq.h library (developed by Heng Li) is used.

  The library is distributed under MIT licence.

  More informations on:
  http://lh3lh3.users.sourceforge.net/parsefastq.shtml

* For supporting Dirent API on Windows environment,
  Dirent API for Microsoft Visual Studio is used.

  More informations on:
  http://softagalleria.net/dirent.php

Download & Source
--------
The binaries and the source code can be downloaded from

https://github.com/snugel/cas-offinder

License
-------
Cas-OFFinder (Except kseq.h and dirent.h) is licensed under the new BSD licence.

Copyright (c) 2013, Jeongbin Park and Sangsu Bae
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the name of the organization nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
