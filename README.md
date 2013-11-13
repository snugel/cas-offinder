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

First, download any target organism's chromosome FASTA files. You can find one in below link:

http://hgdownload.soe.ucsc.edu/downloads.html (UCSC genome sequences library)

or http://ensembl.org/info/data/ftp/index.html (Ensembl sequence library)
  
Untar and ungzip them in a directory.

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

For example, like this:

    /var/chromosomes/human_hg19
    NNNNNNNNNNNNNNNNNNNNNGG
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

    GGCCGACCTGTCGCTGACGCNNN	Chr1	2643000	GGttGACCTGTCGgTGAgcCCAG	+	5
    GGCCGACCTGTCGCTGACGCNNN	Chr1	27772702	GGCCtggCTGTgcCTGACGCAAG	-	5
    CGCCAGCGTCAGCGACAGGTNNN	Chr1	5138727	CGgCAGCtcagGCGACAGGTGAG	-	5
    CGCCAGCGTCAGCGACAGGTNNN	Chr1	7671026	gGCCAGCagCAGCcACAGaTGAG	-	5
    CGCCAGCGTCAGCGACAGGTNNN	Chr1	22212625	CGCagGaGTCAatGACAGGTAAG	-	5
    ACGGCGCCAGCGTCAGCGACNNN	Chr1	689718	ttGGCaCCAGCaTCAGCtACAAG	+	5
    ACGGCGCCAGCGTCAGCGACNNN	Chr1	4445225	ACtGaGCCAGtGTaAGCcACAAG	+	5
    ACGGCGCCAGCGTCAGCGACNNN	Chr1	8874145	ACGGCaCCAGCtTCAGtaAaCAG	+	5
    ACGGCGCCAGCGTCAGCGACNNN	Chr1	16784268	AgaGCtCCAGCtTCAGCGAtAAG	-	5
    ACGGCGCCAGCGTCAGCGACNNN	Chr1	24062269	AtGGCtCCAGCGgCAGCatCAAG	+	5
    ACGGCGCCAGCGTCAGCGACNNN	Chr1	28209719	cCGGCGtCAtCGTCAtCGAtGAG	+	5
    ACGGCGCCAGCGTCAGCGACNNN	Chr1	29108473	ACaGCtCCAGCtTCAGCcAtAAG	-	5
    GTCGCTGACGCTGGCGCCGTNNN	Chr1	222638	GTCGCTGACGCaGcCcgCaTGAG	-	5
    GTCGCTGACGCTGGCGCCGTNNN	Chr1	3640018	GTCGtTGcCGCcGtCGgCGTGAG	-	5
    ...


Installation
----------------

* Cas-OFFinder requires OpenCL device to run properly.

  OpenCL is supported on varoius platforms, including recent Intel/AMD CPUs and NVidia/AMD graphic cards!
  Before installing Cas-OFFinder, please check whether your device is supported on the vendor's website.
  
  Now, download OpenCL SDK. If you know your device's vendor name, it is enough that to install only your vendor's one.
  - AMD: http://developer.amd.com/tools-and-sdks/heterogeneous-computing/amd-accelerated-parallel-processing-app-sdk/downloads/
  - Intel: http://software.intel.com/en-us/vcsource/tools/opencl-sdk
  - NVidia: https://developer.nvidia.com/cuda-downloads

* Download Cas-OFFinder binary, or compile it from its source code (below section).

Compile
----------------
* After install OpenCL library, now you can compile Cas-OFFinder.
  
  To support cross-platform compilation on various operating systems,
  CMake build system is used (more informations on http://www.cmake.org).
  
  First, download CMake here (http://www.cmake.org/cmake/resources/software.html).
  If you use Ubuntu linux, you can also install it via apt-get.
  (apt-get install cmake)
  
  Checkout the source code of Cas-OFFinder with Git client,
  or download it manually on github website.
  
  Now, type the following to build Cas-OFFinder.

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
  The kseq.h library is distributed under MIT licence.

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
Cas-OFFinder is licensed under the new BSD licence.

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
