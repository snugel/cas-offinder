import subprocess
import os
import sys

def download_test_set():
    if not os.path.exists("_test_data"):
        os.mkdir("_test_data")
    if not os.path.exists("_test_data/hg38.fa"):
        if not os.path.exists("_test_data/hg38.fa.gz"):
            print("downloadin fa dataset....",file=sys.stderr)
            subprocess.run(["wget","ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"])
            os.rename("hg38.fa.gz","_test_data/hg38.fa.gz")
        print("unpacking fa dataset....",file=sys.stderr)
        subprocess.run(["gzip", "-d", "_test_data/hg38.fa.gz"])
    else:
        print("using cached fa files")
    if not os.path.exists("_test_data/hg38.2bit"):
        print("downloadin 2bit dataset....",file=sys.stderr)
        subprocess.run(["wget","ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit"])
        os.rename("hg38.2bit","_test_data/hg38.2bit")
    else:
        print("using cached 2bit files")

if __name__ == "__main__":
    download_test_set()
