import subprocess
from pprint import pprint
import difflib
import tempfile
import sys
import random

def line_cmp(line):
    els = line.split(b'\t') 
    if len(els) != 9:
        print(line,file=sys.stderr)
        return (line,)
    pid,mtcty,rna,dna,chrom,loc,dir,mis,blgs = els[:9]
    return (chrom, int(loc), pid, dir, mtcty, int(blgs),int(mis), rna, dna)

def sort_data(bdata):
    lines = bdata.splitlines(keepends=True)
    headers = lines[:2]
    lines = lines[2:]
    lines.sort(key=line_cmp)
    return b"".join(headers + lines)

def compare_on_input(input_filename, device):
    args = [input_filename, device, '-']
    out2 = subprocess.run(["./build/cas-offinder-3"] + args,stdout=subprocess.PIPE)
    print("finishedone",flush=True)
    out1 = subprocess.run(["./bin/cas-offinder-3"] + args,stdout=subprocess.PIPE)
    print("finishedboth",flush=True)    


    with tempfile.NamedTemporaryFile() as file1, \
       tempfile.NamedTemporaryFile() as file2:
        file1.write(sort_data(out1.stdout))
        file1.flush()
        file2.write(sort_data(out2.stdout))
        file2.flush()
        out = subprocess.run(["diff",file1.name, file2.name])
        if out.returncode == 0:
            print("files are the same")
        else:
            print("files differed")     


if __name__ == "__main__":
    assert len(sys.argv) == 3, "requires 2 arguments: input file, device"
    compare_on_input(sys.argv[1], sys.argv[2])
