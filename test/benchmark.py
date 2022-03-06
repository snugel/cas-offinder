import subprocess
from pprint import pprint
import difflib
import tempfile
import sys
import random
import time
import os
import subprocess


def generate_guide(guide_length, pam_len, mismatches):
    letters = ['A','C','G','T']
    guide = "".join([random.choice(letters) for i in range(guide_length)])
    pam = "N"*pam_len
    return f"{guide}{pam} {mismatches}"


def generate_input(num_guides, guide_length, pam, mismatches):
    endl = "\n"
    guides = [generate_guide(guide_length, len(pam), mismatches) for i in range(num_guides)]
    return f"""_test_data/hg38.2bit
{'N'*guide_length}{pam}
{endl.join(guides)}
"""


def compare_performance(binary, num_guides, guide_length, pam, device, mismatches):
    args = [binary, '-', device, '-']
    open("benchmark.log",'a').write("starting benchmark\n"+' '.join(args)+'\n')
    input = generate_input(num_guides, guide_length, pam, mismatches)

    start = time.time()
    subprocess.run(args, input=input.encode('utf-8'), stdout=subprocess.DEVNULL, stderr=open("benchmark.log",'a'))
    end = time.time()

    return end - start


def build_benchmark_table():
    # clear log file if it exists
    open("benchmark.log",'w').write("")

    guide_lengths = [20, 25]
    mismatches_options = [3, 6]
    pams = ["NRG", "NNGRRT", 'TTTN']
    devices = ['G','C']
    num_guides_list = [1,100]
    binaries = ["./build/cas-offinder-2", "./bin/cas-offinder-2"]
    print("Device | Num Guides | Mismatches | PAM | Guide Length | cas-offinder 2.4.1 | new code ", flush=True)
    print("--- | --- | --- | --- | --- | --- | --- ", flush=True)
    for device in devices:
        for num_guides in num_guides_list:
            for mismatches in mismatches_options:
                for pam in pams:
                    for guide_len in guide_lengths:
                        perfs = []
                        for binary in binaries:
                            perfs.append(compare_performance(binary, num_guides, guide_len, pam, device, mismatches))
                        bolded_perfs = [str(p) for p in perfs]
                        best_idx = perfs.index(min(perfs))
                        bolded_perfs[best_idx] = f"**{bolded_perfs[best_idx]}**"
                        print(f"{device} | {num_guides} | {mismatches} | {pam} | {guide_len} | {bolded_perfs[0]} | {bolded_perfs[1]}", flush=True)

if __name__ == "__main__":
    if not os.path.exists("_test_data/hg38.2bit"):
        print("must call `python download.py` before calling benchmark.py to download the data")
    build_benchmark_table()
