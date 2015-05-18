import os
import glob
from concurrent.futures import ProcessPoolExecutor

def printContigs(filepath):
    with open(filepath) as handle:
        for line in handle:
            if line.startswith('>'):
                a = 'a'
                #print('{0}: {1}'.format(os.path.basename(filepath), line))

def main():
    for x in map(printContigs, glob.glob("./benchmarks/Ecoli/external/*.frankenfasta")):
        pass

#timeit.timeit(stmt='pass', setup='pass', timer=<default timer>, number=1000000)

if __name__ == "__main__":
    import timeit
    print(timeit.repeat("main()", setup="from __main__ import main", repeat=5, number=1))
