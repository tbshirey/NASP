#!/usr/bin/env python

"""filter a raw ISG matrix for a
list of genomes"""

from optparse import OptionParser
from collections import deque

def filter_genomes(genomes, in_matrix):
    #genomes_file = open(genomes, "rU")
    in_matrix = open(in_matrix, "rU")
    #all_genomes = [ ]
    firstLine = in_matrix.readline()
    first_fields = firstLine.split()
    last=first_fields.index("#SNPcall")
    all_genomes=first_fields[:last]
    genomes_file = open(genomes, "r").read().splitlines()
    to_keep = [ ]
    for x in all_genomes[1:]:
        if x in genomes_file:
            to_keep.append(all_genomes.index(x))
    return to_keep
    in_matrix.close()

def filter_matrix(to_keep, in_matrix, prefix):
    matrix = open(in_matrix, "rU")
    outfile = open("%s_genomes.matrix" % prefix, "w")
    for line in matrix:
        fields = line.split()
        deque((list.pop(fields, i) for i in sorted(to_keep, reverse=True)), maxlen=0)
        print >> outfile, "\t".join(fields)
    outfile.close()

def main(in_matrix, prefix, genomes):
    to_keep=filter_genomes(genomes, in_matrix)
    filter_matrix(to_keep, in_matrix, prefix)

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-m", "--in_matrix", dest="in_matrix",
                      help="/path/to/file [REQUIRED]",
                      action="store", type="string")
    parser.add_option("-p", "--out_prefix", dest="prefix",
                      help="/path/to/file [REQUIRED]",
                      action="store", type="string")
    parser.add_option("-g", "--genomes", dest="genomes",
                      help="/path/to/genomes_file [REQUIRED]",
                      action="store", type="string")
    parser.add_option("-a", "--actions", dest="action",
                      help="action to perform (keep, remove), defaults to keep",
                      action="store", type="string", default="keep")
    options, args = parser.parse_args()

    mandatories = ["in_matrix", "prefix", "genomes"]
    for m in mandatories:
        if not options.__dict__[m]:
            # FIXME(jtravis): replace print with print()q
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.in_matrix,options.prefix,options.genomes,options.action)
