#!/usr/bin/env python3

__author__ = "David Smith"
__version__ = "0.9.1"
__email__ = "dsmith@tgen.org"


def _parse_args():
    'Parses command line arguments.'
    import argparse
    parser = argparse.ArgumentParser( description="Meant to be called from the pipeline automatically." )
    parser.add_argument( "--nucmerpath", default="nucmer", help="Path to the 'nucmer' executable." )
    parser.add_argument( "--reference", required=True, help="Path to the reference fasta file." )
    parser.add_argument( "--outputfile", required=True, help="Path to the dups file to output." )
    return( parser.parse_args() )

def run_nucmer_on_reference( nucmer_path, reference_path ):
    import subprocess
    import sys
    return_code = subprocess.call( [ nucmer_path, "--prefix=reference", "--maxmatch", "--nosimplify", reference_path, reference_path ] )
    if return_code > 0:
        sys.stderr.write( "NASP WARNING: nucmer may have encountered errors during reference duplicates checking, proceeding anyway\n" )

def _parse_delta_line( line_from_delta_file, dups_data, current_contigs ):
    import re
    line_match = re.match( '^>([^ ]+) ([^ ]+) (\d+) (\d+)\s*$', line_from_delta_file )
    if line_match:
        current_contigs = ( line_match.group(1), line_match.group(2) )
        contig_0_start = int( line_match.group(3) )
        contig_1_start = int( line_match.group(4) )
        dups_data.add_contig( current_contigs[0], False )
        dups_data.add_contig( current_contigs[1], False )
        if contig_0_start > dups_data.get_contig_length( current_contigs[0] ):
            dups_data.append_contig( ( "0" * ( contig_0_start - dups_data.get_contig_length( current_contigs[0] ) ) ), current_contigs[0] )
        if contig_0_start > dups_data.get_contig_length( current_contigs[1] ):
            dups_data.append_contig( ( "0" * ( contig_1_start - dups_data.get_contig_length( current_contigs[1] ) ) ), current_contigs[1] )
    else:
        line_match = re.match( '^(\d+) (\d+) (\d+) (\d+) \d+ \d+ \d+\s*$', line_from_delta_file )
        if line_match:
            contig_0_start = int( line_match.group(1) )
            contig_0_end = int( line_match.group(2) )
            contig_1_start = int( line_match.group(3) )
            contig_1_end = int( line_match.group(4) )
            if ( current_contigs[0] != current_contigs[1] ) or ( contig_0_start != contig_1_start ):
                if contig_0_end < contig_0_start:
                    contig_0_end, contig_0_start = contig_0_start, contig_0_end
                if contig_1_end < contig_1_start:
                    contig_1_end, contig_1_start = contig_1_start, contig_1_end
                dups_data.set_value( ( "1" * ( contig_0_end - contig_0_start + 1 ) ), current_contigs[0], contig_0_start, "?" )
                dups_data.set_value( ( "1" * ( contig_1_end - contig_1_start + 1 ) ), current_contigs[1], contig_1_start, "?" )
    return current_contigs

def parse_delta_file( delta_filename, dups_data ):
    current_contigs = ( "", "" )
    delta_handle = open( delta_filename, 'r' )
    for line_from_delta_file in delta_handle:
        current_contigs = _parse_delta_line( line_from_delta_file, dups_data, current_contigs )
    delta_handle.close()

def write_dups_file( output_filename, dups_data ):
    dups_data.write_to_file( output_filename )

def main():
    commandline_args = _parse_args()
    from nasp_objects import GenomeStatus
    run_nucmer_on_reference( commandline_args.nucmerpath, commandline_args.reference )
    dups_data = GenomeStatus()
    parse_delta_file( "reference.delta", dups_data )
    write_dups_file( commandline_args.outputfile, dups_data )

if __name__ == "__main__": main()


