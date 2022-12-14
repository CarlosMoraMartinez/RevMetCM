# CM 221001: Code checked, some output checked, seems ok but needs refafctor

#!/usr/bin/env python
import argparse
import sys
import os

# ----- command line parsing -----
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='A tool for uniquely binning reads based on percentage coverage')
parser.add_argument('pc_file', type=str, help='Concatenated percentage coverage file to be binned.')
parser.add_argument('nanopore_read_ids', type=str, help='The total list of read ids for the barcode.')

def main():
    args = parser.parse_args()
    # ----- end command line parsing -----
    
    pc_file = open(args.pc_file)
    barcode_ids = open(args.nanopore_read_ids)
    
    dict_res = {}
    
    for id in barcode_ids:
        id = id.strip()
        dict_res[id] = {'unassigned': 0}
    
    for line in pc_file:
        line = line.strip()
        illumina_ref,minion_read_id,pc = line.split(' ')
        if minion_read_id == "none":
            continue
        pc = float(pc)
        if pc > dict_res[minion_read_id][list(dict_res[minion_read_id].keys())[0]]:
            old_ref_key = list(dict_res[minion_read_id].keys())[0]
            dict_res[minion_read_id][illumina_ref] = dict_res[minion_read_id][old_ref_key]
            del dict_res[minion_read_id][old_ref_key]
            dict_res[minion_read_id][illumina_ref] = pc
    
    for key in dict_res:
        illumina_ref_out = list(dict_res[key].keys())[0]
        pc_out = float(dict_res[key][illumina_ref_out])
        sys.stdout.write("{0} {1} {2:.3f}\n".format(key, illumina_ref_out, pc_out))
    
if __name__ == "__main__":
    main()