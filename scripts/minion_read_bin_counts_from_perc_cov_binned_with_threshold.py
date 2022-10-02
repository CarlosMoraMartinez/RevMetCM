# CM 221001: Revised and modified; needs more refactor
#!/usr/bin/env python
import argparse
import sys
import os

# ----- command line parsing -----
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='Count occurrences of each skim reference in binned percent coverage files')
parser.add_argument('binned_pc_file', type=str, help='Binned percentage coverage file.')
parser.add_argument('illumina_reference_list', type=str, help='List of Illumina reference IDs including "unassigned".')
parser.add_argument('min_pc', type=float, default=15, help='Minimum percent coverage.')
parser.add_argument('max_pc', type=float, default=100, help='Maximum percent coverage.')
#parser.add_argument('outfile', type=float, default=100, help='Output file.')
# ----- end command line parsing -----

def main():
    args = parser.parse_args()
    
    
    min_pc = args.min_pc
    max_pc = args.max_pc
    
    illumina_reference_list = open(args.illumina_reference_list, 'r')
    dict_res = {id.strip():0 for id in illumina_reference_list}
    dict_res.setdefault("unassigned", 0)
    
    illumina_reference_list.close()
    
    binned_pc_file = open(args.binned_pc_file, 'r')
    pcs = [line.strip().split(' ') for line in binned_pc_file]
    binned_pc_file.close()
    
    for minion_read_id, illumina_ref, pc in pcs:
        if float(pc) < min_pc or float(pc) > max_pc:
            dict_res["unassigned"] += 1
        else:
            if illumina_ref in dict_res.keys():
                dict_res[illumina_ref] += 1
            else:
                dict_res[illumina_ref] = 1

    
    for illumina_ref_out in dict_res:
        count_out = int(dict_res[illumina_ref_out])
        sys.stdout.write("{0} {1}\n".format(illumina_ref_out, count_out))
    
if __name__ == "__main__":
    main()