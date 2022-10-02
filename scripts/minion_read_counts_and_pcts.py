# CM 221002: Revised and refactored; merged with convert_counts_to_percentages
#!/usr/bin/env python
import argparse
import sys
import os
import pandas as pd

def get_illumina_ref(illumina_reference_list):
    """Reads a list with the names of all Illumina references."""
    illumina_reference_list = open(illumina_reference_list, 'r')
    dict_res = {id.strip():0 for id in illumina_reference_list}
    dict_res.setdefault("unassigned", 0)
    illumina_reference_list.close()
    return dict_res    

def read_binned_pcs(binned_pc_file):
    """Reads barcode .binned file, which contains the Illumina reference 
    assigned to each nanopore read and the % of the nanopore read covered"""
    binned_pc_file = open(binned_pc_file, 'r')
    pcs = [line.strip().split(' ') for line in binned_pc_file]
    binned_pc_file.close()
    return pcs

def calc_counts(pcs, dict_res, min_pc, max_pc):   
    """Number of nanopore reads assigned to each Illumina reference"""   
    for minion_read_id, illumina_ref, pc in pcs:
        if float(pc) < min_pc or float(pc) > max_pc:
            dict_res["unassigned"] += 1
        else:
            dict_res[illumina_ref] = dict_res[illumina_ref] +1 if illumina_ref in dict_res.keys() else 1 
    return dict_res

def calc_percents(dict_res, threshold):
    total = float(sum(c for c in dict_res.values()))
    total = total if total > 0 else 1
    total_pct = [c/total for c in dict_res.values()]
    
    total_assigned = total - dict_res["unassigned"]
    total_assigned = total_assigned if total_assigned > 0 else 1
    total_assigned_pct = [c/total_assigned if k != "unassigned" else 0 for k, c in dict_res.items()]
    
    #If less than threshold_count reads are assigned to some Illumina reference, they will count as 0
    threshold_count = threshold*total_assigned/100
    threshold_corrected_counts = [c if c > threshold_count and k != "unassigned" else 0 for k, c in dict_res.items()]
    threshold_corrected_total = sum(threshold_corrected_counts)
    threshold_corrected_total = threshold_corrected_total if threshold_corrected_total > 0 else 1
    threshold_corrected_pct = [c/threshold_corrected_total for c in threshold_corrected_counts]
    
    new_res = pd.DataFrame.from_dict({"illumina_ref":list(dict_res.keys()),
                        "ont_reads_assigned":list(dict_res.values()),
                        "percent_total": total_pct,
                        "percent_assigned":total_assigned_pct,
                        "count_threshold_corrected":threshold_corrected_counts,
                        "percent_assigned_threshold":threshold_corrected_pct})
    return new_res
    
        
def write_percents(res: pd.DataFrame, outfile):
    ont_name = outfile.split("_")[0]
    res.insert(0,"NanoporeSample", [ont_name for i in range(res.shape[0])])
    res.to_csv(outfile, sep="\t")
    
# ----- command line parsing -----
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, 
                                 description='Count occurrences of each skim reference in binned percent coverage files')
parser.add_argument('binned_pc_file', type=str, help='Binned percentage coverage file.')
parser.add_argument('illumina_reference_list', type=str, help='List of Illumina reference IDs including "unassigned".')
parser.add_argument('output_file', type=str, help='Output file')
parser.add_argument('min_pc', type=float, default=15, help='Minimum percent coverage.')
parser.add_argument('max_pc', type=float, default=100, help='Maximum percent coverage.')
parser.add_argument('threshold', type=float, default=1,
                    help='Bin threshold percentage. Species with fewer than this percentage of the total assigned reads will be set to 0.')



def main():
    args = parser.parse_args()
        
    min_pc = args.min_pc
    max_pc = args.max_pc
    threshold = args.threshold
    args.illumina_reference_list
    illumina_reference_list = args.illumina_reference_list
    binned_pc_file = args.binned_pc_file
    outfile = args.output_file
    
    dict_res = get_illumina_ref(illumina_reference_list)
    pcs = read_binned_pcs(binned_pc_file)
    dict_res = calc_counts(pcs, dict_res, min_pc, max_pc)
    df_res = calc_percents(dict_res, threshold)
    write_percents(df_res, outfile)
      
if __name__ == "__main__":
    main()