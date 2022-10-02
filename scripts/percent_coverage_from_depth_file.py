# CM 221001: Checked code and some numbers; no bugs; needs refactor.

#!/usr/bin/env python
import argparse
import sys

# ----- command line parsing -----
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='Calculates % coverage for each nanopore read from a samtools depth file.')
parser.add_argument('depth_file', type=str, help='Samtools depth file.')

def perc_cov_calc(non_zero,read_length):
    if read_length == 0:
        return 0
    return 100 * (float(non_zero) / read_length)

def main():
    args = parser.parse_args()
    # ----- end command line parsing -----

    skim_id = args.depth_file.split('_')[1].split(".")[0]
    non_zero = 0
    read_length = 0
    col1_read_id_old = "none"

    open_depth_file = open(args.depth_file)
    for line in open_depth_file:
        col1_read_id, col2_position, col3_depth = line.split('\t')
        col2_position = int(col2_position)
        col3_depth = int(col3_depth)

        if col1_read_id_old == col1_read_id:
            read_length += 1
            if col3_depth != 0:
                non_zero += 1

        else:
            if read_length > 0:
                perc_cov = perc_cov_calc(non_zero, read_length)
                sys.stdout.write("{0} {1} {2:.3f}\n".format(skim_id, col1_read_id_old, perc_cov))
                read_length = 1
                non_zero = 0
                if col3_depth != 0:
                    non_zero += 1
                col1_read_id_old = col1_read_id

            else:
                col1_read_id_old = col1_read_id
                read_length = 1
                if col3_depth != 0:
                    non_zero += 1

    perc_cov = perc_cov_calc(non_zero, read_length)
    sys.stdout.write("{0} {1} {2:.3f}\n".format(skim_id, col1_read_id_old, perc_cov))
    
if __name__ == "__main__":
    main()