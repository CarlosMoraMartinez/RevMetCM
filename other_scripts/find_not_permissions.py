import os
import glob

INDIR = "/scratch/groups/assembly/shared/projects/pollergen/data/illumina_genome_skim_EASI_48/*/*gz"

def copy2outdir(sample, files):
    fo = open("bad_files.txt", 'w')
    for f in files:
        try:
            cmd = f"zcat {f} | head -1 > err.txt"
            os.system(cmd)
        except:
            fo.write(f"{sample}\t{f}") 
    fo.close()
    
def main():
    fnames = glob.glob(INDIR)
    samples = set([os.path.basename(i).split(".")[0] for i in fnames])
    samples = {s:[f for f in fnames if os.path.basename(f).startswith(s + '.')] for s in samples}
    for k, v in samples.items():
        copy2outdir(k, v)
      
if __name__ == "__main__":
    main()