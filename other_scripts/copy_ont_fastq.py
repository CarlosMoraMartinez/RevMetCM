import os
import glob


#INDIR = "/scratch/groups/assembly/shared/projects/pollergen/data/ont_mock_x3_EASI_49/*/*gz"
#OUTDIR = "/scratch/groups/assembly/shared/cmora/easi_ont_mock"
INDIR = "/scratch/groups/assembly/shared/projects/pollergen/data/ont_aerial_pollen_EASI_50/*/*fastq.gz"
OUTDIR = "/scratch/groups/assembly/shared/cmora/easi_ont_airsamples"
LIM_SIZE = 2000000

def copy2outdir(files, outdir):
    for f in files:
        try:
            fsize = os.path.getsize(f)
            if fsize > LIM_SIZE:
                cmd = f"cp -v {f} {outdir}"
                os.system(cmd)
            else:
                print(f"Not copying {f}. Size is {fsize}")
        except:
            print(f"Unable to copy {f}") 
    
    
def main():
    try:
        os.system(f"mkdir {OUTDIR}")
    except:
        print("Unable to create output directory")
    fnames = glob.glob(INDIR)
    copy2outdir(fnames, OUTDIR)

      
if __name__ == "__main__":
    main()