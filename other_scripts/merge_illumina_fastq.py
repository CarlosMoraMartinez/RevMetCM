import os
import glob


INDIR = "/scratch/groups/assembly/shared/projects/pollergen/data/illumina_genome_skim_EASI_48/*/*gz"
OUTDIR = "/scratch/groups/assembly/shared/cmora/easi_illumina"

def copy2outdir(sample, files, outdir):
    for f in files:
        try:
            cmd = f"cp -v {f} {outdir}"
            os.system(cmd)
        except:
            print(f"Unable to copy {f}") 

def mergeall(samples, indir, outdir):
    for s in samples:
        print(f"\n***\n{s}: {indir}\n***\n")
        r1 = glob.glob(f"{indir}/{s}.*idt-UMI.1.*gz")
        r2 = [i.replace("idt-UMI.1.", "idt-UMI.2.") for i in r1]
        namer1 = os.path.basename(r1[0]).replace("idt-UMI.1.", "idt-UMI.R1.").replace(".", "_").replace("_fastq_gz", ".fastq")
        namer2 = namer1.replace("_R1", "_R2")
        cmd1 = f"zcat {' '.join(r1)} > {outdir}{namer1}"
        cmd2 = f"zcat {' '.join(r2)} > {outdir}{namer2}"
        cmd3 = f"gzip {outdir}{namer1}"
        cmd4 = f"gzip {outdir}{namer2}"
        for cmd in [cmd1, cmd2, cmd3, cmd4]:
            try:
                print(cmd)
                os.system(cmd)
            except:
                print(f"Error executing:\n{cmd}")
        
        
    
    
def main():
    try:
        os.system(f"mkdir {OUTDIR}")
    except:
        print("Unable to create output directory")
    try:
        os.system(f"mkdir {OUTDIR}/merged")
    except:
        print("Unable to create output directory")
    fnames = glob.glob(INDIR)
    samples = set([os.path.basename(i).split(".")[0] for i in fnames])
    samples = {s:[f for f in fnames if os.path.basename(f).startswith(s + '.')] for s in samples}
    #for k, v in samples.items():
    #    copy2outdir(k, v, OUTDIR)
    mergeall(samples.keys(), OUTDIR, OUTDIR + "/merged/")
      
if __name__ == "__main__":
    main()