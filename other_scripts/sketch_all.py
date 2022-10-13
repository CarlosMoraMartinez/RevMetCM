import os
import subprocess
import glob


paramsont = "/home/carmoma/projects/pollen/revmet_paper_data/ena_files_nanopore/ERR3132323/mock_mixes_nanopore_reads/mock_mixes/barcode*.fastq"
paramsillumina = "/home/carmoma/projects/pollen/revmet_paper_data/ena_files_illumina/*.fastq.gz"

ont_out = "/home/carmoma/projects/pollen/sketch1/ont_sketch/"
illumina_out = "/home/carmoma/projects/pollen/sketch1/illumina_sketch/"
out_dist = "/home/carmoma/projects/pollen/sketch1/dist_all/"

def run(cmd):
    print(cmd)
    try:
        p = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
        p.wait()
    except:
        print(f"Error running {cmd}")
    
def getSamplesDict(directory):
    samples = dict()
    for i in glob.iglob(directory):
        print(i)
        sample = i.split("/")[-1].split("_")[0]
        if sample in samples.keys():
            samples[sample].append(i)
        else:
            samples.setdefault(sample, [i])
    return samples

def sketch(samples: dict, outdir, rm=True):
    for k, v in samples.items():
        vv = " ".join(v)
        newfname = f"{k}_join.fastq"
        print(f"Sketching {k}...")
        if v[0].endswith(".gz"):
            cmd = f"zcat {vv} > {newfname}"
            run(cmd)
        else:
            if len(v) > 1:
                cmd = f"cat {vv} > {newfname}"
                run(cmd)
            else:
                newfname = vv    
        outname = outdir + k
        cmd = f"mash sketch {newfname} -o {outname}"
        run(cmd)
        if rm:
            cmd = f"rm {k}_join.fastq"
            run(cmd)

def dist_all(out1, out2):
    f1 = glob.glob(out1 + "/*.msh")
    f2 = glob.glob(out2 + "/*.msh")
    for f in f1:
        sample = f.split("/")[-1].split("_")[0]
        cmd = f"mash dist {f} {out1}*.msh > {out_dist}{sample}.same.dist"
        run(cmd)
    for f in f2:
        sample = f.split("/")[-1].split("_")[0]
        cmd = f"mash dist {f} {out2}*.msh > {out_dist}{sample}.same.dist"    
        run(cmd)
    for f in f2:
        sample = f.split("/")[-1].split("_")[0]
        cmd = f"mash dist {f} {out1}*.msh > {out_dist}{sample}.other.dist"        
        run(cmd)
def main():
    #ss = getSamplesDict(paramsont)
    #print(ss)
    #sketch(ss, ont_out)
    #ss = getSamplesDict(paramsillumina)
    #print(ss)
    #sketch(ss, illumina_out)
    dist_all(illumina_out, ont_out)
if __name__ == "__main__":
    main()
            
    
