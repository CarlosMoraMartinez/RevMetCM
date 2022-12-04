import argparse
import glob
from functools import partial
from typing import List, Callable
import os
from multiprocessing import Pool

import pandas as pd

DEFAULT_MIN_PERCENTAGE: float = 15.0
INPUT_EXTENSION: str = '.depth.gz'
DEFAULT_PROCESSES: int = 12
UNASSIGNED: str = 'Unassigned'


def read_filelist(txtfile: str, prefix: str) -> List[str]:
    fnames: List[str] = [l.strip() for l in open(txtfile, 'r')]
    fnames = [f for f in fnames if os.path.basename(f).startswith(prefix)]
    return fnames

def add_unassigned_reads(fname, df: pd.DataFrame) -> pd.DataFrame:
    print(f"Reading {fname} to get read list")
    df_case : pd.DataFrame = read_table(fname, -1)
    reads: List = df_case['#rname'].tolist()
    assigned_reads: List[str] = df['#rname'].tolist()
    unassigned_reads: List[str] = [i for i in reads if not i in assigned_reads]
    df_case_ind: pd.DataFrame = df_case.set_index('#rname').loc[unassigned_reads].reset_index()
    df_case_ind = df_case_ind[df_case.columns]
    df_case_ind.ILLUMINA = UNASSIGNED
    for col in df_case_ind.columns[5:]:
        df_case_ind[col] = 0
    df_case_ind = pd.concat([df, df_case_ind])    
    return df_case_ind
  
def read_table(fname: str, min_percentage: float = DEFAULT_MIN_PERCENTAGE, extension:str = INPUT_EXTENSION) -> pd.DataFrame:
    print(f"Reading {fname}")
    ont: str
    illumina: str
    ont, illumina, *_ = os.path.basename(fname).replace(extension, "").split("_")
    ont = ont.split(".")[0]
    illumina = illumina.split(".")[0]   
    df : pd.DataFrame = pd.read_csv(fname, compression="gzip", sep="\t")
    if min_percentage > 0.0:
        df = df.loc[df.coverage >= min_percentage]
    df.insert(0, 'ONT', [ont]*df.shape[0])
    df.insert(1, 'ILLUMINA', [illumina]*df.shape[0])
    print(f"Read {fname} with {df.shape[0]} rows after filtering.")
    return df

def read_all(fnames: List[str], min_percentage: float = DEFAULT_MIN_PERCENTAGE, n_processes: int = DEFAULT_PROCESSES, extension:str = INPUT_EXTENSION) -> pd.DataFrame:
    with Pool(min([n_processes, len(fnames)])) as p:
        res: List[pd.DataFrame] = p.map(partial(read_table, min_percentage=min_percentage, extension=extension), fnames)
    print("All files read")
    df: pd.DataFrame = pd.DataFrame()
    try:
        df = pd.concat(res)
    except ValueError:
        print(f"Error: No DataFrames to concat. Files were: {fnames}")
    return df

#bin_read: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x.loc[x.coverage == max(x.coverage)] if any(x.coverage) else x.iloc[:1]
def bin_read(x: pd.DataFrame) -> pd.DataFrame:
    return x.loc[x.coverage == max(x.coverage)] if any(x.coverage) else x.iloc[:1]    
    
def bin2species(df: pd.DataFrame, n_cores: int = DEFAULT_PROCESSES) -> pd.DataFrame:
    print("Binning to species")
    dfbinned = df.groupby(['#rname']) # .apply(bin_read).reset_index(drop=True) # Too slow. Parallelize
    dflist = [dfbinned.get_group(g) for g in dfbinned.groups]
    with Pool(n_cores) as p:
        res: List[pd.DataFrame] = p.map(bin_read, dflist)
    print("Binned to species")
    dfbinned = pd.concat(res)
    print("Binned concatenated")
    dfbinned.loc[dfbinned.coverage == 0.0, 'ILLUMINA'] = UNASSIGNED
    return dfbinned

def writeOutput(df: pd.DataFrame, prefix: str, min_percent: float, outdir:str = "") -> None:
    outfile: str = f"{prefix}_minp{int(min_percent)}_all.csv"
    if outdir:
        outfile = outdir + '/' + outfile
        try:
            os.system(f"mkdir {outdir}")
        except FileExistsError:
            print(f"Directory already exists: {outdir}")
    print(f"Printing output: {outfile}")
    df.to_csv(outfile, sep="\t", index=False)
       
# ----- command line parsing -----
parser: argparse.ArgumentParser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, 
                                                          description='A tool for uniquely binning reads based on percentage coverage')
parser.add_argument('-m', '--min_percent', type = float, help='Min percentage of ONT read coverage', default=DEFAULT_MIN_PERCENTAGE)
parser.add_argument('-n', '--n_cores', type = int, help='Number of parallel processes', default=DEFAULT_PROCESSES)
parser.add_argument('-p', '--prefix', type = str, help='Prefix of all input files', default="")
parser.add_argument('-f', '--filelist', type = str, help='Text file with full paths of files to read', default="")
parser.add_argument('-e', '--extension', type = str, help='Extension input files', default=INPUT_EXTENSION)
parser.add_argument('-o', '--outdir', type = str, help='Output directory', default="")


def main():
    args = parser.parse_args()
    min_percent: float = args.min_percent
    n_cores: int = args.n_cores
    prefix: str = args.prefix
    extension: str = args.extension
    filelist: str = args.filelist
    outdir: str = args.outdir
    
    fnames: List[str]
    if not filelist:
        print("Looking for files in the workspace")
        fnames = glob.glob(f"{prefix}*{INPUT_EXTENSION}")
    else:
        print("Reading list of files")
        fnames = read_filelist(filelist, prefix)
    
    df: pd.DataFrame = read_all(fnames, min_percent, n_cores, extension)
    if df.shape[0] == 0:
        exit(0)
    df = bin2species(df, n_cores)  
    #df = add_unassigned_reads(fnames[0], df)
    writeOutput(df, prefix, min_percent, outdir)
  
if __name__ == "__main__":
    main()