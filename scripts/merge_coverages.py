import argparse
import glob
from functools import partial
from typing import List, Callable
from multiprocessing import Pool

import pandas as pd

INPUT_EXTENSION: str = '.depth.gz'
DEFAULT_PROCESSES: int = 12
UNASSIGNED: str = 'Unassigned'

def read_table(fname: str, extension:str = INPUT_EXTENSION) -> pd.DataFrame:
    ont: str
    illumina: str
    ont, illumina = fname.replace(extension, "").split("_")
    ont = ont.split(".")[0]
    illumina = illumina.split(".")[0]   
    df : pd.DataFrame = pd.read_csv(fname, compression="gzip", sep="\t")
    df.insert(0, 'ONT', [ont]*df.shape[0])
    df.insert(1, 'ILLUMINA', [illumina]*df.shape[0])
    return df

def read_all(fnames: List[str], n_processes: int = DEFAULT_PROCESSES, extension:str = INPUT_EXTENSION) -> pd.DataFrame:
    with Pool(n_processes) as p:
        res: List[pd.DataFrame] = p.map(partial(read_table, extension=extension), fnames)
    return pd.concat(res)

bin_read: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x.loc[x.coverage == max(x.coverage)] if any(x.coverage) else x.iloc[:1]

def bin2species(df: pd.DataFrame) -> pd.DataFrame:
    dfbinned = df.groupby(['ONT', '#rname']).apply(bin_read).reset_index(drop=True)
    dfbinned.loc[dfbinned.coverage == 0.0, 'ILLUMINA'] = UNASSIGNED
    return dfbinned
    
# ----- command line parsing -----
parser: argparse.ArgumentParser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, 
                                                          description='A tool for uniquely binning reads based on percentage coverage')
parser.add_argument('-n', '--n_cores', type = int, help='Number of parallel processes', default=DEFAULT_PROCESSES)
parser.add_argument('-p', '--prefix', type = str, help='Prefix of all input files', default="")
parser.add_argument('-e', '--extension', type = str, help='Extension input files', default=INPUT_EXTENSION)


def main():
    args = parser.parse_args()
    n_cores: int = args.n_cores
    prefix: str = args.prefix
    extension: str = args.extension
    
    fnames: List[str] = glob.glob(f"{prefix}*{INPUT_EXTENSION}")
    df: pd.DataFrame = read_all(fnames, n_cores, extension)
    df = bin2species(df)
    df.to_csv(f"{prefix}_all.csv", sep="\t", index=False)
  
if __name__ == "__main__":
    main()