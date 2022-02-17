try:
    import h5py
except ModuleNotFoundError:
    print("h5py module not found, activate appropiate enviroment or install with conda")
import os
import glob
import sys
import multiprocessing as mp
import argparse
from itertools import repeat

# Parse arguments
# Instantiate the parser
parser = argparse.ArgumentParser(description='Extracts fastq sequences from MULTI fast5 files to a file named the same as the fast5 (ont_fast5_api only takes SINGLE fast5...)')
# Required positional arguments
parser.add_argument('fast5',
                    help='Path to folder containing fast5 files')
parser.add_argument('fastq',
                    help='Path to folder in which to output fastq files')
# Optional arguments
parser.add_argument('--batch_size', nargs=1, type=int, default=4000,
                    help="Number of reads to process at a time from each multi fast5 file")
parser.add_argument('--jobs', nargs=1, type=int, default=100,
                    help="Number of fast5 files to process in parallel")
# Optinal switches
parser.add_argument('--overwrite', action='store_true',
                    help="Overwrite already existing fastq files")
# Load arguments
args = parser.parse_args()

def extract_fastq(f5_file, fq_path, batch_size = args.batch_size[0]):
    # Check if fastq has already been produced
    if os.path.exists(fq_path) & args.overwrite:
        print("WARNING! " + str(os.path.basename(fq_path)) + " already exists and will be overwritten.")
    elif os.path.exists(fq_path):
        print("WARNING! " + str(os.path.basename(fq_path)) + " already exists and this file is skipped")
        return None
    
    # Create or overwrite fastq file
    fq_file = open(fq_path, "w")
    fq_file.close()
    
    # open fast5 file
    try:
        F = h5py.File(f5_file, 'r')
    except FileNotFoundError:
        print("  Unable to read fast5 file")
    
    # Ensure that the batch size match the number of reads in file
    assert(len(list(F.keys())) % batch_size == 0), "Please set batch_size as a divisor of " + str(en(list(F.keys())))
    
    for i in  range(0, int(4000/batch_size)):
        print("Reading in " + os.path.basename(f5_file) + ". Read " + str(i * batch_size) + " to " + str((i + 1) * batch_size))
        # Extract the fastq files from multi fast5
        try:
            H = [F[read]["Analyses"]['Basecall_1D_000']['BaseCalled_template']["Fastq"][()].decode() for read in list(F.keys())[i*batch_size:((i+1)*batch_size)]]
        except KeyError:
            print("  One of the keys in the Fast5 file was not identified (Make sure the fast5 contains basecalled fastq)")
        except:
            print("  File identified, but something else went wrong")
        
        print("Writing read to fastq")
        # Write the fastqs to fastq file
        try:
            fq_file = open(fq_path, "a")
        except IsADirectoryError:
            print("  The specified fastq file is a directionary (This really shouldn't happen)")
        except PermissionError:
            print("  You do not have writing permission to the fastq file: " + fq_path)
        except:
            print("  Couldn't open output fastq file: " + fq_path)
        
        for read in H:
            fq_file.write(read)
        fq_file.close()
    # Close the fast5 file when done
    F.close()
    return None

def get_fq_path(x, y):
    out = y + "/" + os.path.basename(x) + ".fastq"
    return out

def main():
    f5_dir = args.fast5
    fq_dir = args.fastq
    with mp.Pool(args.jobs[0]) as p:
        try:
            os.mkdir(fq_dir)
        except FileExistsError:
            print("  Fastq output folder already exists")
    
        f5_paths = glob.glob(f5_dir + "/*")
        fq_paths = p.starmap(get_fq_path, zip(f5_paths, repeat(fq_dir)))
        p.starmap(extract_fastq, zip(f5_paths, fq_paths))

if __name__ == "__main__":
    main()
    quit()