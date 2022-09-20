import sys, getopt, os, time, re, math, glob, gzip, shutil

try:
    sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[:-2]))
except:
    sys.path.insert(1, "/".join(os.getcwd().split("/")[:-1]))

import numpy as np
import pandas as pd

from concurrent.futures import ProcessPoolExecutor

help_message = '''demultiplex.py
Wrong flags on input.
'''


###################
# turn fastqs into fasta
###################
def fastq_to_fasta(input_directory, master_dir):
    # list the fastq files in the directory
    fastqs = [file for file in glob.glob(os.path.join(input_directory, '*')) if bool(re.search("fastq|fq", file))]

    # define combined fastq
    combined_fastq = os.path.join(master_dir, "all_passed_fastqs.fastq.gz")

    # concatenate fastqs
    for i, fastq in enumerate(fastqs):
        # if it is the first fastq, replace results file, otherwise append
        if i == 0:
            cline = f"cat {fastq} > {combined_fastq}"
        else:
            cline = f"cat {fastq} >> {combined_fastq}"

        os.system(cline)

    # define combined output fasta        
    fasta_out = ".".join(combined_fastq.split(".")[:-2]) + ".fasta"

    # define command line
    cline = f"zcat {combined_fastq}" + ''' | awk '{if(NR%4==1) {printf(">%s\\n",substr($0,2));} else if(NR%4==2) print;}' ''' f"> {fasta_out}"
    os.system(cline)

    return fasta_out, combined_fastq


###################
# separate fastas into smaller fastas
###################
def split_fasta(infile, master_dir):
    # max file size is 100Mb (seems to work well on most systems for speed)
    max_bytes = 100000000

    # number of splits for the 
    splits = math.ceil(os.path.getsize(infile) / max_bytes)

    # create directories for each file
    dirs = list(map(lambda x: os.path.join(master_dir, f"working_{x}"), range(splits)))
    for dir in dirs:
        try:
            os.makedirs(f"{dir}")
        except:
            print(f"Making working directory failed for {dir}")

    # total bytes already read for all files, and the total bytes in the file we are processing
    total_bytes = 0
    file_number = 0

    # initial output file
    outfile = open(os.path.join(master_dir, f"working_{file_number}", (f"reads.fa")), "w")

    # separate fasta into parts by memory size
    for i, line in enumerate(open(infile)):
        if i % 2 == 0 and total_bytes > max_bytes:
            # close the file and open the new one
            outfile.close()
            file_number += 1
            total_bytes = 0
            outfile = open(os.path.join(master_dir, f"working_{file_number}", (f"reads.fa")), "w")

        # write line to new file
        outfile.write(line)
        total_bytes += len(line)

    # close the file
    outfile.close()

    return dirs


###################
# read barcodes from metadata and turn into fastqs in each sub directory
###################
def barcodes_from_metadata(metadata_file, master_dir):
    # load the dataframe
    ref_df = pd.read_csv(metadata_file, sep='\t')
    ref_df.columns = list(map(lambda x: re.sub("'", "", x), ref_df.columns.to_numpy()))
    ref_df.columns = list(map(lambda x: re.sub(" ", "_", x), ref_df.columns.to_numpy()))

    # define the parcode types
    barcode_types = ['5_barcode_sequence', '3_barcode_sequence']
    outfiles = []

    # iterate through barcode for assignment
    for spec_type in barcode_types:
        # create a fasta of 5' sequences
        res_fasta = []

        # iterate through barcodes and create a fasta
        for i, barcode in enumerate(np.unique(ref_df.loc[:, spec_type].to_numpy())):
            to_add = f">barcode_{1 + round(len(res_fasta) / 2)}"
            res_fasta.extend([to_add, barcode])

        with open(os.path.join(master_dir, f"{spec_type}.fa"), "w") as outfile:
            outfile.write("\n".join(res_fasta))

        outfile_path = os.path.join(master_dir, f"{spec_type}.fastq")
        outfiles.append(outfile_path)
        # create a fastq
        with open(outfile_path, "w") as outfile:
            for i, seq in enumerate(res_fasta):
                if i % 2 == 0:
                    seq = list(seq)
                    seq[0] = "@"
                    outfile.write(f"{''.join(seq)}\n")
                else:
                    outfile.write(f"{seq}\n")
                    outfile.write(f"+\n{'#' * len(seq)}\n")

    return outfiles


###################
# create bowtie index
###################
def bowtie_index(working_dir):
    # define arguments to create index
    reads_fasta = os.path.join(working_dir, 'reads.fa')
    index_name = os.path.join(working_dir, 'reads')

    # cline = f'''conda run -n {conda_env_name} bowtie-build {reads_fasta} {index_name}'''
    cline = f'''/home/flavia/Documents/Concordia/Project/Nanopore/download/bowtie2-2.4.5-linux-x86_64/bowtie2-build {reads_fasta} {index_name}'''
    os.system(cline)

    return index_name


###################
# align barcodes to reads
###################
def bowtie(args):
    # read arguments
    barcodes_fastq_1, barcodes_fastq_2, working_dir = args

    # build index
    reads_fasta_index = os.path.join(working_dir, 'reads')
    if not os.path.isfile(reads_fasta_index + ".1.ebwt"):
        print("Index does not exist")

    # run alignment for 5' barcodes
    out_sam_file_1 = os.path.join(working_dir, "alignment_5.sam")
    os.system(f"bowtie2 -x {reads_fasta_index} -q {barcodes_fastq_1} -v 1 -a --sam > {out_sam_file_1}")

    # run alignment for 3' barcodes
    out_sam_file_2 = os.path.join(working_dir, "alignment_3.sam")
    os.system(f"bowtie2 -x {reads_fasta_index} -q {barcodes_fastq_2} -v 1 -a --sam > {out_sam_file_2}")

    return [out_sam_file_1, out_sam_file_2]


###################
# assign barcodes to reads and create a dataframe
###################
def separate_reads(sam_file):
    # column names of sam file
    col_names = ['barcode', 'flag', 'fastq', 'start_position', 'map_quality', 'cigar_string',
                 'next_read_name', 'position_next_read', 'template_length', 'fastq_sequence', 'barcode_quality',
                 'unknown_1', 'unknown_2', 'unknown_3', 'unknown_4']

    # load sam file as df
    try:
        map_df = pd.read_csv(sam_file, sep='\t', comment="@", header=None)
    except:
        map_df = pd.DataFrame(
            [line for line in list(map(lambda x: x.strip("\n").split("\t"), open(sam_file).readlines())) if
             len(line) == 15])

    map_df.columns = col_names
    map_df.index = map_df['fastq'].to_numpy()

    # load sequence lengths dictionary
    seq_lengths_dict = {}

    # add lengths of sequences
    for line in open(sam_file):
        if line.startswith("@SQ"):
            seq_lengths_dict[line.split()[1][3:]] = int(line.split()[2][3:])

    # reads that have multiple alignments have more than 1 line
    multiple_alignments_fastq = np.unique(map_df.loc[:, 'fastq'].to_numpy(), return_counts=True)
    multiple_alignments_fastq = multiple_alignments_fastq[0][multiple_alignments_fastq[1] >= 2]

    # subset df for assignments
    map_df = map_df.loc[np.invert(map_df['fastq'].duplicated().to_numpy())]

    # iteratively identify most likely barcode
    for seq in multiple_alignments_fastq:
        # the alignment with the barcode closer to the end of the read is selected
        barcode = np.argmin(np.minimum(int(map_df.loc[seq, 'start_position']),
                                       np.abs(int(map_df.loc[seq, 'start_position']) - int(seq_lengths_dict[seq]))))

        # assign barcode to the fastq
        map_df.loc[seq, 'barcode'] = barcode

    # save dataframe to file
    # outfile = os.path.join(working_directory, "barcode_assignments.tsv")
    outfile = sam_file + "-barcode_assignment.tsv"
    map_df.to_csv(outfile, sep='\t')

    # seq_lengths_dict = pd.read_table(fastq_length_file, header=None)
    # seq_lengths_dict = dict(zip(list(map(lambda x: x[3:], seq_lengths_dict.iloc[:,1].to_numpy(dtype='str'))), list(map(lambda x: int(re.sub("[^0-9]", "", x)), seq_lengths_dict.iloc[:,2].to_numpy()))))

    return outfile


def concat_dfs(combined_fasta, master_dir, barcode_version):
    # iteratively concatenate dfs
    for i, file in enumerate(combined_fasta):
        if i == 0:
            res_df = pd.read_csv(file, sep='\t', index_col=0)
        else:
            res_df = pd.concat([res_df, pd.read_csv(file, sep='\t', index_col=0)])

    # # write to file
    # res_file = os.path.join(master_dir, (str(barcode_version) + ".tsv"))
    # res_df.to_csv(res_file, sep='\t')

    return res_df


###################
# create directories for final assignments
###################
def create_final_directories(barcode_df, master_dir):
    # three parent directories: barcode_5_assigned_only, barcode_3_assigned_only, fully_assigned
    parent_dir = os.path.join(master_dir, 'read_assignments')

    # make parent directory
    try:
        os.makedirs(parent_dir)
    except:
        print(f"Making working directory failed for parent directory {parent_dir}")

    # rename the "0" assignment to no assignment
    barcode_df = barcode_df.replace("0", "no_assignment")
    barcode_df = barcode_df.replace(0, "no_assignment")

    # for each barcode_5, make a directory in fully_assigned and barcode_5_assigned_only
    for barcode in np.unique(barcode_df.loc[:, 'barcode_5'].to_numpy(dtype='str')):
        # make directory
        dir = os.path.join(parent_dir, barcode)
        try:
            os.makedirs(dir)
        except:
            print(f"Making working directory failed for directory {dir}")

    return barcode_df, parent_dir


###################
# move reads to respective folders
###################
def organize_fastq_reads(combined_fastq, final_parent_dir, barcode_df):
    # iteratively read lines of combined gzipped file
    with gzip.open(combined_fastq, "rb") as infile:
        for i, line in enumerate(infile):
            if i % 4 == 0:
                # first " " delimination is the read name minus the 
                read_name = line.decode("utf-8").strip("@").split()[0]

                # assign outfile 
                try:
                    oufile = os.path.join(final_parent_dir, barcode_df.loc[read_name, 'barcode_5'],
                                          (barcode_df.loc[read_name, 'barcode_3'] + ".fastq.gz"))
                except:
                    oufile = os.path.join(final_parent_dir, "no_assignment", "no_assignment.fastq.gz")

            # write read to outfile
            with gzip.open(oufile, "ab") as outfile:
                outfile.write(line)

    return i


###################
# cleanup results directory
###################
def cleanup_results_dir(master_dir, final_parent_dir):
    # list all files to remove
    files = np.array(glob.glob(os.path.join(master_dir, "*")))

    # exclude the final parent directory
    files = files[np.invert(files == final_parent_dir)]
    for file in files:
        try:
            shutil.rmtree(file)
        except:
            os.remove(file)


def main():
    print('demultiplex.py')
    t1 = time.time()

    # try making the output directory
    try:
        os.mkdir(master_dir)
    except:
        pass

    ###################
    # turn fastqs into fasta
    ###################
    # turn reads fastqs into fastas
    combined_fasta, combined_fastq = fastq_to_fasta(input_directory, master_dir)

    ###################
    # separate fastas into smaller fastas
    ###################
    # arguments for fasta into smaller fasta for faster processing
    individual_directories = split_fasta(combined_fasta, master_dir)

    ###################
    # read barcodes from metadata and turn into fastqs in each sub directory
    ###################
    # file for meta_data
    metadata_file = os.path.join(input_directory, "meta_data.tsv")

    # barcode fastqs
    barcode_fastqs = barcodes_from_metadata(metadata_file, master_dir)

    ###################
    # create index from each reads.fa
    ###################
    args = [(working_dir) for working_dir in individual_directories]

    with ProcessPoolExecutor(max_workers=threads) as executor:
        # directories for parallelized processing
        _null = list(executor.map(bowtie_index, individual_directories))

    ###################
    # align barcodes to reads
    ###################
    # arguments for mapping 5' barcode
    args = [(barcode_fastqs[0], barcode_fastqs[1], working_dir) for working_dir in individual_directories]

    with ProcessPoolExecutor(max_workers=threads) as executor:
        # directories for parallelized processing
        sam_files = list(executor.map(bowtie, args))

    ###################
    # assign barcodes to reads and create a dataframe describing read assignments
    ###################
    args = list(map(lambda x: x[0], sam_files))

    # process sam files for barcode assignment
    with ProcessPoolExecutor(max_workers=threads) as executor:
        # directories for parallelized processing
        barcodes_5_tsv = list(executor.map(separate_reads, args))

    args = list(map(lambda x: x[1], sam_files))

    # process sam files for barcode assignment
    with ProcessPoolExecutor(max_workers=threads) as executor:
        # directories for parallelized processing
        barcodes_3_tsv = list(executor.map(separate_reads, args))

    # concatenate files
    barcode_5_df = concat_dfs(barcodes_5_tsv, master_dir, "barcode_5")
    barcode_3_df = concat_dfs(barcodes_3_tsv, master_dir, "barcode_3")

    # concatenate data frames
    barcode_df = barcode_5_df.loc[:, ['barcode']].copy()
    barcode_df.columns = ['barcode_5']

    barcode_3_df = barcode_3_df.loc[:, ['barcode']]
    barcode_3_df.columns = ['barcode_3']

    barcode_df = pd.concat([barcode_df, barcode_3_df])
    barcode_df = barcode_df.fillna('no_assignment')

    ###################
    # create directories for final assignments
    ###################
    barcode_df, final_parent_dir = create_final_directories(barcode_df, master_dir)

    ###################
    # move reads into correct directories and files
    ###################
    reads_assigned = organize_fastq_reads(combined_fastq, final_parent_dir, barcode_df)
    print(f"Number of reads assigned: {reads_assigned}")

    ###################
    # save dataframe to output
    ###################
    outfile = os.path.join(master_dir, "read_assignments.tsv")
    barcode_df.to_csv(outfile, sep='\t')

    ###################
    # cleanup results directory
    ###################
    cleanup_results_dir(master_dir, final_parent_dir)

    print(time.time() - t1)
    print('demultiplex.py COMPLETE')


if __name__ == "__main__":
    # number of threads for multi-threaded processes
    threads = 1

    # name of the conda env
    conda_env_name = "venv"

    # the necessary input files and described output file
    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:b:c:d:",
                                   ["input_directory=", "master_dir=", "threads=", "conda_env_name="])
    except getopt.GetoptError:
        print(help_message)
        sys.exit(2)

    for opt, arg in opts:
        # assumes all files in the input directory are gzipped fastq's
        # assumes the file describing barcodes is called "meta_data.tsv" and is in the input_directory
        if opt in ("-a", "--input_directory"):
            input_directory = str(arg)
        if opt in ("-b", "--master_dir"):
            master_dir = str(arg)
        if opt in ("-c", "--threads"):
            threads = int(arg)
        if opt in ("-d", "--conda_env_name"):
            conda_env_name = str(arg)

    main()