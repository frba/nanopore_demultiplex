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
def fastq_to_fasta(input_directory, output_directory):
    # list the fastq files in the directory
    fastqs = [file for file in glob.glob(os.path.join(input_directory, '*')) if bool(re.search("fastq|fq", file))]

    # define combined fastq
    combined_fastq = os.path.join(output_directory, "all_passed_fastqs.fastq.gz")

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

    # split fastq file into smaller fastqs for future processing reads into places
    regex = os.path.join(output_directory, "temporary_fastq")
    cline = f"zcat {combined_fastq} | split -l 100000 - {regex}"
    os.system(cline)

    # list plit fastqs
    fastqs_out = glob.glob(regex + "*")

    return fasta_out, fastqs_out


###################
# separate fastas into smaller fastas
###################
def split_fasta(infile, output_directory):
    # max file size is 100Mb (seems to work well on most systems for speed)
    max_bytes = 100000000

    # number of splits for the
    splits = math.ceil(os.path.getsize(infile) / max_bytes)

    # create directories for each file
    dirs = list(map(lambda x: os.path.join(output_directory, f"working_{x}"), range(splits)))
    for dir in dirs:
        try:
            os.makedirs(f"{dir}")
        except:
            print(f"Making working directory failed for {dir}")

    # total bytes already read for all files, and the total bytes in the file we are processing
    total_bytes = 0
    file_number = 0

    # initial output file
    outfile = open(os.path.join(output_directory, f"working_{file_number}", (f"reads.fa")), "w")

    # separate fasta into parts by memory size
    for i, line in enumerate(open(infile)):
        if i % 2 == 0 and total_bytes > max_bytes:
            # close the file and open the new one
            outfile.close()
            file_number += 1
            total_bytes = 0
            outfile = open(os.path.join(output_directory, f"working_{file_number}", (f"reads.fa")), "w")

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
        # for i, barcode in enumerate(np.unique(ref_df.loc[:, spec_type].to_numpy())):
        #     to_add = f">barcode_{1 + round(len(res_fasta) / 2)}"
        #     res_fasta.extend([to_add, barcode])

        #Code change to add the barcode in same order as the input file
        for i, barcode in ref_df[spec_type].drop_duplicates().iteritems():
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
# read templates from metadata and turn into fastqs in each sub directory
###################
def templates_from_metadata(metadata_file, master_dir):
    # load the dataframe
    ref_df = pd.read_csv(metadata_file, sep='\t')
    ref_df.columns = list(map(lambda x: re.sub("'", "", x), ref_df.columns.to_numpy()))
    ref_df.columns = list(map(lambda x: re.sub(" ", "_", x), ref_df.columns.to_numpy()))

    res_fasta = []
    df_template = ref_df.drop_duplicates(subset=['Full_amplicon_SEQUENCE'])
    #Code change to add the barcode in same order as the input file
    outpath = os.path.join(master_dir, 'templates', 'fasta')
    try:
        os.makedirs(outpath)
    except:
        pass

    for i, row in df_template.iterrows():
        template_name = f">{row['Amplicon_name']}"
        res_fasta.extend([template_name, row['Full_amplicon_SEQUENCE'].upper()])

        with open(os.path.join(outpath, f"{row['Amplicon_name']}.fa"), "w") as outfile:
            outfile.write("\n".join(res_fasta))


###################
# create bowtie index
###################
def bowtie_index(working_dir):
    # define arguments to create index
    reads_fasta = os.path.join(working_dir, 'reads.fa')
    index_name = os.path.join(working_dir, 'reads')
    # try:
    #     cline = f'''conda run -n {conda_env_name} bowtie-build {reads_fasta} {index_name}'''
    #     os.system(cline)
    # except:
    cline = f'''/home/flavia/Documents/Concordia/Project/Nanopore/download/bowtie-1.3.1-linux-x86_64/bowtie-build {reads_fasta} {index_name}'''
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

    # try:
    #     # run alignment for 5' barcodes
    #     out_sam_file_1 = os.path.join(working_dir, "alignment_5.sam")
    #     os.system(f"conda run -n {conda_env_name} bowtie -x {reads_fasta_index} -q {barcodes_fastq_1} -v 1 -a --sam > {out_sam_file_1}")
    #     # os.system(f"bowtie -x {reads_fasta_index} -q {barcodes_fastq_1} -v 1 -a --sam > {out_sam_file_1}")
    #
    #     # run alignment for 3' barcodes
    #     out_sam_file_2 = os.path.join(working_dir, "alignment_3.sam")
    #     os.system(f"conda run -n {conda_env_name} bowtie -x {reads_fasta_index} -q {barcodes_fastq_2} -v 1 -a --sam > {out_sam_file_2}")
    #     # os.system(f"bowtie -x {reads_fasta_index} -q {barcodes_fastq_2} -v 1 -a --sam > {out_sam_file_2}")
    # except:
    # run alignment for 5' barcodes
    out_sam_file_1 = os.path.join(working_dir, "alignment_5.sam")
    os.system(f"/home/flavia/Documents/Concordia/Project/Nanopore/download/bowtie-1.3.1-linux-x86_64/bowtie -x {reads_fasta_index} -q {barcodes_fastq_1} -v 1 -a --sam > {out_sam_file_1}")

    # run alignment for 3' barcodes
    out_sam_file_2 = os.path.join(working_dir, "alignment_3.sam")
    os.system(f"/home/flavia/Documents/Concordia/Project/Nanopore/download/bowtie-1.3.1-linux-x86_64/bowtie -x {reads_fasta_index} -q {barcodes_fastq_2} -v 1 -a --sam > {out_sam_file_2}")

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
             len(line) == 15]) #????????

    map_df.columns = col_names
    map_df.index = map_df['fastq'].to_numpy()
    # Removed lines with no sequences
    map_df = map_df[map_df['flag'] != 4]

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
def create_final_directories(barcode_5, output_directory):
    # # three parent directories: barcode_5_assigned_only, barcode_3_assigned_only, fully_assigned
    # parent_dir = os.path.join(output_directory, 'read_assignments')

    # # make parent directory
    # try:
    #     os.makedirs(parent_dir)
    # except:
    #     print(f"Making working directory failed for parent directory {parent_dir}")

    # for each barcode_5, make a directory in fully_assigned and barcode_5_assigned_only
    for barcode in barcode_5:
        # make directory
        dir = os.path.join(output_directory, barcode)
        try:
            os.makedirs(dir)
        except:
            pass
            # print(f"Making working directory failed for directory {dir}")

    return output_directory


###################
# move reads to respective folders
###################
def organize_fastq_reads(args):
    infile_fastq, output_directory, barcode_df, iteration = args

    # create working directory for this file
    working_dir = os.path.join(output_directory, ("read_assignment_" + str(iteration)))
    try:
        os.mkdir(working_dir)
    except:
        print(f"Making dir failed for {working_dir}")

    # create directory tree
    working_dir = create_final_directories(barcode_df.loc[:, 'barcode_5'].to_numpy(dtype='str'), working_dir)

    # iteratively read lines of combined gzipped file
    with open(infile_fastq) as infile:
        for i, line in enumerate(infile):
            if i % 4 == 0:
                # first " " delimination is the read name minus the
                read_name = line.strip("@").split()[0]

                # assign outfile
                try:
                    oufile = os.path.join(working_dir, barcode_df.loc[read_name, 'barcode_5'],
                                          (barcode_df.loc[read_name, 'barcode_3'] + ".fastq.gz"))
                except:
                    oufile = os.path.join(working_dir, "no_assignment", "no_assignment.fastq.gz")

            # write read to outfile
            with gzip.open(oufile, "ab") as outfile:
                outfile.write(line.encode('utf-8'))

    return working_dir


def combine_output_files(args):
    # decompose args
    regex, dir_5, file_3, final_parent_dir = args

    # go through all possible dir and file combinations
    files_of_interest = glob.glob(regex)

    if len(files_of_interest) > 0:
        # define final output file
        final_output = os.path.join(final_parent_dir, dir_5, (file_3 + ".fastq.gz"))

        # iteratively concatenate file to final output
        for file in files_of_interest:
            cline = f"cat {file} >> {final_output}"
            os.system(cline)

    return


###################
# count assigned reads
###################
def counts_read_assignemnts(barcode_df):
    count_series = barcode_df.groupby(['barcode_5', 'barcode_3']).size()
    count_df = count_series.to_frame(name='count').reset_index()

    return count_df


###################
# cleanup results directory
###################
def cleanup_results_dir(output_directory, mask_arr):
    # list all files to remove
    files = np.array(glob.glob(os.path.join(output_directory, "*")))

    # exclude the final parent directory
    files = files[np.invert(np.isin(files, mask_arr))]
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
        os.mkdir(output_directory)
    except:
        pass

    ###################
    # turn fastqs into fasta
    ###################
    # turn reads fastqs into fastas
    combined_fasta, fastqs_out = fastq_to_fasta(input_directory, output_directory)

    ###################
    # separate fastas into smaller fastas
    ###################
    # arguments for fasta into smaller fasta for faster processing
    individual_directories = split_fasta(combined_fasta, output_directory)

    ###################
    # read barcodes from metadata and turn into fastqs in each sub directory
    ###################
    # file for meta_data
    metadata_file = os.path.join(input_directory, "meta_data.tsv")

    # barcode fastqs
    barcode_fastqs = barcodes_from_metadata(metadata_file, output_directory)
    templates_from_metadata(metadata_file, output_directory)

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
    barcode_5_df = concat_dfs(barcodes_5_tsv, output_directory, "barcode_5")
    barcode_3_df = concat_dfs(barcodes_3_tsv, output_directory, "barcode_3")

    # concatenate data frames
    barcode_df = barcode_5_df.loc[:, ['barcode']].copy()
    barcode_df.columns = ['barcode_5']

    barcode_3_df = barcode_3_df.loc[:, ['barcode']]
    barcode_3_df.columns = ['barcode_3']

    barcode_df = pd.concat([barcode_df, barcode_3_df])

    ###################
    # Change barcode_df # Fixing? Concatenating reads with same name
    ###################
    barcode_df = barcode_df.reset_index(level=0)
    barcode_df = barcode_df.groupby('index').agg({'barcode_5': 'first', 'barcode_3': 'first'}, )

    # rename the "0" assignment to no assignment
    barcode_df = barcode_df.fillna('no_assignment')
    barcode_df = barcode_df.replace("0", "no_assignment")
    barcode_df = barcode_df.replace(0, "no_assignment")

    barcode_df = barcode_df.sort_values(by=["barcode_5", "barcode_3"], ascending=True)
    ###################
    # save dataframe to output
    ###################
    outfile = os.path.join(output_directory, "read_assignments.tsv")
    barcode_df.to_csv(outfile, sep='\t')


    ###################
    # move reads into correct directories and files
    ###################
    # args = [(file, final_parent_dir, np.unique(barcode_df.loc[:,'barcode_5'].to_numpy()), np.unique(barcode_df.loc[:,'barcode_3'].to_numpy()), i) for i, file in enumerate(fastqs_out)]
    args = [(file, output_directory, barcode_df, i) for i, file in enumerate(fastqs_out)]

    # separate files into
    with ProcessPoolExecutor(max_workers=threads) as executor:
        # directories for parallelized processing
        read_assignment_dirs = list(executor.map(organize_fastq_reads, args))

    # create directories for final assignments
    final_parent_dir = create_final_directories(np.unique(barcode_df.loc[:, 'barcode_5'].to_numpy()), output_directory)

    # define directories and files
    dirs = list(np.unique(barcode_df.loc[:, 'barcode_5'].to_numpy()))
    dirs.append("no_assignment")

    files = list(np.unique(barcode_df.loc[:, 'barcode_3'].to_numpy()))
    files.append("no_assignment")

    # list regex files we need to search for
    args = []
    for dir_5 in dirs:
        for file_3 in files:
            args.append(
                [os.path.join(output_directory, "*", dir_5, (file_3 + ".fastq.gz")), dir_5, file_3, final_parent_dir])

    # combined separate files
    with ProcessPoolExecutor(max_workers=threads) as executor:
        # directories for parallelized processing
        _null = list(executor.map(combine_output_files, args))

    ###################
    # cleanup results directory
    ###################
    # a list of possible directories and the final description tsv to keep
    # dirs = list(map(lambda x: os.path.join(output_directory, x), dirs))
    # dirs.append(outfile)
    #
    # cleanup_results_dir(output_directory, np.array(dirs))

    ###################
    # count assigned reads
    ###################
    count_df = counts_read_assignemnts(barcode_df)
    out_count = os.path.join(output_directory, "count_assigned_reads.tsv")
    count_df.to_csv(out_count, sep='\t')

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
                                   ["input_directory=", "output_directory=", "threads=", "conda_env_name="])
    except getopt.GetoptError:
        print(help_message)
        sys.exit(2)

    for opt, arg in opts:
        # assumes all files in the input directory are gzipped fastq's
        # assumes the file describing barcodes is called "meta_data.tsv" and is in the input_directory
        if opt in ("-a", "--input_directory"):
            input_directory = str(arg)
        if opt in ("-b", "--output_directory"):
            output_directory = str(arg)
        if opt in ("-c", "--threads"):
            threads = int(arg)
        if opt in ("-d", "--conda_env_name"):
            conda_env_name = str(arg)

    main()