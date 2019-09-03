import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import gzip
import time
import os
import sys
import logging
from difflib import SequenceMatcher
from functools import reduce
import subprocess
import shlex
import argparse
from datetime import datetime

SAMPLE_NAME_KEY = "sample_name"
FILETYPE_KEY = "filetype"
FILENAME_KEY = "filename"
LIBRARY_LAYOUT_KEY = "library_layout"
FASTQ_TYPE_KEY = "fastq"
BAM_TYPE_KEY = "bam"
VCF_TYPE_KEY = "vcf"
LIB_PAIRED_KEY = "Paired"
FILESIZE_KEY = "filesize"
POLLING_PERIOD = 1  # in seconds
MAX_NUM_PROCESS = 10
COMPLETE = 1
WAITING = -1
RUNNING = 0


BAM_FILE_END_OF_BLOCK = b'\x1f\x8b\x08\x04\x00\x00\x00' \
                        b'\x00\x00\xff\x06\x00BC' \
                        b'\x02\x00\x1b\x00\x03\x00\x00' \
                        b'\x00\x00\x00\x00\x00\x00\x00'


########################################################################################################
# Manifest level functions
########################################################################################################


def collect_sample_level_info(df_manifest):
    info = {}
    group = df_manifest.groupby(SAMPLE_NAME_KEY)
    for sample_name in group.groups.keys():
        df_sample = group.get_group(sample_name)
        sample_info = dict()
        sample_info["common_filename_prefix"] = get_common_filename_prefix(df_sample[FILENAME_KEY])
        sample_info["sample_name_in_filename_prefix"] = (sample_name in sample_info["common_filename_prefix"])
        sample_info.update(get_sample_filemodel(df_sample[FILETYPE_KEY]))
        filetype_group = df_sample.groupby(FILETYPE_KEY)
        for filetype in filetype_group.groups.keys():
            df_filetype = filetype_group.get_group(filetype)
            sample_info[f"{filetype}_extension"] = get_common_file_extension(df_filetype[FILENAME_KEY])
        info[sample_name] = sample_info
    df_result = pd.DataFrame.from_dict(info, orient="index")
    df_result.index.name = SAMPLE_NAME_KEY
    return df_result


def add_filesize(df_manifest, datadirpath):
    df_manifest[FILESIZE_KEY] = [
        os.path.getsize(os.path.join(datadirpath, filename)) for filename in df_manifest[FILENAME_KEY]
    ]
    return df_manifest


def check_fastq_pairing(df_manifest):
    sample_group = df_manifest.groupby(by=[SAMPLE_NAME_KEY, FILETYPE_KEY])
    samples = set(df_manifest[SAMPLE_NAME_KEY])
    check_pass = True
    for sample_name in samples:
        df_sample_fastq = sample_group.get_group((sample_name, FASTQ_TYPE_KEY))
        if df_sample_fastq.shape[0] % 2 != 0:
            logging.warning(f"WARNING: Number of FastQ files for sample {sample_name} is not a multiple of 2.")
        paired_fastq = df_sample_fastq.groupby(LIBRARY_LAYOUT_KEY).get_group(LIB_PAIRED_KEY).loc[:, FILENAME_KEY]
        for fastq in paired_fastq:
            has_pair = False
            for possible_pair in paired_fastq:
                if _is_fastq_pair(fastq, possible_pair):
                    has_pair = True
            if not has_pair:
                logging.warning(f"ERROR: Missing partner FastQ file for paired fastq file {fastq} "
                                f"in sample {sample_name}.")
                check_pass = False
    return check_pass


########################################################################################################
# Worker Functions
########################################################################################################


def get_common_filename_prefix(filelist):
    # ensure that the shortest string is matched first.
    prefixes = sorted([filename[: filename.index('.')] for filename in filelist], key=len)
    return reduce(lambda prefix1, prefix2: longest_common_substring(prefix1, prefix2), prefixes)


def get_common_file_extension(filelist):
    extensions = map(lambda filename: filename[filename.index('.'):], filelist)
    return reduce(lambda ext1, ext2: longest_common_substring(ext1, ext2), extensions)


def get_sample_filemodel(sample_filetype_list):
    model = {}
    for ftype in sample_filetype_list:
        if ftype not in model:
            model[ftype] = 0
        model[ftype] += 1
    return model


def get_missing_files(filelist, datadirpath):
    missing_files = []
    for filename in filelist:
        filepath = os.path.join(datadirpath, filename)
        if not os.path.isfile(filepath):
            missing_files.append(filename)
    return missing_files


def longest_common_substring(string1, string2):
    matcher = SequenceMatcher(isjunk=None, a=string1, b=string2)
    match = matcher.find_longest_match(0, len(string1), 0, len(string2))
    return string1[match.a: match.a + match.size]


def check_all_files_exist(filelist, datadirpath):
    return len(get_missing_files(filelist, datadirpath)) == 0


def check_samplename_in_filename(sample_name, sample_file_list):
    for filename in sample_file_list:
        if sample_name not in filename:
            return False
    return True


def check_filetype_in_filename(filetype, file_list):
    for filename in file_list:
        if filetype not in filename:
            return False
    return True


def check_bam_file_eof(filename):
    logging.info(f"Checking BAM EOF for `{filename}`")
    n = len(BAM_FILE_END_OF_BLOCK)
    with open(filename, 'rb') as infile:
        infile.seek(0, os.SEEK_END)
        infile.seek(infile.tell() - n, 0)
        byte_string = infile.read(n)
    return byte_string == BAM_FILE_END_OF_BLOCK


def generate_sample_bam_vcf_pairs(df_manifest):
    group = df_manifest.groupby(by=[SAMPLE_NAME_KEY, FILETYPE_KEY])
    for sample_name in sorted(df_manifest.loc[:, SAMPLE_NAME_KEY]):
        df_bam = group.get_group((sample_name, BAM_TYPE_KEY))
        df_vcf = group.get_group((sample_name, VCF_TYPE_KEY))
        n_bam, n_vcf = df_bam.shape[0], df_vcf.shape[0]
        if n_bam > 1 or n_vcf > 1:
            logging.warning(f"WARNING: Unexpected number of BAM ({n_bam}) and VCF ({n_vcf}).")
        yield (sample_name, df_bam.loc[:, FILENAME_KEY].iloc[0], df_vcf.loc[:, FILENAME_KEY].iloc[0])


def get_sample_vcfpath(df_manifest, sample_name):
    df_sample_vcf = df_manifest.groupby(by=(SAMPLE_NAME_KEY, FILETYPE_KEY)).get_group((sample_name, VCF_TYPE_KEY))
    if df_sample_vcf.shape[0] < 1:
        logging.error(f"ERROR: Sample {sample_name} has no VCF file.")
        return
    if df_sample_vcf.shape[0] != 1:
        logging.warning(f"WARNING: Sample {sample_name} has more than one VCF files.")
    return list(df_sample_vcf[FILENAME_KEY])[0]


def get_fastq_file_pair(filenames):
    pairs = []
    num_file = len(filenames)
    for i in range(num_file - 1):
        for j in range(i + 1, num_file):
            name1 = filenames[i]
            name2 = filenames[j]
            if _is_fastq_pair(name1, name2):
                pairs.append((name1, name2))
    return pairs


def _is_fastq_pair(name1, name2):
    if len(name1) != len(name2):
        return False
    name1 = name1.replace('R1', '')
    name1 = name1.replace('R2', '')
    name2 = name2.replace('R1', '')
    name2 = name2.replace('R2', '')
    return name1 == name2


########################################################################################################
# Utilities
########################################################################################################


def is_gzip(filepath: str) -> bool:
    """
    Check if a the file specified by filepath
    is a gzipped file. We do this by checking
    whether the first two bytes is the magic
    numbers included in the file header.
    We check them bytewise to avoid issue with
    endianness.

    Note that this wouldn't detect if the gzip
    file is a concatenation of multiple "members".
    See gzip specification at:
     https://www.ietf.org/rfc/rfc1952.txt
    """
    if not os.path.isfile(filepath):
        logging.error("ERROR: The file %s does not exist" % filepath)
        return False

    with open(filepath, 'rb') as filehandle:
        byte1 = filehandle.read(1)
        byte2 = filehandle.read(1)
    if byte1 == b'\x1f' and byte2 == b'\x8b':
        return True
    else:
        return False


def open_handler(file_path, filemode='rt'):
    if is_gzip(file_path):
        return gzip.open(file_path, filemode)
    else:
        return open(file_path, filemode)

########################################################################################################
# Multi-process handling
########################################################################################################


def cmd_runner(cmd, stdin=None, stdout=None):
    args = shlex.split(cmd)
    return subprocess.Popen(args, stdin=stdin, stdout=stdout)


def dumb_scheduler(param_list, max_process_num=10, polling_period=0.5):
    num_jobs = len(param_list)
    logging.info(f"Scheduler starting with {num_jobs} jobs.")
    remaining_job_index = list(range(num_jobs - 1, -1, -1))
    job_states = np.full(num_jobs, WAITING, dtype=np.int8)
    processes = []
    while not np.all(job_states == COMPLETE):
        while np.sum(job_states == RUNNING) < max_process_num and WAITING in job_states:
            job_index = remaining_job_index.pop()
            cmd, stdin, stdout = param_list[job_index]
            logging.info(f"Running cmd: {cmd}")
            processes.append((job_index, cmd_runner(cmd, stdin=stdin, stdout=stdout)))
            job_states[job_index] = RUNNING

        # Poll current active processes
        for job_index, proc in processes:
            if job_states[job_index] != COMPLETE:
                poll_state = proc.poll()
                if poll_state is not None:
                    return_code = poll_state
                    cmd, stdin, stdout = param_list[job_index]
                    if return_code != 0:
                        logging.warning(f"WARNING: Command ```{cmd}``` return with nonzero return code: {return_code}.")
                    job_states[job_index] = COMPLETE
                    logging.info(f"Command ```{cmd}``` completed with return code {return_code}.")
                    logging.info(f"Number of remaining jobs : {np.sum(job_states != COMPLETE)}")
                    logging.info(f"Number of completed jobs : {np.sum(job_states == COMPLETE)}")
                    logging.info(f"Number of running jobs   : {np.sum(job_states == RUNNING)}")
        time.sleep(polling_period)
    for cmd, stdin, stdout in param_list:
        if stdin is not None:
            stdin.close()
        if stdout is not None:
            stdout.close()
    logging.info("Exiting job scheduler.")
    return processes


########################################################################################################
# Main(s)
########################################################################################################


def start_log(log="stderr", level=logging.DEBUG):
    """
    Initiate program logging. If no log file is specified,
    then log output goes to the default log file: stdout
    """
    if log.lower() == "stdout":
        log = sys.stdout
    elif log.lower() == "stderr":
        log = sys.stderr
    else:
        log = open(log, 'w')
    logging.basicConfig(stream=log,
                        level=level,
                        filemode='w',
                        format="%(asctime)s %(message)s",
                        datefmt="[%m/%d/%Y %H:%M:%S] ")
    logging.info("Program started")
    logging.info("Command line: {0}\n".format(' '.join(sys.argv)))
    return


def parse_args():
    parser = argparse.ArgumentParser(
        description="Bioinformatics file checks."
    )
    parser.add_argument('--manifest', metavar="PATH", type=str, required=True,
                        help='Path to manifest file.')
    parser.add_argument('--datadir', metavar="PATH", type=str, required=True,
                        help='Path to data directory')
    parser.add_argument('--outdir', metavar='PATH', type=str, required=True,
                        help="Name of output directory.")
    parser.add_argument('--logfile', metavar="FILENAME", type=str, required=False,
                        default="stderr",
                        help="Log file. Default: stderr")
    parser.add_argument('--max_num_process', metavar="NUM", type=int, required=False,
                        default=MAX_NUM_PROCESS,
                        help=f"Maximum number of subprocess spawned. Default: {MAX_NUM_PROCESS}")
    parser.add_argument('--polling_period', metavar="SECONDS", type=float, required=False,
                        default=POLLING_PERIOD,
                        help=f"Polling period. Default {POLLING_PERIOD}")
    return parser.parse_args()


def main():
    start_time = datetime.now()
    user_options = parse_args()
    max_num_process = user_options.max_num_process
    polling_period = user_options.polling_period
    manifestpath = user_options.manifest
    datadirpath = user_options.datadir
    outdirpath = user_options.outdir
    logfile = user_options.logfile

    if not os.path.isdir(outdirpath):
        os.mkdir(outdirpath)
    if logfile.lower() not in ["stderr", "stdout"]:
        logfile = os.path.join(outdirpath, logfile)
    start_log(log=logfile)

    # Manifest file checks
    df_manifest = pd.read_csv(manifestpath, sep='\t')
    sample_group = df_manifest.groupby(SAMPLE_NAME_KEY)
    missing_files = get_missing_files(df_manifest[FILENAME_KEY], datadirpath)
    if missing_files:
        logging.error("FATAL: These files are missing: {0}.".format('\n'.join(missing_files)))
        sys.exit(1)
    if not check_fastq_pairing(df_manifest):
        logging.error("FATAL: Fastq paring error. ")
        sys.exit(1)

    # Fastq pair check
    fastq_files = list(df_manifest.groupby(FILETYPE_KEY).get_group(FASTQ_TYPE_KEY).loc[:, FILENAME_KEY])
    fastq_pairs = get_fastq_file_pair(fastq_files)
    fastq_pair_check_cmd = """fastq_pair_check {filepath1} {filepath2}"""
    fastq_pair_check_params = []
    for name1, name2 in fastq_pairs:
        print(name1, name2)
        filepath1 = os.path.join(datadirpath, name1)
        filepath2 = os.path.join(datadirpath, name2)
        cmd = fastq_pair_check_cmd.format(filepath1=filepath1, filepath2=filepath2)
        fastq_pair_check_params.append((cmd, None, None))
    dumb_scheduler(fastq_pair_check_params, max_process_num=os.cpu_count(), polling_period=5)

    bam_files = list(df_manifest.groupby(FILETYPE_KEY).get_group(BAM_TYPE_KEY).loc[:, FILENAME_KEY])
    for name in bam_files:
        filepath = os.path.join(datadirpath, name)
        check_result = check_bam_file_eof(filepath)
        if not check_result:
            logging.error(f"ERROR: BAM file EOF check failed for `{name}`.")

    df_sample_info = collect_sample_level_info(df_manifest)
    df_manifest = add_filesize(df_manifest, datadirpath)

    img_dir = os.path.join(outdirpath, "images")
    if not os.path.isdir(img_dir):
        os.mkdir(img_dir)

    df_manifest.to_csv(os.path.join(outdirpath, "manifest.tsv"), sep='\t', index=False)
    df_sample_info.to_csv(os.path.join(outdirpath, "sample_info.tsv"), sep='\t')
    filetype_group = df_manifest.groupby(FILETYPE_KEY)
    n_filetype_group = len(filetype_group.groups)
    fig, ax = plt.subplots(n_filetype_group, 1, figsize=(10, 5 * n_filetype_group))
    c = 0
    for filetype in filetype_group.groups.keys():
        df = filetype_group.get_group(filetype)
        df = df.set_index(SAMPLE_NAME_KEY)
        df.loc[:, FILESIZE_KEY].plot(kind='bar', ax=ax[c])
        ax[c].set_xticklabels(ax[c].get_xticklabels(), rotation=90)
        ax[c].set_title(f"File sizes: {filetype}")
        c += 1
    fig.savefig(os.path.join(img_dir, "filesize_by_filetype.png"))

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    df = df_sample_info.reset_index().loc[:, [SAMPLE_NAME_KEY, FASTQ_TYPE_KEY, BAM_TYPE_KEY, VCF_TYPE_KEY]]
    sns.barplot(x=SAMPLE_NAME_KEY, y="value", hue="variable", data=pd.melt(df, id_vars=SAMPLE_NAME_KEY), ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    fig.savefig(os.path.join(img_dir, "filemodel_by_sample.png"))

    # Run external commands
    fastqc_cmd = """fastqc {filepath} --outdir {outdirpath} """
    fastq_validator_cmd = ("fastQValidator --noeof "
                           "--printableErrors 100 "
                           "--auto "
                           "--file {filepath}")
    samtools_stats_cmd = "samtools stats --remove-dups --required-flag 1 {filepath}"
    validatesamfile_cmd = "java -jar /jarfiles/picard.jar ValidateSamFile I={filepath} MODE=VERBOSE"
    bcftools_stats_cmd = "bcftools stats {filepath}"
    samtools_index_cmd = "samtools index -@ 4 {filepath} {samtools_index_outfile}"
    verify_reference_cmd = "bash verify_reference.sh {bamfilepath} {vcffilepath}"
    verifybamid_cmd = ("verifyBamID "
                       "--bam {bamfilepath} "
                       "--vcf {vcffilepath} "
                       "--bai {bamindexfilepath} "
                       "--precise "
                       "--best "
                       "--chip-mix "
                       "--ignoreRG "
                       "--noPhoneHome "
                       "--out {outprefix}")

    processes_params = []
    for sample_name in sorted(sample_group.groups.keys()):
        df_sample = sample_group.get_group(sample_name)
        for i, row in df_sample.iterrows():
            filetype = row[FILETYPE_KEY]
            filepath = os.path.join(datadirpath, row[FILENAME_KEY])
            basename = os.path.basename(filepath).split('.')[0]
            # Run fastqc
            if filetype in [FASTQ_TYPE_KEY, BAM_TYPE_KEY]:
                cmd = fastqc_cmd.format(filepath=filepath, outdirpath=outdirpath)
                processes_params.append((cmd, None, None))
            # Run FastqValidator
            if filetype == FASTQ_TYPE_KEY:
                cmd = fastq_validator_cmd.format(filepath=filepath)
                stdout = open(os.path.join(outdirpath, f"{basename}_FastQValidator.txt"), 'w')
                processes_params.append((cmd, None, stdout))
            # Run `samtools stats`
            if filetype == BAM_TYPE_KEY:
                cmd = samtools_stats_cmd.format(filepath=filepath)
                stdout = open(os.path.join(outdirpath, f"{basename}_samtools_stats.txt"), 'w')
                processes_params.append((cmd, None, stdout))
            # Run ValidateSamFile
            if filetype == BAM_TYPE_KEY:
                cmd = validatesamfile_cmd.format(filepath=filepath)
                stdout = open(os.path.join(outdirpath, f"{basename}_ValidateSamFile.txt"), 'w')
                processes_params.append((cmd, None, stdout))
            # Run `bcftools stats`
            if filetype == VCF_TYPE_KEY:
                cmd = bcftools_stats_cmd.format(filepath=filepath)
                stdout = open(os.path.join(outdirpath, f"{basename}_bcftools_stats.txt"), 'w')
                processes_params.append((cmd, None, stdout))
    dumb_scheduler(processes_params, max_process_num=max_num_process, polling_period=polling_period)

    # Run verify_reference.sh
    for sample_name, bamfilename, vcffilename in generate_sample_bam_vcf_pairs(df_manifest):
        bamfilepath = os.path.join(datadirpath, bamfilename)
        vcffilepath = os.path.join(datadirpath, vcffilename)
        args = shlex.split(verify_reference_cmd.format(bamfilepath=bamfilepath, vcffilepath=vcffilepath))
        return_code = subprocess.call(args)
        if return_code != 0:
            logging.warning(f"WARNING: BAM and VCF reference assembly for sample {sample_name} does not match.")

    # Run `samtools index` before verifybamid
    samtools_index_params = []
    df_bam = df_manifest.groupby(FILETYPE_KEY).get_group(BAM_TYPE_KEY)
    bamfilelist = []
    bamindexfilelist = []
    for i, row in df_bam.iterrows():
        filepath = os.path.join(datadirpath, row[FILENAME_KEY])
        samtools_index_outfile = os.path.join(outdirpath, f"{os.path.basename(filepath)}.bai")
        bamfilelist.append(filepath)
        bamindexfilelist.append(samtools_index_outfile)
        cmd = samtools_index_cmd.format(filepath=filepath, samtools_index_outfile=samtools_index_outfile)
        stdout = None
        samtools_index_params.append((cmd, None, stdout))
    dumb_scheduler(samtools_index_params, max_process_num=max_num_process, polling_period=polling_period)

    verifybamid_params = []
    for bamfilepath, bamindexfilepath in zip(bamfilelist, bamindexfilelist):
        bam_name = os.path.basename(bamfilepath)
        sample_name = df_manifest.loc[df_manifest[FILENAME_KEY] == bam_name, SAMPLE_NAME_KEY]
        vcf_name = df_manifest.groupby(by=[FILETYPE_KEY, SAMPLE_NAME_KEY]).get_group((VCF_TYPE_KEY, sample_name))
        vcffilepath = ' '.join([os.path.join(datadirpath, name) for name in vcf_name])
        outprefix = os.path.join(outdirpath, f"{sample_name}_VerifyBamID_output")
        cmd = verifybamid_cmd.format(bamfilepath=bamfilepath,
                                     vcffilepath=vcffilepath,
                                     bamindexfilepath=bamindexfilepath,
                                     outprefix=outprefix)
        verifybamid_params.append((cmd, None, None))
    dumb_scheduler(verifybamid_params, max_process_num=max_num_process, polling_period=polling_period)

    # Running MultiQC
    date = datetime.now().strftime("%d-%m-%Y-%H:%M")
    multiqc_cmd = f"""multiqc --force \
                              --dirs-depth 1 \
                              --fullnames \
                              --title "MULTIQC {date}" \
                              --module samtools \
                              --module bcftools \
                              --module picard \
                              --module fastqc \
                              --module verifybamid \
                              --outdir {outdirpath} \
                              {outdirpath}"""
    subprocess.call(shlex.split(multiqc_cmd))
    time_taken = (datetime.now() - start_time).total_seconds()
    logging.info(f"Total time taken: {time_taken // 3600}:{ time_taken // 60}:{time_taken % 60}")
    return


if __name__ == "__main__":
    main()

