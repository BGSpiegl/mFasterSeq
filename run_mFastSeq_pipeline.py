#!/usr/bin/python3

# ----------------------------------------------------------------------------------------------------------------------
# Python imports
from pathlib import Path
from datetime import datetime
from natsort import humansorted
from re import match as re_match
from sys import exit as sys_exit
from re import DOTALL as re_DOTALL
from argparse import ArgumentParser, Namespace
from typing import List, Dict, Tuple
from csv import reader as csv_reader
from multiprocessing import Pool as mproc_Pool
from shutil import which as sh_which, move as sh_move
from os.path import join as pth_join, isfile as pth_isfile
from os import walk as os_walk, chdir as os_chdir, makedirs as os_mkdrs
from os import environ as os_environ, pathsep as os_pathsep
from subprocess import run as subp_run, Popen as subp_Popen, check_output as subp_check_out, PIPE as subp_PIPE, \
    STDOUT as subp_STDOUT
from os.path import basename as pth_basename, isdir as pth_isdir, split as pth_split, sep as pth_sep, \
    abspath as pth_abspath, dirname as pth_dirname
# ----------------------------------------------------------------------------------------------------------------------
# version
MAJOR_RELEASE = 2
MINOR_RELEASE = 1
PATCH_NUMBER = 2
PIPE_VERSION = f'v{MAJOR_RELEASE}.{MINOR_RELEASE}.{PATCH_NUMBER}'
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# functions
def get_commandline_args() -> Namespace:
    # WARNING: newlines ('\n') in help will not be output to the command line!
    cmdln_prsr = ArgumentParser(add_help=True)
    # versioning
    cmdln_prsr.add_argument('--version', action='version', version=PIPE_VERSION)
    # required arguments
    cmdln_prsr.add_argument('-i', '--input-directory', dest='in_dir', required=True, type=Path, metavar='File',
                            help='Directory containing your sample FastQs (sequencing data). Use absolute paths!')
    # optional arguments
    cmdln_prsr.add_argument('-r', '--reference-genome', dest='ref_genome_fasta', type=Path, metavar='File',
                            help='Path to your local indexed reference genome file. If the index ' +
                                 'files are not found, the pipeline will attempt to create them using BWA, samtools ' +
                                 'and picard. This argument is required if you run the mFastSeq pipeline from ' +
                                 'the beginning (including read alignment)')
    cmdln_prsr.add_argument('-b', '--bam-input', dest='input_is_bam', action='store_true',
                            help='Flag: if set, the input directory will be searched for BAM files ' +
                                 'instead of FastQ files. Use this to skip the alignment step if your ' +
                                 "samples have already been aligned. Only searches for BAM files with name matching " +
                                 "'(?P<NAME>).*'.sorted.bam. If another name pattern should be used, please specify "
                                 "via the '--alternative-fastq-name-pattern' parameter. "
                                 "Should also recognize SAM files.")
    cmdln_prsr.add_argument('-x', '--sexs', dest='sex_file', type=Path, metavar='File',
                            help="Path to comma-separated text file (*.csv) containing the biological sexs of your "
                                 "samples. Sex can be left unspecified. One sample per line: '<sample ID>," +
                                 "[optional sex specification: either 'f' or 'm']'. A valid example for such a file "
                                 "would be: P1_S18,m<next line>P2_S7,<next line>P3_S9,f<end of file>" +
                                 "If no sex code (either m or f) is given, the sex will be guessed based " +
                                 "on the X-chromosomal portion of all reads. The name you specify must follow the " +
                                 "name pattern used in sample ID extraction. For instance, the default naming " +
                                 "pattern would take 'W-P791' from the file " +
                                 "'W-P791_S2_R1.fastq.gz' which would be the Illumina sample name up to '_S?'.")
    cmdln_prsr.add_argument('-afc', '--alternative-female-control', dest='alternative_female_control',
                            type=Path, metavar='File',
                            help="If you don't want to use the old initial controls on which the " +
                                 "old plasmaDB and pipeline were based on, provide an alternative " +
                                 "FEMALE control here. Per default, the reference files created on te 11th of " +
                                 "January 2022 will be used (i.a. '2022-01-11_FastSeq_controls_female.txt')\n")
    cmdln_prsr.add_argument('-amc', '--alternative-male-control', dest='alternative_male_control',
                            type=Path, metavar='File',
                            help="If you don't want to use the old initial controls on which the " +
                                 "old plasmaDB and pipeline were based on, provide an alternative " +
                                 "MALE control here.\nPer default, the reference files created on te 11th of " +
                                 "January 2022 will be used (i.a. '2022-01-11_FastSeq_controls_male.txt')\n")
    cmdln_prsr.add_argument('-p', '--parallel-processes', dest='n_parallel', default=4, type=int,
                            help="Number of processes that should be used in computation. " +
                                 "The more the better. Do not use more than you have! [default: 4]")
    cmdln_prsr.add_argument('-o', '--output-directory', dest='out_dir', type=Path, metavar='File',
                            help="OPTIONAL: if provided, pipeline output will be placed in this directory. " +
                                 "Otherwise, output will be placed in parent directory of input FASTQ directory " +
                                 "under '<DATE>_mFastSeq_output'. The first <DATE> part follows the international " +
                                 "standard date and time notation: https://www.cl.cam.ac.uk/~mgk25/iso-time.html. " +
                                 "Existing files will be overwritten!")
    # plotting arguments
    cmdln_prsr.add_argument('-vp', '--vertical-plot', dest='create_vertical_plot', action='store_true',
                            help='Flag: if this flag is set, also a vertical plot of the chromosome arm profiles ' +
                                 'will be created.')
    # arguments for troubleshooting
    cmdln_prsr.add_argument('-ape', '--add-paths-to-environment', action='append', dest='add_paths',
                            help="TROUBLESHOOTING: if you get an error message stating an executable for " +
                                 "samtools, picard tools, or bwa was not found, you probably " +
                                 "must add a path to one of your local directories containing those to " +
                                 "the 'pythonpath' so that the pipeline can find them. Provide these here " +
                                 "(multiple paths possible by re-declaring parameter like this: " +
                                 "'-ape <path_1> -ape <path_2>' etc.).")
    cmdln_prsr.add_argument('-afnp', '--alternative-fastq-name-pattern', dest='alternative_r1_fastq_pattern',
                            help="TROUBLESHOOTING: if the pipeline tells you that no FastQs were found, " +
                                 "check the fastq folder input you provided. If the provided path is correct, " +
                                 "your FastQ files might be named differently than the default search pattern " +
                                 "expects: r'(?P<NAME>.*)_S.+_R1.*.fastq.gz' (default Illumina file naming). " +
                                 "In this case, you must rename your files or provide an alternative RegEx pattern " +
                                 "via this parameter. You might need to double-escape backslashes. " +
                                 "WARNING: MUST contain an Identifier for the 'R1' read if R1 and R2 read fastq " +
                                 "files are present in input directory! Can be used to specify BAM file name pattern.")
    cmdln_prsr.add_argument('-sd', '--search-subdirectories', dest='search_subdirs', action='store_true',
                            help="Flag: if set, the input directory's subdirectories will also be searched for " +
                                 "fastQ files matching the (default/provided) name pattern.")
    cmdln_prsr.add_argument('-k', '--keep-sam', action='store_true', dest='keep_sam',
                            help='Do not delete the SAM file!')
    return cmdln_prsr.parse_args()


def get_folders(parent_directory: Path) -> list:
    return [dr[0] for dr in os_walk(parent_directory)]  # returns full paths!


def get_custom_files(search_dir: Path, name_pattern: str, search_subdirectories=True) -> List[Dict[str, str]]:
    # supports auto annotations
    cust_files = list()
    output_data = list()
    if '(?P<NAME>' not in name_pattern:
        raise AttributeError("name_patterns must include '(?P<NAME>'")
    name_pattern = name_pattern.strip('"').strip("'")
    if search_subdirectories:
        sub_dirs = get_folders(parent_directory=search_dir)
        for sub_d in sub_dirs:
            sub_d = Path(sub_d)
            if not sub_d:
                continue
            search_result = list(sub_d.glob("*." + name_pattern.split('.')[-1]))
            if search_result:
                cust_files.extend([sub_d / res for res in humansorted(search_result)])
    else:
        search_result = list(search_dir.glob("*." + name_pattern.split('.')[-1]))
        cust_files.extend([search_dir / res for res in humansorted(search_result)])
    for idx, cust_file in enumerate(cust_files):
        try:
            name_dict = re_match(name_pattern, cust_file.name, re_DOTALL).groupdict()
            if len(name_dict) > 1:
                annotations = list(name_dict.keys())
                del annotations[annotations.index('NAME')]
                out_dict = {'name': name_dict['NAME'], 'file': cust_file}
                for anno in annotations:
                    out_dict.update({anno.lower(): name_dict[anno]})
                output_data.append(out_dict)
            else:
                output_data.append({'name': name_dict['NAME'], 'file': cust_file})
        except AttributeError:
            pass  # exclude sample
    if len(output_data) == 0:
        raise IOError("provided name pattern yielded 0 files. " +
                      "Check name pattern and if you would actually want to search in subdirectories!")
    return output_data


def parallel_fast_seq_worker(data_tuple: Tuple):
    sn, ref_genome, fq1, out_dir = data_tuple
    subp_run(['bash', str(fast_seq_script_path), sn, str(ref_genome), str(fq1), str(out_dir)])
    return pth_join(out_dir, sn + '.chrarm.txt')


def parallel_fast_seq_bam_worker(data_tuple: Tuple):
    sn, bam, chrom_arm_file, out_dir = data_tuple
    python3_path = sh_which('python3')
    print(f"    creating results table for {sn} ..")
    cmnd = [python3_path, str(fast_seq_read_counting_script_path), bam, out_dir / f'{sn}.chrarm.txt',
            chrom_arm_file]  # data_tuple: (fq1_data['name'], ref_genome_fasta, fq1_data['file'], out_dir_sex)
    # receiving: aln_file = sys_argv[1]
    #     results_file = sys_argv[2]
    #     chr_arms_file = sys_argv[3]
    _ran_read_counting = subp_run(cmnd)
    return pth_join(out_dir, sn + '.chrarm.txt')


def parallel_zscores_worker(data_tuple: tuple):
    sex = None
    try:
        out_dir, sn, chrm_arm_file, sex, exclude_arms_file = data_tuple
    except ValueError:
        out_dir, sn, chrm_arm_file, exclude_arms_file = data_tuple
    if sex is None or sex == 'Sex_guessed':
        if USE_ALTERNATIVE_CONTROLS and pth_isfile(alternative_controls['Male']) and \
                pth_isfile(alternative_controls['Female']):
            subp_run(['python3', str(z_scores_fastseq_script), '-ca', chrm_arm_file, '-gs',
                      '-o', f'{pth_join(out_dir, sn)}.zscoresANDreadcounts.txt',
                      '--female-control', alternative_controls['Female'],
                      '--male-control', alternative_controls['Male'],
                      '--exclude-chromosome-arms-csv', str(exclude_arms_file)])
        else:
            subp_run(['python3', str(z_scores_fastseq_script), '-ca', chrm_arm_file, '-gs',
                      '-o', f'{pth_join(out_dir, sn)}.zscoresANDreadcounts.txt',
                      '--female-control', controls['Female'], '--male-control', controls['Male'],
                      '--exclude-chromosome-arms-csv', str(exclude_arms_file)])
    else:
        if USE_ALTERNATIVE_CONTROLS and pth_isfile(alternative_controls['Male']) and \
                pth_isfile(alternative_controls['Female']):
            subp_run(['python3', str(z_scores_fastseq_script), '-ca', chrm_arm_file, '-s', sex,
                      '-o', f'{pth_join(out_dir, sn)}.zscoresANDreadcounts.txt',
                      '--female-control', alternative_controls['Female'],
                      '--male-control', alternative_controls['Male'],
                      '--exclude-chromosome-arms-csv', str(exclude_arms_file)])
        else:
            subp_run(['python3', str(z_scores_fastseq_script), '-ca', chrm_arm_file, '-s', sex,
                      '-o', f'{pth_join(out_dir, sn)}.zscoresANDreadcounts.txt',
                      '--female-control', controls['Female'], '--male-control', controls['Male'],
                      '--exclude-chromosome-arms-csv', str(exclude_arms_file)])


def ensure_supporting_ref_genome_files_existence(ref_fasta: Path):
    print(f" -| checking supporting files for reference genome ..", flush=True)
    BWA_PATH = sh_which('bwa')
    SAMTOOLS_PATH = sh_which('samtools')
    PICARD_PATH = sh_which('picard')
    dct = ref_fasta.parent / f'{ref_fasta.stem}.dict'
    fai = ref_fasta.parent / f'{ref_fasta}.fai'
    bwt = ref_fasta.parent / f'{ref_fasta}.bwt'
    if not dct.is_file():
        if PICARD_PATH is None:
            raise EnvironmentError(f"Could not determine path to your picardtools executable when " +
                                   f"trying to create sequence dictionary '{dct}' of reference " +
                                   f"genome file '{ref_fasta}' because it was not found.")
        print(" creating PicardTools sequence dictionary ..", flush=True)
        subp_run([PICARD_PATH, "CreateSequenceDictionary", f'REFERENCE="{ref_fasta}"', f'OUTPUT="{dct}"'])
    if not fai.is_file():
        if SAMTOOLS_PATH is None:
            raise EnvironmentError(f"Could not determine path to your samtools executable when " +
                                   f"trying to create fasta index '{fai}' of reference " +
                                   f"genome file '{ref_fasta}' because it was not found.")
        print(" -| creating samtools index ..", flush=True)
        subp_run([SAMTOOLS_PATH, "faidx", ref_fasta])
    if not bwt.is_file():
        if BWA_PATH is None:
            raise EnvironmentError(f"Could not determine path to your bwa executable when " +
                                   f"trying to create bwa mapping index '{bwt}' of reference " +
                                   f"genome file '{ref_fasta}' because it was not found.")
        print(" -| creating bwa index ..", flush=True)
        subp_run([BWA_PATH, "index", "-a", "bwtsw", ref_fasta])


def document_tools(out_folder: Path) -> None:
    BWA_PATH = sh_which('bwa')
    SAMTOOLS_PATH = sh_which('samtools')
    # PYTHON2_PATH = sh_which('python2')
    out_folder.mkdir(parents=True, exist_ok=True)
    with open(out_folder / 'tools_doc.txt', 'wt') as f_tooldoc:
        # no version getting for bwa -> only error output
        bwa_proc = subp_Popen(BWA_PATH, stderr=subp_PIPE, stdout=subp_PIPE)
        f_tooldoc.write(bwa_proc.communicate()[1].decode('utf-8'))
        f_tooldoc.write("-" * 70 + '\n')
        f_tooldoc.write(subp_check_out(f'{SAMTOOLS_PATH} --version',
                                       shell=True, stderr=subp_STDOUT).decode('utf-8'))
        f_tooldoc.write("-" * 70 + '\n')
        # f_tooldoc.write(subp_check_out(f'{PYTHON2_PATH} --version',
        #                                shell=True, stderr=subp_STDOUT).decode('utf-8'))
        f_tooldoc.write(subp_check_out(f'python3 --version',
                                       shell=True, stderr=subp_STDOUT).decode('utf-8'))


def read_sex_data(sex_file: str) -> dict:
    # sex is here either 'm' or 'f'
    sex_dat = dict()
    with open(sex_file, 'rt') as f_sex:
        for line in f_sex:
            try:
                pid, sex = line.strip().split(',')
            except ValueError:
                pid = line.strip().split(',')[0]
                sex = None
            if sex_dat.get(pid):  # insanity check: sex input
                if sex_dat[pid] != sex:
                    raise AssertionError(f"Inconsistent multiple sex entries for ID {pid}!")
            else:
                sex_dat.update({pid: sex})
    return sex_dat


def does_match_sex(cur_sex: str, sex_dat: dict, cur_fq1_data: dict) -> bool:
    # sex_dat: sex_dat.update({pid: gend})
    # cur_fq1_data: output_data.append({'name': name_dict['NAME'], 'file': cust_file})
    #               name_dict depends on name pattern r'(?P<NAME>.*)_S.+_R1.fastq.gz' and
    #               MUST match what was provided in the sex file -> pid
    if cur_sex == 'Sex_guessed':
        if sex_dat is None:
            return True
        if sex_dat[cur_fq1_data['name']] is None:
            return True
        else:
            return False
    cur_sample_sex = remap_sex(sex_dat[cur_fq1_data['name']])
    return cur_sample_sex == cur_sex


def remap_sex(sex: str):
    # if sex is None:
    #     return None
    sex_map = {'Male': 'm',
               'Female': 'f',
               'Sex_guessed': None,
               'm': 'Male',
               'f': 'Female',
               None: 'Sex_guessed'}
    return sex_map[sex]


def plot_fast_seq_parallel(data_tuple: dict):
    replace_na = None  # NA values will be replaced with None -> no contribution in computation
    # [1] "/<SOME-PATH>/lib/R/bin/exec/R"
    # [2] "--slave"
    # [3] "--no-restore"
    # [4] "--file=test_cmd_args.R"
    # [5] "--args"
    # [6] "some_file"
    out_dir.mkdir(parents=True, exist_ok=True)
    file_to_plot = data_tuple['file']
    filled_template = out_dir / f'{file_to_plot.stem}.fast_seq.csv'
    chrarm_coordinates = dict()
    # read in template lines
    with open(plotting_template_path, 'rt') as f_temp:
        chrarm_coordinates.update({'header': f_temp.readline()})  # read header line
        for line in f_temp:
            chrm, arm, strt, stp = line.strip().split('\t')[:-1]
            chrarm_coordinates.update({arm: f'{chrm}\t{arm}\t{strt}\t{stp}\t'})
    # fill template with sample Z-scores
    filed_template_buffer = list()
    with open(file_to_plot, 'rt') as f_zscores:
        with open(filled_template, 'wt') as f_filled:
            filed_template_buffer.append(chrarm_coordinates['header'])  # write header line
            for line in f_zscores:
                # key+"\t"+str(sample_reads[key])+"\t"+str(sample_fract[key])+"\t"+str(sample_z_score[key])+"\n"
                arm, _n_reads, _ratio, z_score = line.strip().split('\t')
                if line[:10] != 'genomewide':  # skip genomewide z-score for plotting
                    if z_score == 'n/a':
                        if replace_na:
                            filed_template_buffer.append(chrarm_coordinates[arm] + f'{replace_na}\n')
                        else:
                            continue  # DON'T WRITE CHROMOSOME ARM IF NO READS MAPPED TO IT! (e.g. WALDO analysis)
                    else:
                        filed_template_buffer.append(chrarm_coordinates[arm] + f'{z_score}\n')
            f_filled.writelines(filed_template_buffer)
    # plot!
    subp_run(['Rscript', str(plot_script_path_horizontal_chrlen), str(filled_template)])
    if create_vertical_plot:
        subp_run(['Rscript', str(plot_script_path_vertical_fastseq), str(filled_template)])


def mv_files(files_data: list, out_dir: str):
    for files_dat in files_data:
        mv_file = files_dat['file']
        sh_move(mv_file, pth_join(out_dir, pth_basename(mv_file)))


def rn_files(files_data: list, postfix=None):
    for files_dat in files_data:
        rn_file = files_dat['file']
        out_file = pth_join(pth_split(rn_file)[0],
                            f"{pth_basename(rn_file).split('.')[0]}{postfix if postfix is not None else ''}.png")
        sh_move(rn_file, out_file)


def convert_sams_parallel(dir: Path, n_parallel: int, keep_sams: bool):
    sam_data = get_custom_files(search_dir=dir, name_pattern=r'(?P<NAME>.*).sam', search_subdirectories=False)
    samtools_path = sh_which('samtools')
    if len(sam_data) == 1:
        sam_converter(sam_data[0]['file'], samtools_path=samtools_path, keep_sam=keep_sams)
    else:
        parallel_args = list()
        for sam in [dat['file'] for dat in sam_data]:
            parallel_args.append((sam, samtools_path, keep_sams))
        with mproc_Pool(n_parallel) as converter_pool:
            converter_pool.starmap(sam_converter, parallel_args)


def sam_converter(sam: Path, samtools_path: str, keep_sam: bool):
    out_bam = f'{sam.stem}.sorted.bam'
    sh_command_convert_sort = f"{samtools_path} view -hb {sam} -o - | {samtools_path} sort --threads 4 " + \
                              f"-o {out_bam} - "
    conv_proc = subp_run(sh_command_convert_sort, shell=True)
    conv_proc.check_returncode()
    index_proc = subp_run([samtools_path, 'index', '-b', str(out_bam), f'{out_bam}.bai'])
    index_proc.check_returncode()
    if not keep_sam:
        sam.unlink()


if __name__ == '__main__':
    # record time ------------------------------------------------------------------------------------------------------
    pipeline_start = datetime.now()
    print(f"\nStarted mFastSeq Pipeline (version {PIPE_VERSION}) on {str(pipeline_start).split(' ')[0]} at " +
          f"{str(pipeline_start).split(' ')[1].split('.')[0]}", flush=True)
    # get commandline arguments and unpack them ------------------------------------------------------------------------
    cmd_args = get_commandline_args()
    input_is_bam = cmd_args.input_is_bam
    keep_sam = cmd_args.keep_sam
    create_vertical_plot = cmd_args.create_vertical_plot
    search_subdirs = cmd_args.search_subdirs
    n_parallel = cmd_args.n_parallel
    in_dir = cmd_args.in_dir
    out_dir = cmd_args.out_dir
    sex_file = cmd_args.sex_file
    ref_genome_fasta = cmd_args.ref_genome_fasta
    add_paths = cmd_args.add_paths
    alternative_r1_fastq_pattern = cmd_args.alternative_r1_fastq_pattern
    alternative_female_control = cmd_args.alternative_female_control
    alternative_male_control = cmd_args.alternative_male_control
    if not input_is_bam and (ref_genome_fasta is None or ref_genome_fasta == ''):
        print("Error: no reference FastA file specified. Please use the '--reference-genome' parameter and try again.")
        sys_exit(1)
    r1_fastq_name_pattern = r'(?P<NAME>.*)_S.+_R1.*.fastq.gz'
    chrm_arm_name_pattern = r'(?P<NAME>.*).chrarm.txt'
    fast_seq_zscores_name_pattern = r'(?P<NAME>.*).zscoresANDreadcounts.txt'  # for plotting: use file with all columns!
    # set and check script paths ---------------------------------------------------------------------------------------
    script_root_path = Path(__file__).parent
    if n_parallel >= 4:
        fast_seq_script_path = script_root_path / 'src/fast_pipeline_adapted_4_threads.sh'
    else:
        fast_seq_script_path = script_root_path / 'src/fast_pipeline_adapted_single_thread.sh'
    if not fast_seq_script_path.is_file():
        raise FileNotFoundError(f"FastSeq script not found under path '{fast_seq_script_path}'")
    fast_seq_read_counting_script_path = script_root_path / 'src/faster_seq_analysis_chr_arm.py'
    if not fast_seq_read_counting_script_path.is_file():
        raise FileNotFoundError(f"Read counting script not found under path '{fast_seq_read_counting_script_path}'")
    chrom_arm_file_hg19 = script_root_path / 'ref/chrom_arm_hg19.txt'
    if not chrom_arm_file_hg19.is_file():
        raise FileNotFoundError(f"Chromosome arms file not found under path '{chrom_arm_file_hg19}'")
    z_scores_fastseq_script = script_root_path / 'src/z_scores.py'
    if not z_scores_fastseq_script.is_file():
        raise FileNotFoundError(f"FastSeq Z-scores script not found under path '{z_scores_fastseq_script}'")
    plot_script_path_horizontal_chrlen = script_root_path / 'src/makeGraph_linear_chrlen_fastseq_adapted.R'
    if not plot_script_path_horizontal_chrlen.is_file():
        raise FileNotFoundError(f"Horizontal plots script not found under path '{plot_script_path_horizontal_chrlen}'")
    plot_script_path_vertical_fastseq = script_root_path / 'src/makeGraph_vertical_fastseq_adapted.R'
    if not plot_script_path_vertical_fastseq.is_file():
        raise FileNotFoundError(f"Vertical plots script not found under path '{plot_script_path_vertical_fastseq}'")
    plotting_template_path = script_root_path / 'src/Template.fast_seq.csv'
    if not plotting_template_path.is_file():
        raise FileNotFoundError(f"Plotting template file not found under path '{plotting_template_path}'")
    controls = {'Female': script_root_path / 'ref/2022-01-12_FastSeq_controls_female.txt',
                'Male': script_root_path / 'ref/2022-01-12_FastSeq_controls_male.txt'}
    if not controls['Female'].is_file():
        raise FileNotFoundError(f"Female controls file not found under path '{controls['Female']}'")
    if not controls['Male'].is_file():
        raise FileNotFoundError(f"Male controls file not found under path '{controls['Male']}'")
    exclude_arms_path = script_root_path / 'ref/exclude_chromosome_arms.csv'
    if not exclude_arms_path.is_file():
        raise FileNotFoundError(f"Scaffold arms exclusion file not found under path '{exclude_arms_path}'")
    with open(exclude_arms_path, 'rt') as f_excl:
        EXCLUDE_LOW_MAPABILITY_ARMS = list(csv_reader(f_excl))[0]  # scaffold names as CSVs (single row)
    # process commandline arguments ------------------------------------------------------------------------------------
    if alternative_r1_fastq_pattern is not None and alternative_r1_fastq_pattern != '':
        r1_fastq_name_pattern = alternative_r1_fastq_pattern
    if add_paths is not None and add_paths != '':
        if not isinstance(add_paths, list):
            add_paths = list(add_paths)
        os_environ["PATH"] += os_pathsep + os_pathsep.join(add_paths)
    if out_dir is None or not out_dir:
        out_dir = in_dir.parent / f'{datetime.now().isoformat()[:10]}_mFastSeq_output'
    alternative_controls = dict()
    USE_ALTERNATIVE_CONTROLS = False
    if alternative_female_control is not None and alternative_female_control != '':
        if not pth_isfile(alternative_female_control):
            raise AttributeError(f'alternative female controls file given but could not access: ' +
                                 f'{alternative_female_control}')
        alternative_controls.update({'Female': alternative_female_control})
        USE_ALTERNATIVE_CONTROLS = True
    if alternative_male_control is not None and alternative_male_control != '':
        if not pth_isfile(alternative_male_control):
            raise AttributeError(f'alternative male controls file given but could not access: ' +
                                 f'{alternative_male_control}')
        alternative_controls.update({'Male': alternative_male_control})
        USE_ALTERNATIVE_CONTROLS = True
    # document software versions ---------------------------------------------------------------------------------------
    document_tools(out_folder=out_dir)  # creates output directory
    # give user feedback about parameters ------------------------------------------------------------------------------
    # check if bwa etc. can be found in the current runtime environment
    tool_paths = (sh_which('bwa'), sh_which('samtools'), sh_which('picard'))
    if any([tlpth is None for tlpth in tool_paths]):
        tls_not_found_str = "could not find path of these tools: " + ' '.join(
            [tl if tl is not None else '' for tl in tool_paths])
        raise EnvironmentError(tls_not_found_str)
    # ensure the presence of BWA index, dict and bwt for reference genome FastA ----------------------------------------
    input_bam_data = []
    if not input_is_bam:
        ensure_supporting_ref_genome_files_existence(ref_fasta=ref_genome_fasta)
        r1_fastq_data = get_custom_files(search_dir=in_dir, name_pattern=r1_fastq_name_pattern,
                                         search_subdirectories=search_subdirs)
    else:
        r1_fastq_data = []
        aln_name_patterns = ("(?P<NAME>.*).sam", "(?P<NAME>.*).bam")
        if alternative_r1_fastq_pattern is not None and alternative_r1_fastq_pattern != '':
            aln_name_patterns = (alternative_r1_fastq_pattern, )
        for aln_name_pattern in aln_name_patterns:
            try:
                input_bam_data.extend(get_custom_files(search_dir=in_dir, name_pattern=aln_name_pattern,
                                                       search_subdirectories=search_subdirs))
                for dat_dict in input_bam_data:
                    while '.sorted' in dat_dict['name']:
                        dat_dict['name'] = dat_dict['name'].replace('.sorted', '')
            except IOError:  # nothing found
                continue
        # check if something was found:
        if not input_bam_data:
            raise FileNotFoundError(f"no input aligment data could be found in input directory '{in_dir}'")
    # process files per sex -----------------------------------------------------------------------------------------
    if n_parallel >= 4:
        fastseq_pool = mproc_Pool(int(n_parallel / 4))
    else:
        fastseq_pool = mproc_Pool(int(n_parallel / 1))
    # read sex mapping
    if sex_file is not None and sex_file != '' and not input_is_bam:
        if pth_isfile(sex_file):
            sex_data = read_sex_data(sex_file=sex_file)
        else:
            print(f" !  WARNING: sample sex file not found at '{sex_file}'. Will guess sex instead", flush=True)
            sex_data = None
    else:
        sex_data = None
    # create required sex output directories
    if input_is_bam:  # sex will always be 'Sex_guessed' for BAM input
        out_dirs = {'Sex_guessed': out_dir / 'Sex_guessed'}
        out_dirs['Sex_guessed'].mkdir(parents=True, exist_ok=True)
    else:
        if sex_data is None:
            sex_iter = (None,)
            all_sex_specs = ('Sex_guessed',)
        else:
            all_sex_specs = [remap_sex(sex=sx) for sx in list(sex_data.values())]
            sex_iter = set(list(sex_data.values()))
        out_dirs = dict.fromkeys(all_sex_specs)
        for sex in sex_iter:
            remapped_sex = remap_sex(sex=sex)
            if out_dirs[remapped_sex] is None:  # remapped_sex out dir not defined
                out_dir_sex = out_dir / remapped_sex
                out_dir_sex.mkdir(parents=True, exist_ok=True)
                out_dirs[remapped_sex] = out_dir_sex
    # process samples per sex ---------------------------------------------------------------------------------------
    for sex in out_dirs.keys():  # sex will always be 'Sex_guessed' for BAM input
        # compute chromosome arm read counts ---------------------------------------------------------------------------
        print(f" -| computing chromosome arm-wide read counts for biological sex '{sex}'..", flush=True)
        parallel_args = list()
        out_dir_sex = out_dirs[sex]
        os_chdir(out_dir_sex)
        if not input_is_bam:  # default - input is FastQ
            for fq1_data in r1_fastq_data:
                if does_match_sex(cur_sex=sex, sex_dat=sex_data, cur_fq1_data=fq1_data):
                    parallel_args.append((fq1_data['name'], ref_genome_fasta, fq1_data['file'], out_dir_sex))
            if not parallel_args:
                print(f" -| INFO: No samples matching sex '{sex}'. Continuing ..", flush=True)
                try:
                    out_dirs[sex].unlink()  # only works if the directory is really empty -> keep older results!
                except Exception:
                    pass
                continue
            fastseq_pool.map(parallel_fast_seq_worker, parallel_args)
        else:  # start only z-scores computation if input is (are) BAM file(s)
            print(" !  WARNING: operating directly on SAM files was specified. Input SAM/BAM files from " +
                  f"{in_dir} MUST come from FastSeq samples, otherwise computation results will be wrong.")
            for bam_data in input_bam_data:
                parallel_args.append((bam_data['name'], bam_data['file'], str(chrom_arm_file_hg19), out_dir_sex))
            if not parallel_args:
                print(f" -| INFO: No samples matching sex '{sex}'. Continuing ..", flush=True)
                try:
                    out_dirs[sex].unlink()  # only works if the directory is really empty -> keep older results!
                except:
                    pass
                continue
            fastseq_pool.map(parallel_fast_seq_bam_worker, parallel_args)
        # compute z-scores ---------------------------------------------------------------------------------------------
        chrom_arm_data = get_custom_files(out_dir_sex, name_pattern=chrm_arm_name_pattern,
                                          search_subdirectories=False)
        zscores_parallel_args = list()
        for chrm_dat in chrom_arm_data:
            sn = chrm_dat['name']
            chrm_fl = chrm_dat['file']
            if sex_data is not None and sex_data[sn] != remap_sex(sex):
                print(f" !  THIS IS A SANITY WARNING (your bioinformatician might be a little crazy): " +
                      f"inconsistent sex information for sample {sn}: {sex_data[sn]} (user specified) and " +
                      f"{sex} (internal program logic). Sex will be guessed.", flush=True)
                sex = None
            sex_for_zscores = remap_sex(sex=sex)  # just use sex from the outer loop above
            zscores_parallel_args.append((out_dir_sex, sn, chrm_fl, sex_for_zscores, exclude_arms_path))
        # convert SAMs to BAMs, sort, index
        if not input_is_bam or any([inp_aln['file'].name.lower().endswith('sam') for inp_aln in input_bam_data]):
            print(" -| converting SAM files to BAM, sorting, and creating index file ..", flush=True)
            convert_sams_parallel(dir=out_dir_sex, keep_sams=keep_sam, n_parallel=n_parallel)
        fastseq_pool.map(parallel_zscores_worker, zscores_parallel_args)
        # remove chromosome arm count files (information is present in the 'zscoresANDreadcounts' file)
        for ch_arm_dat in chrom_arm_data:
            try:
                ch_arm_dat['file'].unlink()
            except (FileNotFoundError, IsADirectoryError):
                continue
        # create chrarm profile plots ----------------------------------------------------------------------------------
        print(" -| creating chromosome arm copy number profile plots ..", flush=True)
        plots_out_dir = out_dir_sex
        analysis_files_data = get_custom_files(search_dir=out_dir_sex, name_pattern=fast_seq_zscores_name_pattern,
                                               search_subdirectories=False)
        n_parallel_plot = n_parallel
        if n_parallel_plot > 8:
            n_parallel_plot = 6  # restrict the number of processes for plotting
        fastseq_plot_pool = mproc_Pool(n_parallel_plot)
        fastseq_plot_pool.map(plot_fast_seq_parallel, analysis_files_data)
        # move filled templates
        print(" -| moving output to designated directory ..", flush=True)
        out_csv_data = get_custom_files(search_dir=out_dir, name_pattern=r'(?P<NAME>.*).fast_seq.csv',
                                        search_subdirectories=False)
        filled_templates_dir = pth_join(out_dir, 'filled_templates')
        if not pth_isdir(filled_templates_dir):
            os_mkdrs(filled_templates_dir)
        mv_files(files_data=out_csv_data, out_dir=filled_templates_dir)
        # rename pngs
        out_png_data = get_custom_files(search_dir=out_dir,
                                        name_pattern=r'(?P<NAME>.*).zscoresANDreadcounts.fast_seq.csv.' +
                                                     'chrscale_fast-seq.png',
                                        search_subdirectories=False)
        rn_files(files_data=out_png_data, postfix='.chrscale_fast-seq')
    # give temporal user feedback --------------------------------------------------------------------------------------
    pipeline_stop = datetime.now()
    print(f"Finished mFastSeq pipeline (version {PIPE_VERSION}) on {str(pipeline_stop).split(' ')[0]} at " +
          f"{str(pipeline_stop).split(' ')[1].split('.')[0]} after " +
          f"{str(pipeline_stop - pipeline_start).split('.')[0]} (hh:mm:ss).\n", flush=True)
