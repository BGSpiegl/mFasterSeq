#! /usr/bin/python2

# Calculate chromosome arm and genomewiede Z-scores

from argparse import ArgumentParser, Namespace
from sys import exit as sys_exit
from natsort import humansorted
from os.path import abspath as pth_abspath, dirname as pth_dirname, join as pth_join, basename as pth_basename
from csv import reader as csv_reader
from typing import Dict


def guess_biological_sex(genome, chrx_counts):
    if float(chrx_counts) / float(genome) > 0.06:  # todo: optimize this threshold!
        print("  i : sample sex guess: female", flush=True)
        this_sex = "f"
    else:
        print("  i : sample sex guess: male", flush=True)
        this_sex = "m"
    return this_sex


def get_commandline_arguments() -> Namespace:
    parser = ArgumentParser(description='Calculate FastSeq_from_nucleos1 Z-scores from Chromosome arm file')
    parser.add_argument('-ca', '--chromosome-arms-file', dest='chrarm_file',
                        help='Chromosome arm read counts', required=True)
    parser.add_argument('-eca', '--exclude-chromosome-arms-csv', dest='exclude_chromosome_arms_file',
                        help='Optional: exclude chromosome arms in this csv file. See reference counts file for '
                             'excluded chromosome arms.')
    parser.add_argument('-o', '--output-file', dest='out_file',
                        help='Z-scores output file', required=True)
    parser.add_argument('-fc', '--female-control', dest='cont_f', required=True,
                        help='Path to controls file for female samples.')
    parser.add_argument('-mc', '--male-control', dest='cont_m', required=True,
                        help='Path to controls file for male samples.')
    parser.add_argument('-s', '--biological-sex', dest='sex',
                        help='sex of sample', default="m", choices=["m", "f"])
    parser.add_argument('-gs', '--guess-sex', dest='guess_sex',
                        help='Guess sex of sample', action="store_true")
    return parser.parse_args()


if __name__ == '__main__':
    # get cmd args
    cmd_args = get_commandline_arguments()
    chrarm_file = cmd_args.chrarm_file
    exclude_chromosome_arms_file = cmd_args.exclude_chromosome_arms_file
    out_file = cmd_args.out_file
    cont_f = cmd_args.cont_f
    cont_m = cmd_args.cont_m
    sex = cmd_args.sex
    guess_sex = cmd_args.guess_sex
    # define analysis constants
    ABERRANT_ZSCORE_THRESHOLD = 5
    script_path = pth_dirname(pth_abspath(__file__))
    with open(exclude_chromosome_arms_file, 'rt') as f_excl:
        EXCLUDE_LOW_MAPPABILITY_ARMS = list(csv_reader(f_excl))[0]  # all in first line
    # i: short chromosome arms with recurrent low-read count (low mappability) by setting values in controls file to 0
    # initialize variables before updating locals with results from the commandline
    sample_name = pth_basename(chrarm_file).split('.chrarm.txt')[0]
    control_means = {}
    control_stdev = {}
    print(f"  i : working on chromosome arms file\n" +
          f"      {pth_basename(chrarm_file)}",
          flush=True)
    sample_reads = {}
    sample_fract = {}
    sample_z_score = {}
    chrarms = list()
    read_count = 0
    with open(chrarm_file, "r") as f_sample:
        for line in f_sample.readlines():
            info = line.rstrip().split("\t")
            chrm_arm = info[0]
            if chrm_arm not in EXCLUDE_LOW_MAPPABILITY_ARMS:
                reads = int(info[1])
                read_count += reads
                sample_reads[info[0]] = reads
    print(f"  i : guessing sex of sample {sample_name}\n" +
          "       (guess is always performed, but if sex was provided, guess is used to check it):", flush=True)
    # guess biological gender of sample based on fraction of reads mapping to chrX
    sex_guess = guess_biological_sex(read_count, (sample_reads["chrXp"] + sample_reads["chrXq"]))
    # define chromosome arms for both sexes
    aut_chrarms = [f'chr{i}{arm}' for i in range(1, 23) for arm in ('p', 'q')]
    male_sex_arms = [f'chr{chrm}{arm}' for chrm in ('X', 'Y') for arm in ('p', 'q')]
    female_sex_arms = [f'chrX{arm}' for arm in ('p', 'q')]
    if guess_sex:
        if sex_guess == "m":
            sex = "m"
            chrarms = aut_chrarms + male_sex_arms
        elif sex_guess == "f":
            sex = "f"
            chrarms = aut_chrarms + female_sex_arms
    elif not guess_sex:
        if sex != sex_guess:
            print(f"  ! Warning: Probably wrong sex specified for sample {sample_name}!\n" +
                  f"             You provided '{sex}' but I guessed '{sex_guess}'.", flush=True)
        else:
            print(f"  i : sex guess is in concordance with provided sex of sample {sample_name}.", flush=True)
    # correct total read count for females:
    if (guess_sex and sex_guess == 'f') or (not guess_sex and sex == 'f'):
        try:
            read_count -= (sample_reads["chrYp"] + sample_reads["chrYq"])
        except KeyError:  # in case chrY is excluded from analysis
            pass
    # choose controls based on sample sex
    controls_path = cont_f
    if sex == "m":
        controls_path = cont_m
    elif sex != "f":
        print(f" !! ERROR: undefined sex for sample {sample_name}!\nExiting ..", flush=True)
        sys_exit(1)
    try:
        with open(controls_path, 'rt') as f_controls:
            if sex == "m":
                print("  i : using male controls", flush=True)
                chrarms = aut_chrarms + male_sex_arms
            elif sex == "f":
                print("  i : using female controls")
                chrarms = aut_chrarms + female_sex_arms
            for line in f_controls.readlines():
                info = line.strip().split("\t")
                control_means[info[0]] = float(info[1])
                control_stdev[info[0]] = float(info[2])
    except FileNotFoundError:
        print(f" !! ERROR: Could not find controls file at expected location\n" +
              f"           '{controls_path}'\nExiting ..", flush=True)
        sys_exit(1)
    chrarms = humansorted(chrarms)  # achieve human-like sorting for chromosome arms
    sum_squares = 0
    for arm in chrarms:
        if arm in EXCLUDE_LOW_MAPPABILITY_ARMS:  # exclude low mappability/high variability chromosome arms
            sample_z_score[arm] = "n/a"
            continue
        sample_fract[arm] = float(sample_reads[arm]) / float(read_count)
        # i: remove short arms with recurrent low-read count by setting values in controls file to 0
        if arm not in control_means or control_means[arm] == 0:
            sample_z_score[arm] = "n/a"
        else:
            sample_z_score[arm] = (sample_fract[arm] - control_means[arm]) / control_stdev[arm]
            if abs(sample_z_score[arm]) > ABERRANT_ZSCORE_THRESHOLD:
                print(f'  i : {arm} aberrant for sample {sample_name} with Z-score: {sample_z_score[arm]}', flush=True)
            sum_squares += sample_z_score[arm] ** 2
    genomewide_z = (sum_squares - control_means["genomewide"]) / control_stdev["genomewide"]
    if genomewide_z > 5:
        print(f'  i : sample {sample_name} is genome-wide aberrant with Z-score: {genomewide_z}', flush=True)
    # create standard out file
    with open(out_file, "wt") as f_out:
        f_out.writelines([f'{arm}\t{sample_reads[arm]}\t{sample_fract[arm]}\t{sample_z_score[arm]}\n'
                          if arm in sample_reads else
                          f'{arm}\tn/a\tn/a\tn/a\n'
                          for arm in chrarms] + [f'genomewide\t{read_count}\t{sum_squares}\t{genomewide_z}'])
    # create reduced out_file
    reduced_out_file = '.'.join(out_file.split('.')[:-2] + ['zscores', 'txt'])
    with open(reduced_out_file, "wt") as f_out:
        f_out.writelines([f'{arm}\t{sample_z_score[arm]}\n' for arm in chrarms] + [f'genomewide\t{genomewide_z}'])
    # create top n aberrant chromosome arms (top 5 and top 10 for easier tracking)
    most_abberant_zscores = sorted([(0.0, chrm_arm, 'n/a') if zscore == 'n/a' else
                                    (abs(float(zscore)), chrm_arm, zscore)
                                    for chrm_arm, zscore in sample_z_score.items()], reverse=True)
    # output most abberant chromosome arms
    reduced_top5_out_file = f"{'.'.join(out_file.split('.')[:-2])}.zscores.top5aberrant.txt"
    with open(reduced_top5_out_file, "wt") as f_out:
        f_out.writelines([f'{chrarm}\t{zscore}\n' for _abs, chrarm, zscore in most_abberant_zscores[:5]] +
                         [f'genomewide\t{genomewide_z}'])
    reduced_top10_out_file = f"{'.'.join(out_file.split('.')[:-2])}.zscores.top10aberrant.txt"
    with open(reduced_top10_out_file, "wt") as f_out:
        f_out.writelines([f'{chrarm}\t{zscore}\n' for _abs, chrarm, zscore in most_abberant_zscores[:10]] +
                         [f'genomewide\t{genomewide_z}'])
