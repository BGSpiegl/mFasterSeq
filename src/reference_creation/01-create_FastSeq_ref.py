#!/usr/bin/python3

# You need to create *.zscoresANDreadcounts.txt files first using the fastseq pipeline script first.
# WARNING: requires a 'exclude_chromosome_arms.csv' file specifying bad arms to exclude (acrocentric ones like chr14p)
#  i :  run from within python3.9 conda env

from glob import glob as glb_glob
from os.path import join as pth_join, sep as pth_sep, basename as pth_basename, \
    isdir as pth_isdir, dirname as pth_dirname, abspath as pth_abspath, isfile as pth_isfile
from os import walk as os_walk, makedirs as os_mkdirs
from re import match as re_match, DOTALL as re_DOTALL
from natsort import humansorted
from numpy import mean as np_mean, std as np_std
from argparse import ArgumentParser
import pandas as pd
import matplotlib.pyplot as plt
from csv import reader as csv_reader
from typing import List
from scipy.stats import shapiro, normaltest, kstest
from time import localtime, strftime
from shutil import copyfile as sh_cpfile


def get_cmdline_args():
    prsr = ArgumentParser()
    prsr.add_argument('-f', '--females-dir', dest='fem_dir', required=True,
                      help="Directory containing female samples to regard in female control creation.")
    prsr.add_argument('-m', '--males-dir', dest='mal_dir', required=True,
                      help="Directory containing male samples to regard in male control creation.")
    prsr.add_argument('-o', '--out-dir', dest='out_dir',
                      help="Directory in which created controls files shall be stored. If not " +
                           "specified, controls will be placed in the parent folder of the females dir.")
    prsr.add_argument('-p', '--prefix', dest='pref', help="Prefix added to the filename of the created reference file.")
    prsr.add_argument('-s', '--search-subdirectories', dest='search_subdirs',
                      help="Flag: if set, subdirectories of '--females-dir' and '--males-dir' will " +
                           "also be searched for *.chrarm.txt files. [DEFAULT: False]", action='store_true')
    return prsr.parse_args()


def get_folders(directory: str) -> list:
    return [dr[0] for dr in os_walk(directory)]  # returns full paths!


def get_custom_files(directory: str, name_pattern: str, search_subdirectories=True):
    # supports auto annotations
    cust_files = list()
    output_data = list()
    if '(?P<NAME>' not in name_pattern:
        raise AttributeError("name_patterns must include '(?P<NAME>'")
    if search_subdirectories:
        sub_dirs = get_folders(directory=directory)
        for sub_d in sub_dirs:
            search_result = glb_glob(pth_join(sub_d, "*." + name_pattern.split('.')[-1]))
            cust_files.extend(humansorted(search_result))
    else:
        search_result = glb_glob(pth_join(directory, "*." + name_pattern.split('.')[-1]))
        cust_files.extend(humansorted(search_result))
    for idx, cust_file in enumerate(cust_files):
        try:
            name_dict = re_match(name_pattern, pth_basename(cust_file), re_DOTALL).groupdict()
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
        raise IOError("provided name pattern yielded 0 files.")
    return output_data


def read_chrarm_data(data: list) -> dict:
    read_count_data = {'total_count': dict()}
    for read_count_dat in data:
        total_count = 0
        sn = read_count_dat['name']
        read_count_file_path = read_count_dat['file']
        file_total_read_count = None
        with open(read_count_file_path, 'rt') as f_read_counts:
            for line in f_read_counts:
                chrm_arm, read_count, *other = line.strip().split('\t')
                if chrm_arm == 'genomewide':
                    try:
                        read_count_data['sum_of_squares'].update({sn: float(other[0])})
                    except KeyError:
                        read_count_data.update({'sum_of_squares': {sn: float(other[0])}})
                    file_total_read_count = int(read_count)
                else:
                    try:
                        total_count += int(read_count)
                    except ValueError:  # ignore this chromosome arm because was likely 'n/a'
                        try:
                            read_count_data[chrm_arm].update({sn: 0})
                        except KeyError:
                            read_count_data.update({chrm_arm: {sn: 0}})
                        continue
                    try:
                        read_count_data[chrm_arm].update({sn: int(read_count)})
                    except KeyError:
                        read_count_data.update({chrm_arm: {sn: int(read_count)}})
        read_count_data['total_count'].update({sn: total_count})
        if file_total_read_count is None:
            print(f" !  WARNING: 'genomewide' read count line not encountered for sample '{sn}'")
        elif total_count != file_total_read_count:
            print(" !  WARNING: total read count from genomewide line in file and counted entries " +
                  f"for chromosme arms was off by {abs(total_count - int(read_count)):,} for sample '{sn}'")
    return read_count_data


def create_ref_read_count_data_and_plot(read_count_data: dict, out_file: str, exclude_arms: List[str], sex: str):
    sample_names = humansorted(read_count_data['total_count'].keys())
    chr_arm_names = list(read_count_data.keys())
    chr_arm_names.remove('total_count')
    chr_arm_names.remove('sum_of_squares')
    # create visualisation data structures
    raw_total_counts = {'raw total read counts': [read_count_data['total_count'][sn]
                                                  for sn in sample_names]}
    aligned_fraction_arm_per_sample_dict = {}
    # all_zscores = {}
    all_zscores_leave1out = {}
    # prepare normal distribution stats
    stats_out_read_fractions = pth_join(pth_dirname(out_file),
                                        f'{sex}_mFastSeq_readcountfractions_chromwise_tests_normal_distribution.tsv')
    stats_out_z_scores = pth_join(pth_dirname(out_file),
                                  f'{sex}mFastSeq_zscores_chromwise_tests_normal_distribution.tsv')
    # doc from https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.normaltest.html:
    # test_statistic: s^2 + k^2, where s is the z-score returned by skewtest and k is the z-score
    #                 returned by kurtosistest.
    # p-value: A 2-sided chi squared probability for the hypothesis test.
    #          (null hypothesis: x comes from a normal distribution)
    stats_alpha = 0.05
    stats_out_buffer_readfracs = [f"mFastSeq: test results for normal distribution of chromosome arm-wise read count " +
                                  f"fractions\n" +
                                  (("WARNING: not enough samples for kurtosistest for normaltest to be valid! (20 " +
                                   f"required but was {len(sample_names):,}\n") if len(sample_names) < 20 else '') +
                                  f"(alpha was {stats_alpha}; created on {strftime('%d/%m/%Y, %H:%M', localtime())})\n",
                                  "arm\tmean_read_count_fraction\tstddev_read_count_fraction\ttest_statistic\t" +
                                  "p_value\tShapiro-Wilk_statistic\tShapiro-Wilk_p_value\t" +
                                  "Kolmogorov-Smirnov_statistic\tKolmogorov-Smirnov_p_value\t" +
                                  "comes_from_normal_distribution\n"]
    stats_out_buffer_zscores = [f"mFastSeq: test results for normal distribution of chromosome arm-wise z-scores\n" +
                                (("WARNING: not enough samples for kurtosistest for normaltest to be valid! (20 " +
                                 f"required but was {len(sample_names):,})") if len(sample_names) < 20 else '') +
                                f"(alpha was {stats_alpha}; created on {strftime('%d/%m/%Y, %H:%M', localtime())})\n",
                                "arm\tmean_z-score\tstddev_z-score\tnormal_test_statistic\tp_value\t" +
                                "Shapiro-Wilk_statistic\tShapiro-Wilk_p_value\t" +
                                "Kolmogorov-Smirnov_statistic\tKolmogorov-Smirnov_p_value\t" +
                                "comes_from_normal_distribution\n"]
    # create output and do test for normal distribution
    read_count_line_buffer = list()
    for charm in chr_arm_names:  # process per arm
        if charm in exclude_arms:  # handle bad arms
            read_count_line_buffer.append(f'{charm}\t0.0\t0.0\n')  # use dummy line for bad chromosome arms
            stats_out_buffer_readfracs.append(f"{charm}\tn/a\tn/a\tn/a\tn/a\tn/a\tn/a\tn/a\tn/a\tn/a\n")
            stats_out_buffer_zscores.append(f"{charm}\tn/a\tn/a\tn/a\tn/a\tn/a\tn/a\tn/a\tn/a\tn/a\n")
            # all_zscores.update({charm: [0.0] * len(sample_names)})
            all_zscores_leave1out.update({charm: [0.0] * len(sample_names)})
            continue
        # compute fraction of alignments of current chrom. arm relative to total sample aignments (bad arms excluded)
        aligned_fraction_arm = [read_count_data[charm][sn] / read_count_data['total_count'][sn]
                                for sn in sample_names]  # for all samples
        aligned_fraction_arm_per_sample_dict.update({charm: aligned_fraction_arm})  # for distribution plots l8r
        mean_arm_fraction = np_mean(aligned_fraction_arm)  # of current chromosome arm
        # leave-one-out version:
        mean_arm_fraction_leave1out = [np_mean(aligned_fraction_arm[:i] + aligned_fraction_arm[i+1:])
                                       for i, _ in enumerate(aligned_fraction_arm)]
        stddev_arm_fraction = np_std(aligned_fraction_arm, ddof=1)  # standard deviation of current chromosome arm
        # degrees of freedom -1 for more accurate results for smaller sample sizes; equates results from EXCEL STDEV
        stddev_arm_fraction_leave1out = [np_std(aligned_fraction_arm[:i] + aligned_fraction_arm[i+1:], ddof=1)
                                         for i, _ in enumerate(aligned_fraction_arm)]
        # compute normal distribution stats for the fraction of reads aligned to current chromosome arm
        large_enough_for_normaltest = len(aligned_fraction_arm) > 19  # kurtosis test used only valid above 19 samples
        if large_enough_for_normaltest:
            test_stat, test_pval = normaltest(aligned_fraction_arm, nan_policy='omit')
        else:
            test_stat = test_pval = 'n/a'
        sh_stat, sh_pval = shapiro(aligned_fraction_arm)  # H0: "sample comes from normal distribution" like for normaltest
        ks_stat, ks_pval = kstest(aligned_fraction_arm, 'norm', alternative='two-sided')
        # reject H0 (data come from same (i.e. normal) distributions if p-value < alpha
        normal_distributed = (test_pval > stats_alpha or sh_pval > stats_alpha or ks_pval > stats_alpha) \
            if large_enough_for_normaltest else (sh_pval > stats_alpha or ks_pval > stats_alpha)
        stats_out_buffer_readfracs.append(f"{charm}\t{mean_arm_fraction}\t{stddev_arm_fraction}\t{test_stat}\t" +
                                          f"{test_pval}\t{sh_stat}\t{sh_pval}\t{ks_stat}\t{ks_pval}\t" +
                                          f"{normal_distributed}\n")
        if stddev_arm_fraction == 0.0:  # handle zero division case
            mean_arm_zscore = 'n/a'
            # all_zscores.update({charm: [0.0] * len(sample_names)})
            all_zscores_leave1out.update({charm: [0.0] * len(sample_names)})
            stddev_arm_zscore = 'n/a'
        else:
            arm_zscores = [(afa - mean_arm_fraction) / stddev_arm_fraction for afa in aligned_fraction_arm]
            arm_zscores_leave1out = [(afa - mean_arm_fraction_leave1out[af_idx]) / stddev_arm_fraction_leave1out[af_idx]
                                     for af_idx, afa in enumerate(aligned_fraction_arm)]
            all_zscores_leave1out.update({charm: arm_zscores_leave1out})
            # compute normal distribution stats for read fractions
            large_enough_for_normaltest = len(arm_zscores_leave1out) > 19  # kurtosis test used only valid above 19 samples
            if large_enough_for_normaltest:
                test_stat, test_pval = normaltest(arm_zscores, nan_policy='omit')
            else:
                test_stat = test_pval = 'n/a'
            sh_stat, sh_pval = shapiro(arm_zscores)
            ks_stat, ks_pval = kstest(arm_zscores, 'norm', alternative='two-sided')
            mean_arm_zscore = np_mean(arm_zscores)
            stddev_arm_zscore = np_std(arm_zscores, ddof=1)
            # why is this 1? -> because we computed a reference and we compared the distribution against itself.
        normal_distributed = (test_pval > stats_alpha or sh_pval > stats_alpha or ks_pval > stats_alpha) \
            if large_enough_for_normaltest else (sh_pval > stats_alpha or ks_pval > stats_alpha)
        stats_out_buffer_zscores.append(f"{charm}\t{mean_arm_zscore}\t{stddev_arm_zscore}\t{test_stat}\t{test_pval}\t" +
                                        f"{sh_stat}\t{sh_pval}\t{ks_stat}\t{ks_pval}\t{normal_distributed}\n")
        # BELOW is a line (one for each chromosome arm) in reference file
        read_count_line_buffer.append(f'{charm}\t{mean_arm_fraction}\t{stddev_arm_fraction}\n')
    # FINALIZE genomewide stats
    # add the total number of reads line
    all_raw_read_counts = raw_total_counts['raw total read counts']
    large_enough_for_normaltest = len(all_raw_read_counts) > 19  # kurtosis test used only valid above 19 samples
    if large_enough_for_normaltest:
        test_stat, test_pval = normaltest(all_raw_read_counts)
    else:
        test_stat = test_pval = 'n/a'
    sh_stat, sh_pval = shapiro(all_raw_read_counts)
    ks_stat, ks_pval = kstest(all_raw_read_counts, 'norm', alternative='two-sided')
    normal_distributed = (test_pval > stats_alpha or sh_pval > stats_alpha or ks_pval > stats_alpha) \
        if large_enough_for_normaltest else (sh_pval > stats_alpha or ks_pval > stats_alpha)
    stats_out_buffer_readfracs.append(f"raw_total_read_counts\t{np_mean(all_raw_read_counts)}\t" +
                                      f"{np_std(all_raw_read_counts, ddof=1)}\t{test_stat}\t{test_pval}\t{sh_stat}\t" +
                                      f"{sh_pval}\t{ks_stat}\t{ks_pval}\t{normal_distributed}\n")
    # compute sum of sample_arm_zscore ** 2 for each control sample: (mean and stddev required)
    # the next is computed on statistics created excluding the respective sample
    # -> MUST DO! otherwise the genomewide sum of squares will always be the number of included chromosome arms!!!
    #    !It is not correct to include a sample into the reference that is used in its own evaluation of how off it is!
    samplewise_sum_of_squared_zscores_leave1out = [sum([all_zscores_leave1out[chrarm][sample_id] ** 2
                                                        for chrarm in chr_arm_names])
                                                   for sample_id, _sample_name in enumerate(sample_names)]
    all_zscores_leave1out.update({'genomewide sum of zscore squares': samplewise_sum_of_squared_zscores_leave1out})
    mean_sum_of_squared_zscores_leave1out = np_mean(samplewise_sum_of_squared_zscores_leave1out)
    stddev_sum_of_squared_zscores_leave1out = np_std(samplewise_sum_of_squared_zscores_leave1out, ddof=1)
    # ad normaltest below: kurtosis test used only valid above 19 samples
    large_enough_for_normaltest = len(samplewise_sum_of_squared_zscores_leave1out) > 19
    if large_enough_for_normaltest:
        test_stat, test_pval = normaltest(samplewise_sum_of_squared_zscores_leave1out, nan_policy='omit')
    else:
        test_stat = test_pval = 'n/a'
    sh_test_l1o, sh_pval_l1o = shapiro(samplewise_sum_of_squared_zscores_leave1out)
    ks_stat_l1o, ks_pval_l1o = kstest(samplewise_sum_of_squared_zscores_leave1out, 'norm', alternative='two-sided')
    normal_distributed = (test_pval > stats_alpha or sh_pval > stats_alpha or ks_pval > stats_alpha) \
        if large_enough_for_normaltest else (sh_pval > stats_alpha or ks_pval > stats_alpha)
    stats_out_buffer_zscores.append(f"genomewide_sum_of_zscore_squares\t{mean_sum_of_squared_zscores_leave1out}\t" +
                                    f"{stddev_sum_of_squared_zscores_leave1out}\t{test_stat}\t{test_pval}\t" +
                                    f"{sh_test_l1o}\t{sh_pval_l1o}\t{ks_stat_l1o}\t{ks_pval_l1o}\t" +
                                    f"{normal_distributed}\n")
    # BELOW is last line in reference file
    read_count_line_buffer.append(f"genomewide\t{mean_sum_of_squared_zscores_leave1out}\t" +
                                  f"{stddev_sum_of_squared_zscores_leave1out}")
    # write outputs
    with open(out_file, 'wt') as f_read_counts_out:
        f_read_counts_out.writelines(read_count_line_buffer)
    with open(stats_out_read_fractions, 'wt') as f_readfrac_stats:
        f_readfrac_stats.writelines(stats_out_buffer_readfracs)
    with open(stats_out_z_scores, 'wt') as f_zscores_stats:
        f_zscores_stats.writelines(stats_out_buffer_zscores)
    # create visualisations of read count distributions
    out_dir = pth_dirname(out_file)
    # NOW: raw total counts
    vis_data_rawread = pd.DataFrame(data=raw_total_counts)
    data_label = 'raw total read counts'
    cur_fname = pth_join(out_dir, f"{sex}_{data_label.replace(' ', '_')}.totalreads_dist.png")
    _n, _bins, _patches = plt.hist(vis_data_rawread[data_label], bins=21)
    plt.ylabel('frequency / 1')
    plt.xlabel('total read count / 1')
    plt.title(f"Distribution of Total Read Count ({sex})")
    plt.grid(True)
    plt.savefig(cur_fname, dpi=150)
    plt.close()
    # NOW: normalized counts
    # pd.DataFrame WARNING: don't mix datatypes!
    vis_data_readfracs = pd.DataFrame(data=aligned_fraction_arm_per_sample_dict)
    for data_label in chr_arm_names:
        if data_label in exclude_arms:
            continue
        try:
            cur_fname = pth_join(out_dir, f"{sex}_{data_label.replace(' ', '_')}.readfraction_dist.png")
            _n, _bins, _patches = plt.hist(vis_data_readfracs[data_label], bins=21)
            plt.ylabel('frequency / 1')
            plt.xlabel('normalized read count / 1')
            plt.title(f"Distribution of {sex} Normalized Read Count for {data_label}")
            plt.grid(True)
            plt.savefig(cur_fname, dpi=150)
            plt.close()
        except KeyError:
            continue
    # NOW: Z-scores
    ordered_vis_labels = chr_arm_names + ['genomewide sum of zscore squares']
    vis_data_zscores = pd.DataFrame(data=all_zscores_leave1out)
    for data_label in ordered_vis_labels:
        if data_label in exclude_arms:
            continue
        cur_fname = pth_join(out_dir, f"{sex}_{data_label.replace(' ', '_')}.zscores_dist.png")
        _n, _bins, _patches = plt.hist(vis_data_zscores[data_label], bins=11)
        plt.ylabel('frequency / 1')
        if data_label == 'genomewide sum of zscore squares':
            plt.xlabel('genomewide sum of squared Z-scores / 1')
            plt.title(f'Distribution of {sex} Sum of Squared Z-scores (genomewide)')
        else:
            plt.xlabel('Z-score / 1')
            plt.title(f"Distribution of {sex} Z-scores for {data_label}")
        plt.grid(True)
        plt.savefig(cur_fname, dpi=150)
        plt.close()


if __name__ == '__main__':
    cmdline_args = get_cmdline_args()
    female_controls_dir = cmdline_args.fem_dir
    if female_controls_dir not in (None, '') and not pth_isdir(female_controls_dir):
        raise AttributeError(" female controls directory is non-existent.")
    male_controls_dir = cmdline_args.mal_dir
    if male_controls_dir not in (None, '') and not pth_isdir(male_controls_dir):
        raise AttributeError(" female controls directory is non-existent.")
    if male_controls_dir in (None, '') and female_controls_dir in (None, ''):
        raise AttributeError(" both female and male control directories were undefined.")
    out_dir = cmdline_args.out_dir
    search_subdirs = cmdline_args.search_subdirs
    ref_prefix = cmdline_args.pref
    script_path = pth_dirname(pth_abspath(__file__))
    with open(pth_join(script_path, '../FastSeq/exclude_chromosome_arms.csv'), 'rt') as f_excl:
        EXCLUDE_LOW_MAPPABILITY_ARMS = list(csv_reader(f_excl))[0]  # all in first line
    if out_dir in (None, ''):
        dir_parts = list(filter(lambda x: x != '', female_controls_dir.split(pth_sep)))
        out_dir_fem = pth_sep.join(dir_parts[:-1])
        dir_parts = list(filter(lambda x: x != '', male_controls_dir.split(pth_sep)))
        out_dir_mal = pth_sep.join(dir_parts[:-1])
    else:
        out_dir_fem = pth_join(out_dir, 'female_control_results')
        out_dir_mal = pth_join(out_dir, 'male_control_results')
    out_dirs = {'female': out_dir_fem, 'male': out_dir_mal}
    created_refs = []
    for (read_counts_parent_directory, sex) in ((female_controls_dir, 'female'), (male_controls_dir, 'male')):
        if not pth_isdir(out_dirs[sex]):
            os_mkdirs(out_dirs[sex])
        print(f"Fetching read count files for {sex}...")
        read_count_files_data = get_custom_files(directory=read_counts_parent_directory,
                                                 name_pattern=r'(?P<NAME>.*).zscoresANDreadcounts.txt',
                                                 search_subdirectories=search_subdirs)
        print(f"-> {len(read_count_files_data):,} files found.")
        new_reference_file = pth_join(out_dirs[sex],
                                      f"{ref_prefix if ref_prefix else ''}FastSeq_controls_{sex}.txt")
        print("----------------------------------------------------")
        print("Reading read count data...")
        chrm_arm_data = read_chrarm_data(data=read_count_files_data)  # only includes not-nan counts
        print(f"-> data for {len(chrm_arm_data['total_count'])} samples created.")
        print("----------------------------------------------------")
        print("creating new reference file...")
        create_ref_read_count_data_and_plot(read_count_data=chrm_arm_data, out_file=new_reference_file,
                                            exclude_arms=EXCLUDE_LOW_MAPPABILITY_ARMS, sex=sex)
        if pth_isfile(new_reference_file):
            created_refs.append(new_reference_file)
        print(f"Successfully created new reference file for {sex}.")
    # copy created reference files to out dir (only if new directories were created in output directory)
    for created_ref in created_refs:
        if out_dir is not None:
            dest_path = pth_join(out_dir, pth_basename(created_ref))
            sh_cpfile(created_ref, dest_path)
