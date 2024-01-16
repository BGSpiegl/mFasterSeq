#!/usr/bin/python3

from matplotlib import pyplot as plt
from numpy import array as np_array, percentile as np_percentile, concatenate as np_concatenate, linspace as np_linspace
from scipy import stats
from typing import Dict, Tuple
from glob import glob as glb_glob
from os import makedirs as os_mkdirs
from os.path import join as pth_join, basename as pth_basename, isdir as pth_isdir
from sys import exit as sys_exit
from natsort import humansorted
from math import ceil as mth_ceil, floor as mth_floor


# path definitions - YOU MUST RUN 01-create_FastSeq_ref.py FIRST!
reference_files = {'previous': {'m': '/<SOME-PATH>/FastSeq/ref/' +
                                     'OLD_MaleControls.txt',
                                'f': '/<SOME-PATH>/FastSeq/ref/' +
                                     'OLD_FemaleControls.txt'},
                   'current': {'m': '/<SOME-PATH>/FastSeq/' +
                                    'NEW_FastSeq_controls_male.txt',
                               'f': '/<SOME-PATH>/FastSeq/' +
                                    'NEW_FastSeq_controls_female.txt'}}

# new, created files
read_fraction_normal_distribution_files = {'m': '/<SOME-PATH>/male_control_results/' +
                                                'male_mFastSeq_readcountfractions_chromwise_tests_normal_' +
                                                'distribution.tsv',
                                           'f': '/<SOME-PATH>/female_control_results/' +
                                                'female_mFastSeq_readcountfractions_chromwise_tests_normal_' +
                                                'distribution.tsv'}
sum_o_squared_zscores_normal_distribution_files = {'m': '/<SOME-PATH>/male_control_results/' +
                                                        'malemFastSeq_zscores_chromwise_tests_normal_distribution.tsv',
                                                   'f': '/<SOME-PATH>/female_control_' +
                                                        'results/femalemFastSeq_zscores_chromwise_tests_normal_' +
                                                        'distribution.tsv'}

analyses_outputs = {'previous': '/<SOME-PATH>/FastSeq_Controls/' +
                                'results_previous_controls/Sex_guessed',
                    'current': '/<SOME-PATH>/FastSeq_Controls/' +
                               'control_samples_with_new_controls/Sex_guessed'}
read_counts_file_pattern = '*.zscoresANDreadcounts.txt'  # WILDCARD MUST INCLUDE * or +!
zscores_file_pattern = '*.zscores.txt'  # WILDCARD MUST INCLUDE * or +!

output_directory = '/<SOME-PATH>/report'


# sexy, ugly functions (a combination of these traits in one entity is only possible in coding it seems ..)
def read_reference_data(data: str) -> Dict[str, Dict[str, Dict[str, Tuple[float, float]]]]:
    ref_data = {}.fromkeys(data)
    for analysis_group, group_based_data in data.items():
        ref_data[analysis_group] = {}.fromkeys(data[analysis_group])
        for bio_logical_sex, ref_path in data[analysis_group].items():
            with open(ref_path, 'rt') as f_ref:
                # file looks like this (arm, arm_aln_fraction, std dev arm_aln_fraction):
                # chr1p	0.03195813313974712	0.000556366721912164
                # chr1q	0.03627571135611149	0.0007578353342833798
                # ..
                # and last line:
                # region,
                # mean(sum_over_all(squared(arm Z-scores of individual sample))),
                # std_dev(sum_over_all(squared(arm Z-scores of individual sample)))
                # genomewide	41.00000000000001	20.657322022516226
                # below yields Tuple[str, Tuple[float, float]] for dictionary creation
                arm_reference_data = dict(tuple(map(lambda x: (x[0], (float(x[1]), float(x[2]))),
                                                    filter(lambda x: x[1] != '0.0',
                                                           map(lambda x: x[:3],
                                                               [str_line.strip().split('\t')
                                                                for str_line in f_ref.readlines()]
                                                               )))))  # 'n/a' arms removed
            ref_data[analysis_group][bio_logical_sex] = arm_reference_data
    return ref_data


def get_aligned_reads(in_dirs: Dict[str, str], file_pattern: str) \
        -> Tuple[Dict[str, Dict[str, Dict[str, np_array]]]]:
    # As a programmer, who had to read a lot of other peoples code, I want to apologize in advance for this code.
    # Here, I tried to process data as fast and in as few lines as possible without spending hours on the matter.
    if '*' not in file_pattern and '+' not in file_pattern:
        print("ERROR: file pattern did not include '*' or '+'. Try again with different name pattern.")
        sys_exit(1)
    if '*' in file_pattern:
        split_char = '*'
    else:
        split_char = '+'
    # local function defs
    name_extractor_function = lambda s: pth_basename(s).split(file_pattern.split(split_char)[1]
                                                              if file_pattern[0] == split_char else
                                                              file_pattern.split(split_char)[0])[0]

    aligned_counts = {}.fromkeys(in_dirs)
    aligned_fractions = {}.fromkeys(in_dirs)
    # gender_splitter_function = lambda t: (aligned_counts[analysis_group]['f'].update({t[0], t[1]})
    #                                       if sample_name[0].lower() == 'f' else
    #                                       aligned_counts[analysis_group]['m'].update({t[0], t[1]}))
    # structure: {analysis_group: {biological_sex: {chromosome_arm: np_array(list_of_values)}}}}
    # 2 iterations - split biological sex based on first char!
    exclude_strings = ('n/a', 'na', 'none')
    for analysis_group, group_based_dir_path in in_dirs.items():
        aligned_counts[analysis_group] = {'f': None, 'm': None}
        aligned_fractions[analysis_group] = {'f': None, 'm': None}
        in_dir = pth_join(group_based_dir_path, file_pattern)
        all_count_files = glb_glob(in_dir)
        if not all_count_files:
            if pth_isdir(in_dir):
                print(f"ERROR: no files matching pattern '{file_pattern}' found in directory '{in_dir}'.\n" +
                       "       Check your pattern!", flush=True)
            else:
                print(f"ERROR: no files matching pattern '{file_pattern}' found in directory '{in_dir}'.\n" +
                      f"       Directory '{in_dir}' does not exist!", flush=True)
            sys_exit(1)
        try:
            all_sample_names = [name_extractor_function(cur_pth) for cur_pth in all_count_files]
        except (IndexError,):
            print("ERROR: there was a problem extracting the sample name from your files. Check your name pattern!")
            sys_exit(1)
        if '' in all_sample_names:
            print("ERROR: empty sample name encountered (''). Check your name pattern!")
            sys_exit(1)
        for sample_name, aln_path in zip(all_sample_names, all_count_files):
            be_logical_sex = sample_name[0].lower()
            if be_logical_sex not in ('f', 'm'):
                print(f"ERROR: sample '{sample_name}' could not be identified as either female or male based " +
                      "on first letter.")
                sys_exit(1)
            with open(aln_path, 'rt') as f_counts:
                # file looks like this (arm, count, arm_aln_fraction, z-score):
                # chr1p	13851	0.032577863704999930	1.3169136618369781
                # chr1q	15498	0.036451644769337155	0.3376005000841651
                # ..
                input_lines = f_counts.readlines()
                sample_arm_counts = dict(tuple(map(lambda z: (z[0], int(z[1])),
                                                   filter(lambda y: y[1].lower() not in exclude_strings,
                                                          map(lambda x: x[:2],
                                                              [str_line.strip().split('\t')
                                                               for str_line in input_lines]
                                                              )))))  # 'n/a' arms removed
                sample_arm_fractions = dict(tuple(map(lambda z: (z[0], float(z[2])),
                                                      filter(lambda y: y[1].lower() not in exclude_strings,
                                                             map(lambda x: x[:3],
                                                                 [str_line.strip().split('\t')
                                                                  for str_line in input_lines]
                                                                 )))))  # 'n/a' arms removed
            if aligned_counts[analysis_group][be_logical_sex] is None:  # initialize once
                aligned_counts[analysis_group][be_logical_sex] = {}.fromkeys(sample_arm_counts)
                _initialize_empty = [None for cur_arm, aligned_counts[analysis_group][be_logical_sex][cur_arm] in
                                     zip(aligned_counts[analysis_group][be_logical_sex].keys(),
                                         ([] for _i in range(len(sample_arm_counts))))]  # initialize with empty list

                aligned_fractions[analysis_group][be_logical_sex] = {}.fromkeys(sample_arm_fractions)
                _initialize_empty = [None for cur_arm, aligned_fractions[analysis_group][be_logical_sex][cur_arm] in
                                     zip(aligned_fractions[analysis_group][be_logical_sex].keys(),
                                         ([] for _i in range(len(sample_arm_fractions))))]  # initialize with empty list
            _fill_in_counts = [aligned_counts[analysis_group][be_logical_sex][cur_arm].append(cur_count)
                               for cur_arm, cur_count in sample_arm_counts.items()]
            _fill_in_fractions = [aligned_fractions[analysis_group][be_logical_sex][cur_arm].append(cur_frac)
                                  for cur_arm, cur_frac in sample_arm_fractions.items()]
        # cast content of lists to float numpy arrays (after processing all samples of differnet biological sexes):
        for b_s in ('f', 'm'):
            for arm in aligned_counts[analysis_group][b_s]:
                aligned_counts[analysis_group][b_s][arm] = np_array(aligned_counts[analysis_group][b_s][arm])
                aligned_fractions[analysis_group][b_s][arm] = np_array(aligned_fractions[analysis_group][b_s][arm])
    return aligned_counts, aligned_fractions


def create_freedman_diacony_rule_bins(data=None, n_times=1.):
    if data is None:
        raise AttributeError("you must provide an iterable object to calculate number of bins.")
    data_max = max(data)
    data_min = min(data)
    interquartile_range = np_percentile(data, q=75, interpolation='midpoint') - \
                          np_percentile(data, q=25, interpolation='midpoint')
    n = len(data)
    binwidth = 2 * interquartile_range / (n ** (1 / 3))
    is_integer_data = isinstance(data[0], int)
    if not is_integer_data and interquartile_range < 2:  # compute float bins
        ten_power = 1
        while binwidth * 10**ten_power < 1:
            ten_power += 1
            nicer_binwidth = round((binwidth * 10**ten_power)) / 10**ten_power
        binwidth = nicer_binwidth / n_times
        # alternate_x_start = mth_floor(data_min*10**ten_power)/10**ten_power
        # alternate_x_stop = mth_ceil(data_max*10**ten_power)/10**ten_power
        data_min = mth_floor(data_min / binwidth) * binwidth
        data_max = mth_ceil(data_max / binwidth) * binwidth
        n_bins = mth_ceil((data_max - data_min) / binwidth)
        bins_to_return = [data_min + i * binwidth for i in range(n_bins + 1)]
    else:  # compute integer bins
        if binwidth == 0:
            binwidth = 1
        else:
            ten_power = 1
            while binwidth / 10**ten_power > 10:
                ten_power += 1
            nicer_binwidth = round((binwidth / 10**ten_power)) * 10**ten_power
            data_min = mth_floor(data_min / nicer_binwidth) * nicer_binwidth
            data_max = mth_ceil(data_max / nicer_binwidth) * nicer_binwidth
            binwidth = nicer_binwidth
        binwidth = int(binwidth / n_times)
        bins_to_return = list(range(data_min, data_max + binwidth, binwidth))
    return bins_to_return


def load_norm_test_results(paths: Dict[str, str], alpha=0.05) -> Dict[str, Tuple[str, bool]]:
    # skip to line starting with 'arm'; get values for column 'Shapiro-Wilk_p_value'
    # -> if value is > alpha=0.05, data for chrom. arm was drawn from a normally distributed random variable
    return_data = {'m': None, 'f': None}
    for sex, f_path in paths.items():
        with open(f_path, 'rt') as f_in:
            for line in f_in:
                line_content = line.strip().split('\t')
                if line_content[0] != 'arm':
                    continue
                shapiro_wilk_pvalue_column = line_content.index('Shapiro-Wilk_p_value')
                return_data[sex] = dict(map(lambda l: (l[0], float(l[shapiro_wilk_pvalue_column]) > alpha),
                                            filter(lambda l: l[shapiro_wilk_pvalue_column] != 'n/a',
                                                   [line.strip().split('\t') for line in f_in.readlines()]))
                                        )
    return return_data


if __name__ == '__main__':
    # create a bar chart as follows:
    # chart = alt.Chart(data).mark_bar().encode(
    #                                           x='a',
    #                                           y='average(b)'
    #                                           )
    # saving as follows:
    # chart.save('chart.json')

    # plan: create figures with old and new read fraction distributions (per chrom arm for all control samples;
    #       afterwards for sum of squared residuals -> is genomewide actually valid?)
    #       -> theoretical dists and real + info about whether KS-test says data is sample coming from a
    #          normal distribution; make intermediate plots for presentation -> only old ; then only new
    # ! we don't have definitive raw read counts for old reference creation!

    # PREPARE:
    os_mkdirs(output_directory, exist_ok=True)
    # 1) load mFastSeq reference files
    reference_values = read_reference_data(data=reference_files)  # create theoretical distributions from these vals
    # 2) load "aligned read fraction originating from normal distribution" test results
    #    (True means data is sample drawn from a normal distribution)
    normal_disted_aln_fractions = load_norm_test_results(paths=read_fraction_normal_distribution_files, alpha=0.05)

    # 3) load "sum of squared z-scores originating from normal distribution" test results
    #    (True means data is sample drawn from a normal distribution)
    normal_disted_sumozscores = load_norm_test_results(paths=sum_o_squared_zscores_normal_distribution_files, alpha=0.05)

    # PART 1: investigate the number of aligned reads per chromosome arm (and biological sex)
    # read in all read counts
    aligned_reads_reads_per_arm, aligned_reads_fraction_per_arm = get_aligned_reads(
        in_dirs=analyses_outputs, file_pattern=read_counts_file_pattern)
    # aligned_*_per_arm is of structure: {analysis_group: {biological_sex: {chromosome_arm: [list_of_values]}}}}

    # create bar chart for each biological sex and each arm that shows increase in % over previous results
    # -> in ascending order
    num_aln_perc_diff_all_sex = {'f': None, 'm': None}
    for bio_logical_sex in ('f', 'm'):
        all_current_arms_ordered = humansorted(aligned_reads_reads_per_arm['current'][bio_logical_sex].keys())
        for arm in all_current_arms_ordered:
            for data_label, (previous_data, current_data) in (
                    ('Aligned Read Counts', (aligned_reads_reads_per_arm['previous'][bio_logical_sex][arm],
                                             aligned_reads_reads_per_arm['current'][bio_logical_sex][arm])),
                    ('Aligned Read Fractions', (aligned_reads_fraction_per_arm['previous'][bio_logical_sex][arm],
                                                aligned_reads_fraction_per_arm['current'][bio_logical_sex][arm]))
            ):
                if arm == 'genomewide':
                    if 'Aligned Read Counts':
                        data_label = 'Sum of Squared Z-Scores'
                    elif data_label == 'Aligned Read Fractions':
                        data_label = 'Standard Deviation of Sum of Squared Z-Scores'
                current_hist_path = pth_join(output_directory, f"distribution_{data_label.lower().replace(' ', '_')}_" +
                                             f"{bio_logical_sex}_{arm}.png")
                bins = create_freedman_diacony_rule_bins(data=np_concatenate((previous_data, current_data)), n_times=2.)
                fig, ax = plt.subplots(figsize=(16, 9))
                # plot histogram (previous number of aligned reads (and, thus, fractions) is teh same since aligment is
                # not depending on refernce! only Z-scores are!
                n_cur, cur_bins, cur_patches = ax.hist(current_data, bins=bins, label='current reference',
                                                       color='blue', rwidth=1., alpha=0.5)
                # plot theoretical distributions if fractions:
                if data_label == 'Aligned Read Fractions':
                    mu_prev, sigma_prev = reference_values['current'][bio_logical_sex][arm]
                    dist_space = np_linspace(min(bins), max(bins), 501)  # from where to where, how many datapoints
                    prev_dist = stats.norm.pdf(dist_space, mu_prev, sigma_prev)
                    prev_dist /= (prev_dist.max() / n_cur.max())
                    ax.plot(dist_space, prev_dist, '--', color='darkorange', label='theoretical (previous reference)')
                    mu_cur, sigma_cur = reference_values['current'][bio_logical_sex][arm]
                    cur_dist = stats.norm.pdf(dist_space, mu_cur, sigma_cur) * sum(n_cur)
                    cur_dist /= (cur_dist.max() / n_cur.max())
                    ax.plot(dist_space, cur_dist, '--', color='darkblue', label='theoretical (current reference)')
                fig.suptitle(f'Distribution of Chromosome Arm {data_label} for ' +
                             f"{'Male' if bio_logical_sex[0].lower() == 'm' else 'Female'} Samples on {arm}",
                             fontsize=20)
                try:
                    if data_label == 'Aligned Read Fractions':
                        cur_test_result = normal_disted_aln_fractions[bio_logical_sex][arm]
                        plt.title(f"""(sample {"was" if cur_test_result else "wasn't"} drawn from a """ +
                                  """normal distribution)""")
                    elif data_label == 'Sum of Squared Z-Scores':
                        cur_test_result = normal_disted_sumozscores
                        plt.title(f"""(sample {"was" if cur_test_result else "wasn't"} drawn from a """+
                                  """normal distribution)""")
                except:
                    plt.title(f'(no test data for arm)')
                # don't add title for read counts
                ax.grid(axis='y', alpha=0.75)
                ax.legend(loc='upper right', fontsize=14)
                plt.ylabel('counts', fontsize=16)
                plt.xlabel(f'{data_label.lower()}', fontsize=16)
                plt.savefig(current_hist_path, dpi=150, format='png')
                plt.close(fig=fig)
