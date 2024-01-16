#!/usr/bin/env python3

from sys import argv as sys_argv
from typing import Tuple, TextIO
from pysam import AlignmentFile
from itertools import islice as itr_islice
from csv import reader as csv_reader
from pathlib import Path
from typing import Union as OneOf


def read_genomic_interval(interval_size: int, sam_handle: OneOf[AlignmentFile, TextIO]) -> Tuple[Tuple, bool]:
    sam_content = list(itr_islice(sam_handle, interval_size))
    empty_line_encountered = len(sam_content) != interval_size
    sam_content = tuple(filter(
        lambda x: not x.is_unmapped and x.mapping_quality > 14 and not x.is_qcfail and not x.is_secondary, sam_content))
    return sam_content, empty_line_encountered


if __name__ == '__main__':
    aln_file = Path(sys_argv[1])
    results_file = sys_argv[2]
    chr_arms_file = sys_argv[3]
    # read in bad chromosome arms
    script_parent_path = Path(__file__).parent
    with open(script_parent_path.parent / 'ref/exclude_chromosome_arms.csv', 'rt') as f_excl:
        exclude_low_mappability_arms = list(csv_reader(f_excl))[0]  # all in the first line
    # read in chromosome arm coordinates -------------------------------------------------------------------------------
    chromosome_arms = {}
    with open(chr_arms_file, "rt") as f_in:
        # looks something like the following:
        # chr1	0	125000000	chr1p
        # chr1	125000000	249250621	chr1q
        chr_arm_content = [line.strip().split('\t') for line in f_in.readlines()]
        try:
            for chrm, arm_start, arm_stop, arm_label in chr_arm_content:
                try:
                    chromosome_arms[chrm].update({arm_label: (int(arm_start), int(arm_stop))})
                except KeyError:
                    chromosome_arms.update({chrm: {arm_label: (int(arm_start), int(arm_stop))}
                                            })  # create entry for chromosome
        except ValueError:  # not enough values to unpack (expected 4, got 1) -> empty line
            print(f" !i : got arms for {len(chromosome_arms)} chromosomes but experienced an error.\n" +
                  "       Continuing with current chromosome arms ..")
            pass
    # process alignments -----------------------------------------------------------------------------------------------
    chromosome_arm_counts = {}
    all_arms = [f'chr{i}{arm}' for i in range(1, 23) for arm in ('p', 'q')] + \
               [f'chr{chrm}{arm}' for chrm in ('X', 'Y') for arm in ('p', 'q')]
    _1 = [None for idx, chromosome_arm_counts[all_arms[idx]] in enumerate([0] * len(all_arms))]  # populate dict with 0s
    strange_reads = 0
    if not any([aln_file.name.lower().endswith(valid_ending) for valid_ending in ('sam, bam')]):
        print(f"WARNING - the input alignment file '{aln_file.name}' does not end in an expected alignment file string")
    with AlignmentFile(aln_file, "r") as f_aln:  # figures out if is BAM or is SAM itself
        sam_fully_processed = False
        while not sam_fully_processed:
            chunk_content, incomplete_chunk_flag = read_genomic_interval(interval_size=100000, sam_handle=f_aln)
            sam_fully_processed = incomplete_chunk_flag
            for aligned_segment in chunk_content:  # TEST THIS
                chrom = aligned_segment.reference_name
                mp_start = aligned_segment.reference_start
                try:
                    short_start, short_stop = chromosome_arms[chrom][f'{chrom}p']
                except KeyError:  # e.g. unplaced contig 'chrUn_gl000219'
                    continue
                long_start, long_stop = chromosome_arms[chrom][f'{chrom}q']
                if long_start <= int(mp_start) < long_stop:  # mapping start inside short arm
                    if f'{chrom}q' in exclude_low_mappability_arms:
                        continue
                    chromosome_arm_counts[f'{chrom}q'] += 1
                elif short_start <= int(mp_start) < short_stop:
                    if f'{chrom}p' in exclude_low_mappability_arms:
                        continue
                    chromosome_arm_counts[f'{chrom}p'] += 1
                else:
                    strange_reads += 1
    if strange_reads != 0:
        print(f" !i : {strange_reads:,} strange R1 reads encountered.\n" +
              "       Continuing ..")
    # create output ----------------------------------------------------------------------------------------------------
    out_lines = [f"{chrm_arm}\t{chromosome_arm_counts[chrm_arm]}\n" for chrm_arm in all_arms]
    with open(results_file, "wt") as f_out:
        f_out.writelines(out_lines)
