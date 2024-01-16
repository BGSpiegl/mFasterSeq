
# FASTerSeq Analysis

## v2.1.2

-------------------
## TODO: FINISH THIS README.md!!!

## Content
* Purpose
* Contributors
* Installation
* Usage
* Input
* Output
* Required Files

-------------------

chrom_arm_hg19.txt: provides coordinate definitions of chromosome arms

-------------------

fast_seq_analysis_chr_arm.py: Python script to count reads per chromosome arm 

## Purpose
Computing chromosome-arm-wide copy numbers from shallow sequenced WGS cfDNA liquid biopsy data 
(single-end or paired-end).
Currently, the pipeline **as-is** supports only hg19.
Code for creating a new reference is available. The reference (biological sex-matched controls) has to match the genome 
build used by mFastSeq! Comparing results between samples makes only sense if they were computed with THE SAME reference.
If you change the reference, you must reanalyze your entire cohort!

## Contributors
* original pipeline: Peter Ulz
* upgraded pipeline: Benjamin Spiegl

## Installation

Please use the conda environment on a UNIX machine. Tested only on Ubuntu (Debian).

## USAGE

### Synopsis:

`run_mFastSeq.py [-h] [--version] -i IN_DIR [-r REF_GENOME_FASTA] [-b] [-x SEX_FILE] [-afc ALTERNATIVE_FEMALE_CONTROL]
[-amc ALTERNATIVE_MALE_CONTROL] [-p N_PARALLEL] [-o OUT_DIR] [-vp] [-ape ADD_PATHS] [-afnp ALTERNATIVE_R1_FASTQ_PATTERN]
[-sd] [-k]`


### FULL OPTIONS:
`
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -i File, --input-directory File
                        Directory containing your sample FastQs (sequencing data). Use absolute paths!
  -r File, --reference-genome File
                        Path to your local indexed reference genome file. If the index files are not found, the pipeline will attempt to create them using BWA, samtools and picard. This argument is required if you run the mFastSeq pipeline
                        from the beginning (including read alignment)
  -b, --bam-input       Flag: if set, the input directory will be searched for BAM files instead of FastQ files. Use this to skip the alignment step if your samples have already been aligned. Only searches for BAM files with name matching
                        '(?P<NAME>).*'.sorted.bam. If another name pattern should be used, please specify via the '--alternative-fastq-name-pattern' parameter. Should also recognize SAM files.
  -x File, --sexs File  Path to comma-separated text file (*.csv) containing the biological sexs of your samples. Sex can be left unspecified. One sample per line: '<sample ID>,[optional sex specification: either 'f' or 'm']'. A valid
                        example for such a file would be: P1_S18,m<next line>P2_S7,<next line>P3_S9,f<end of file>If no sex code (either m or f) is given, the sex will be guessed based on the X-chromosomal portion of all reads. The name you
                        specify must follow the name pattern used in sample ID extraction. For instance, the default naming pattern would take 'W-P791' from the file 'W-P791_S2_R1.fastq.gz' which would be the Illumina sample name up to
                        '_S?'.
  -afc File, --alternative-female-control File
                        If you don't want to use the old initial controls on which the old plasmaDB and pipeline were based on, provide an alternative FEMALE control here. Per default, the reference files created on te 11th of January 2022
                        will be used (i.a. '2022-01-11_FastSeq_controls_female.txt')
  -amc File, --alternative-male-control File
                        If you don't want to use the old initial controls on which the old plasmaDB and pipeline were based on, provide an alternative MALE control here. Per default, the reference files created on te 11th of January 2022
                        will be used (i.a. '2022-01-11_FastSeq_controls_male.txt')
  -p N_PARALLEL, --parallel-processes N_PARALLEL
                        Number of processes that should be used in computation. The more the better. Do not use more than you have! [default: 4]
  -o File, --output-directory File
                        OPTIONAL: if provided, pipeline output will be placed in this directory. Otherwise, output will be placed in parent directory of input FASTQ directory under '<DATE>_mFastSeq_output'. The first <DATE> part follows the
                        international standard date and time notation: https://www.cl.cam.ac.uk/~mgk25/iso-time.html. Existing files will be overwritten!
  -vp, --vertical-plot  Flag: if this flag is set, also a vertical plot of the chromosome arm profiles will be created.
  -ape ADD_PATHS, --add-paths-to-environment ADD_PATHS
                        TROUBLESHOOTING: if you get an error message stating an executable for samtools, picard tools, or bwa was not found, you probably must add a path to one of your local directories containing those to the 'pythonpath'
                        so that the pipeline can find them. Provide these here (multiple paths possible by re-declaring parameter like this: '-ape <path_1> -ape <path_2>' etc.).
  -afnp ALTERNATIVE_R1_FASTQ_PATTERN, --alternative-fastq-name-pattern ALTERNATIVE_R1_FASTQ_PATTERN
                        TROUBLESHOOTING: if the pipeline tells you that no FastQs were found, check the fastq folder input you provided. If the provided path is correct, your FastQ files might be named differently than the default search
                        pattern expects: r'(?P<NAME>.*)_S.+_R1.*.fastq.gz' (default Illumina file naming). In this case, you must rename your files or provide an alternative RegEx pattern via this parameter. You might need to double-escape
                        backslashes. WARNING: MUST contain an Identifier for the 'R1' read if R1 and R2 read fastq files are present in input directory! Can be used to specify BAM file name pattern.
  -sd, --search-subdirectories
                        Flag: if set, the input directory's subdirectories will also be searched for fastQ files matching the (default/provided) name pattern.
  -k, --keep-sam        Do not delete the SAM file!`
