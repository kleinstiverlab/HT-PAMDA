# REQUIRED ARGUMENTS
RUN_NAME = 'example_PAMDA_data'
BARCODE_CSV = 'barcode_csv/example_PAMDA_barcodes.csv'
FASTQ_DIR = 'fastqs/example_PAMDA_data'
TIMEPOINT_FASTQ = {
                   'expRW086_pool_10_S10': 0,
                   'expRW086_pool_11_S11': 1, 
                   'expRW086_pool_12_S12': 2
                  }
PAM_ORIENTATION = 'three_prime' # or 'five_prime'
PAM_LENGTH = 4
PAM_START = 0
CONTROL_RAW_COUNT_CSV = 'output/example_PAM_library/PAMDA_1_raw_counts.csv.gz'
CONTROL_SAMPLE = 'control_sample'

# ADDITIONAL ARGUMENTS (with default values entered)
CONTROL_SAMPLE_TIMEPOINT_FASTQ = None
TIMEPOINTS = [0, 60, 480, 1920]
MAX_PAM_LENGTH = 8
SPACERS = {'SPACER1':'GGGCACGGGCAGCTTGCCGG', 
           'SPACER2':'GTCGCCCTCGAACTTCACCT'}
P5_SAMPLE_BARCODE_START = 2
P7_SAMPLE_BARCODE_START = 2
USE_TIMEPOINTS = None
TOP_N_NORMALIZE = 5
INIT_RATE_EST = [0.0001, 0.001, 0.01]
READ_SUM_MIN = 4
TPS_SUM_MIN = 1
PAM1_NT_RANK = {1:'A',2:'C',3:'G',4:'T'}
PAM2_NT_RANK = {1:'A',2:'C',3:'G',4:'T'}
PAM1_INDEX_RANK = None
PAM2_INDEX_RANK = None
AVERAGE_SPACER = True
HEATMAP_FIXED_MIN = -5.0
HEATMAP_FIXED_MAX = -1.5
LOG_SCALE_HEATMAP = True
CSV_INPUT_RAW_COUNT = None
CSV_INPUT_NORM_COUNT = None
CSV_INPUT_RATES = None

"""
Argument definitions:

`RUN_NAME`: Provide a name for the run (no spaces or special characters).

`BARCODE_CSV`: Filepath to a CSV with the barcode information and the following headers:

|sample	            |description|P5_sample_barcode|P7_sample_barcode|
|:------------------|:----------|:---------------:|:---------------:|
|unique_sample_id_01|WT_SpCas9  |CCTG             |TGAC             |
|unique_sample_id_02|SpG        |ATTC             |TACT             |
|unique_sample_id_03|SpRY       |CTAC             |ACCG             |

Barcodes should be provided 5' to 3' as they appear in the primer sequences that encode them. Avoid the use of special characters.

`FASTQ_DIR`: Folder directory contianing all fastq files.

`TIMEPOINT_FASTQ`: A Python dictionary to indicate the timepoint of each set of fastq files. The fastq determines the timepoint (except for the untreated library control). Indexing of timepoints must start at zero (first timepoint = 0, second timepoint = 1, etc.). Timepoint fastq names must be abbreviated as demonstrated in the following example: 

|fastq complete file names	        |key in TIMEPOINT_FASTQ dictionary|timepoint value|
|:----------------------------------|:-------------------------------:|:-------------:|
|timepoint1_S1_L001_R1_001.fastq.gz |timepoint1_S1                    |0              |
|timepoint1_S1_L001_R2_001.fastq.gz |                                 |               |
|timepoint1_S1_L002_R1_001.fastq.gz |                                 |               |
|timepoint1_S1_L002_R2_001.fastq.gz |                                 |               |
|timepoint1_S1_L003_R1_001.fastq.gz |                                 |               |
|timepoint1_S1_L003_R2_001.fastq.gz |                                 |               |
|timepoint1_S1_L004_R1_001.fastq.gz |                                 |               |
|timepoint1_S1_L004_R2_001.fastq.gz |                                 |               |
|                                   |                                 |               |
|timepoint2_S2_L001_R1_001.fastq.gz |timepoint2_S2                    |1              |
|timepoint2_S2_L001_R2_001.fastq.gz |                                 |               |
|timepoint2_S2_L002_R1_001.fastq.gz |                                 |               |
|timepoint2_S2_L002_R2_001.fastq.gz |                                 |               |
|timepoint2_S2_L003_R1_001.fastq.gz |                                 |               |
|timepoint2_S2_L003_R2_001.fastq.gz |                                 |               |
|timepoint2_S2_L004_R1_001.fastq.gz |                                 |               |
|timepoint2_S2_L004_R2_001.fastq.gz |                                 |               |
|                                   |                                 |               |
|timepoint3_S3_L001_R1_001.fastq.gz |timepoint3_S3                    |2              |
|timepoint3_S3_L001_R2_001.fastq.gz |                                 |               |
|timepoint3_S3_L002_R1_001.fastq.gz |                                 |               |
|timepoint3_S3_L002_R2_001.fastq.gz |                                 |               |
|timepoint3_S3_L003_R1_001.fastq.gz |                                 |               |
|timepoint3_S3_L003_R2_001.fastq.gz |                                 |               |
|timepoint3_S3_L004_R1_001.fastq.gz |                                 |               |
|timepoint3_S3_L004_R2_001.fastq.gz |                                 |               |

Here, `TIMEPOINT_FASTQ` = `{ 'timepoint1_S1': 0, 'timepoint2_S2': 1, 'timepoint3_S3': 2 }`

`PAM_ORIENTATION`: Orientation of the PAM relative to the spacer (`three_prime` (Cas9) or `five_prime` (Cas12)). 

`PAM_LENGTH`: Number of bases to include in the PAM.

`PAM_START`: Position relative to the spacer to begin defining the PAM. Enter `0` to start immediately adjacent to the spacer. Many CRISPR-Cas nucleases have extended PAM preferences that may make it advantageous to consider more distal positions as the "start" of the PAM. For example, for the canonical *Staphylococcus aureus* Cas9 PAM of "NNGRRT", it may be desirable to define `PAM_START` = 2, to consider a 4 nucleotide "GRRT" PAM, rather than the full 6 nucleotides.

Examples for specifying the PAM sequence:

Specifying the `PAM_LENGTH` and `PAM_START` arguments for an example spacer and `three_prime` PAM orientation: 

    5'-[SPACER][PAM]-3'    
    5'-[SPACER]CGGTA-3'
               01234  (PAM_START)
Examples:

|PAM_LENGTH|PAM_START|PAM sequence|
|:--------:|:-------:|:----------:|
|3         |0        |CGG         |
|4         |0        |CGGT        |
|3         |1        |GGT         |

Specifying the `PAM_LENGTH` and `PAM_START` arguments for an example spacer and `five_prime` PAM orientation: 

    5'-[PAM][SPACER]-3'    
    5'-CTTTA[SPACER]-3'
       43210 (PAM_START)
Examples:

|PAM_LENGTH|PAM_START|PAM sequence|
|:--------:|:-------:|:----------:|
|5         |0        |CTTTA       |
|4         |0        |TTTA        |
|4         |1        |CTTT        |

`CONTROL_RAW_COUNT_CSV`: Filepath to the "PAMDA_1_raw_counts_csv.gz" output from the randomized PAM library quality control analysis. If the quality control analysis was not performed and the untreated PAM library control sample is instead pooled with unique sample barcodes in one of the timepoint fastqs, set `CONTROL_RAW_COUNT_CSV` = `None` and instead indicate the timepoint fastq containing the control sample in the `CONTROL_SAMPLE_TIMEPOINT_FASTQ` parameter.

`CONTROL_SAMPLE`: Unique sample ID of the untreated randomized PAM library control sample.

`CONTROL_SAMPLE_TIMEPOINT_FASTQ`: Integer index of the timepoint fastq that contains the control sample. If the quality control analysis was not performed and the untreated PAM library control sample is instead pooled with unique sample barcodes in one of the timepoint fastqs, set `CONTROL_RAW_COUNT_CSV` = `None` and instead indicate the timepoint fastq containing the control sample in the `CONTROL_SAMPLE_TIMEPOINT_FASTQ` parameter. If the control were pooled in the timepoint 3 fastq, then the value of `CONTROL_SAMPLE_TIMEPOINT_FASTQ` would be `3`.

`TIMEPOINTS`: List of integer timepoints used, including the zeroth timepoint. Use a consistent unit. Default is `[0, 60, 480, 1920]` (seconds).

`MAX_PAM_LENGTH`: Integer. The maximum PAM length that will be considered, starting immediately adjacent to the spacer. PAMs of this length will be enumerated from fastq files. PAMs will later be truncated to the bases indicated by `PAM_START` and `PAM_LENGTH`. Default value is `8`. 


*The following parameters should be changed if target site design, amplicon design, or library preparation were altered from those used in Walton et al. (Science, 2020).*

`SPACERS`: Dictionary with key = spacer name, value = spacer sequence. Default spacers: `{'SPACER1':'GGGCACGGGCAGCTTGCCGG', 'SPACER2':'GTCGCCCTCGAACTTCACCT'}`
           
`P5_SAMPLE_BARCODE_START`: Integer, location of the 5' end of the P5 sample barcode. Default = `2`.
        
`P7_SAMPLE_BARCODE_START`: Integer, location of the 5' end of the P7 sample barcode. Default = `2`.

Example P5 and P7 barcodes in an example amplicon for a 3' PAM library below. Barcodes are the capitalized "NNNN" sequences, read from 5' to 3' on their respective strands:

                 P7 BARCODE                                                
    POSITION:  012
            5' --NNNN-----------------[SPACER][PAM]-----------------nnnn-- 3'
            3' --nnnn-----------------[SPACER][PAM]-----------------NNNN-- 5'
                                                        POSITION:      210
                                                                    P5 BARCODE


*Addional arguments for normalization, rate constant calculation, and heat map generation:*

`USE_TIMEPOINTS`: Default (`None`) will use all timepoints to fit rate constants. To use only a subset, provide a list specifing by the indices of the `TIMEPOINTS` list to use. For example, to use only timepoints 0 and 3, enter: `[0,3]`

`TOP_N_NORMALIZE`: For each sample, PAM enrichment is normalized to reflect no depletion. Default is to average the top `5` most enriched PAMs to define a per-timepoint normalization. To change, provide the desired integer value.

`INIT_RATE_EST`: Nonlinear rate calculations require initial rate estimates. Default is `[0.0001, 0.001, 0.01]`. To change, provide a list of floats.

`READ_SUM_MIN`: Minimum number of total raw reads for a PAM across all timepoints to calculate a rate. Otherwise rate value will be NaN. Default is `4`.

`TPS_SUM_MIN`: Minimum number of timepoints with non-zero raw read count for a PAM to calculate a rate. Otherwise rate value will be NaN. Default is `1`.


*The following options are for formatting heatmaps. The heatmap representation of PAM preference splits the bases of the PAM across the y- and x-axis. For example, an ACGT PAM could be represented with AC on the y-axis and GT on the x-axis. The following options are intended to allow customization of this representation. Rate constants of PAM targeting are used to color the heatmap where darker blue indcates faster PAM targeting.*

`PAM1_NT_RANK`: Dictionary describing the ordering of nucleotides on the heatmap y-axis. An example change might group purines (A and G) and pyrimidines (C and T) together: `{1:'A',2:'G',3:'C',4:'T'}`. Default is `{1:'A',2:'C',3:'G',4:'T'}` 

`PAM2_NT_RANK`: Dictionary describing the ordering of nucleotides on the heatmap x-axis. An example change might group purines (A and G) and pyrimidines (C and T) together: `{1:'A',2:'G',3:'C',4:'T'}`. Default is `{1:'A',2:'C',3:'G',4:'T'}`

`PAM1_INDEX_RANK` and `PAM2_INDEX_RANK` determine how the heatmap axes are defined and how positions of the PAM are ordered. The bases of the PAM will be split across the axes. Both are lists of integers. The length of each list will determine the number of bases that are represented on each axis. The position of the list corresponds to the position in the first or second half of the PAM and the integer value at each position in the list corresponds to the priority of the position. A higher integer value indicates a higher priority.

`PAM1_INDEX_RANK`: Ranking of nucleotide positions on the heatmap y-axis (first half of the PAM), represented as a list of integers (determines order of rows). Default `None` will use a default organization. To change, provide a list of integers like `[2,1]` where the first position is given a priority of 2 and the second position is given a priority of 1. If specified, the length of `PAM1_INDEX_RANK` and `PAM2_INDEX_RANK` must add to the total `PAM_LENGTH`.

`PAM2_INDEX_RANK`: Ranking of nucleotide positions on the heatmap x-axis (second half of the PAM), represented as a list of integers (determines order of columns). Default `None` will use a default organization. To change, provide a list of integers like `[2,1]` where the first position is given a priority of 2 and the second position is given a priority of 1. If specified, the length of `PAM1_INDEX_RANK` and `PAM2_INDEX_RANK` must add to the total `PAM_LENGTH`.

`AVERAGE_SPACER`: Average rates across spacers for the heatmap. Default is `True`. If `False`, separate heatmaps are generated for each spacer.

`HEATMAP_FIXED_MIN`: Minimum value of heatmap. Default is `-1.5` (log scale). Set to `False` to automatically choose the minimum value for each sample from rate values.

`HEATMAP_FIXED_MAX`: Maximum value of heatmap. Default is `-5.0` (log scale). Set to `False` to automatically choose the maximum value for each sample from rate values.

`LOG_SCALE_HEATMAP`: Plot log10(rates) or rates. `True` or `False`. Default `True` plots log10(rates).
"""
