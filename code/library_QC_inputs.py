# REQUIRED ARGUMENTS
RUN_NAME = 'example_PAM_library'
CONTROL_BARCODE_CSV = 'barcode_csv/example_PAM_library_barcodes.csv'
CONTROL_FASTQ_DIR = 'fastqs/example_PAM_library'
CONTROL_FASTQ = 'expRW086_pool_03_S3'
PAM_ORIENTATION = 'three_prime' # or 'five_prime'
PAM_LENGTH = 4
PAM_START = 0

# ADDITIONAL ARGUMENTS (with default values entered)
MAX_PAM_LENGTH = 8
SPACERS = {'SPACER1':'GGGCACGGGCAGCTTGCCGG', 
           'SPACER2':'GTCGCCCTCGAACTTCACCT'}
P5_SAMPLE_BARCODE_START = 2
P7_SAMPLE_BARCODE_START = 2

"""
Argument definitions:

`RUN_NAME`: Provide a name for the run (no spaces or special characters).

`BARCODE_CSV`: Filepath to a CSV with the control barcode information and the following headers:

|sample	            |description|P5_sample_barcode|P7_sample_barcode|
|:------------------|:----------|:---------------:|:---------------:|
|unique_sample_id_01|control    |ATGC             |TCGC             |

Barcodes should be provided 5' to 3' as they appear in the primer sequences that encode them. Avoid the use of special characters.

`FASTQ_DIR`: Folder directory contianing fastq files.

`CONTROL_FASTQ`: Abbreviated name of the fastq files for the library. Fastq names must be abbreviated as follows:

|fastq complete file names       |abbreviated name|
|:-------------------------------|:--------------:|
|control_S1_L001_R1_001.fastq.gz |timepoint1_S1   |
|control_S1_L001_R2_001.fastq.gz |                |
|control_S1_L002_R1_001.fastq.gz |                |
|control_S1_L002_R2_001.fastq.gz |                |
|control_S1_L003_R1_001.fastq.gz |                |
|control_S1_L003_R2_001.fastq.gz |                |
|control_S1_L004_R1_001.fastq.gz |                |
|control_S1_L004_R2_001.fastq.gz |                |

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


`CONTROL_SAMPLE`: Unique sample ID of the untreated randomized PAM library control sample.

`MAX_PAM_LENGTH`: Integer. The maximum PAM length that will be considered, starting immediately adjacent to the spacer. PAMs of this length will be enumerated from fastq files. PAMs will later be truncated to the bases indicated by `PAM_START` and `PAM_LENGTH`. Default value is `8`. 

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
"""
