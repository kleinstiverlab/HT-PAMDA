<p align="center">
<a href="https://www.kleinstiverlab.org/" target="_blank"> 
    <img src="https://static.wixstatic.com/media/d66ab3_7c653ab4dc6b44319f112fefef60fafe~mv2.png/v1/fill/w_99,h_80,al_c,q_85,usm_0.66_1.00_0.01/190409%20-%20KleinstiverLabLogo-05.webp" alt="Kleinstiver Lab" /> 
</a>
</p>

<h1 align="center"> High-Throughput PAM Determination Assay (HT-PAMDA)</h1>

# Overview

HT-PAMDA is used to comprehensively profile the protospacer-adjacent motif (PAM) preferences of CRISPR-Cas nucleases. This repository contains code to perform analysis of data that has been generated using the HT-PAMDA method as described in [Walton et al. (Science, 2020)](https://science.sciencemag.org/content/368/6488/290).

|SpCas9|SpRY |
|:----:|:---:|
|<img src="https://github.com/kleinstiverlab/HT-PAMDA/blob/master/figures/github_example_figures/PAMDA_HEATMAP_sample_01_WT_SpCas9.png" alt="PAM preference of SpCas9" width="450"/>|<img src="https://github.com/kleinstiverlab/HT-PAMDA/blob/master/figures/github_example_figures/PAMDA_HEATMAP_sample_03_SpRY.png" alt="PAM preference of SpRY" width="450"/>|

## Installation

Download the repository (on Github use the green "Clone or download" button and select "Download Zip"). Expand the zip file in a convenient location.

In Terminal, go to the repository directory and create a Python 3 virtual environment:

```
python3 -m venv venv
```

This will create a virtual environment called `venv`. Next, the commands in `install.sh` will add the required packages to the `venv` virtual environment. The installation should take no more than a few minutes.

```
sh install.sh
```

The required dependencies are now installed to the virtual environment. Activate the virtual environment:

```
source venv/bin/activate
```

Perform all steps of the analysis within this virtual environment.

## Folder structure

The repository has the following folder structure:
* `barcode_csv`: directory containing example `BARCODE_CSV` files indicating sample barcoding
* `code`: directory containing the HT-PAMDA scripts
* `fastqs`: directory containing example fastq files
* `figures`: directory where figures generated from the analysis are saved
* `output`: directory where the HT-PAMDA output files are saved

# Analysis

There are two main sections of the analysis:
1) **Randomized PAM library quality control**: Assess the distribution of PAM sequences in an untreated randomized library control.
2) **HT-PAMDA data analysis**: Analyze HT-PAMDA data to define PAM preferences of CRISPR-Cas nucleases.

Example data has been provided to run both sections of the analysis.

## Randomized PAM library quality control

The following function performs the complete randomized PAM library quality control analysis:

**library_QC**: Takes fastq files as input and determines the distribution of PAMs in the randomized PAM library. Raw read counts are saved to the `output` folder. Figures describing the library distribution are saved to the `figures` folder.

**To run the library QC analysis:**

Navigate to the `code` directory, open the `library_QC_inputs.py` file, and enter the required parameters. Run the analysis:
```
python3 library_QC.py
```

On the provided example data, this step should only take a few seconds.

## Input parameters
### Required parameters

*These parameters are likely to change with any given implementation of the assay. Please enter values for all the following parameters. Leave default values to run an example dataset.*

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
|control_S1_L001_R1_001.fastq.gz |control_S1      |
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

### Additional parameters

*The following arguments have default values consistent with our standard HT-PAMDA implementation. They are defined here with the intention of providing flexibility for changes to assay design. These values do not need to be changed to run the analysis for our standard HT-PAMDA protocol. Changes to assay design will likely require changes to these parameters.*

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

## HT-PAMDA data analysis

The analysis consists of three steps plus figure generation to represent PAM preference as a heatmap. 

The following function performs the complete analysis:

**PAMDA_complete**: Run the complete pipeline from fastq files to rate constants and heatmaps of each sample. Takes fastq files and CSV with sample barcodes as input. Outputs CSV files for raw counts, normalized counts, and rates. Also outputs heatmaps as PDF and CSV files.

|function          |input                               |outputs                    |
|:-----------------|:-----------------------------------|:--------------------------|
|PAMDA_complete    | fastq files                        |`PAMDA_1_raw_counts.csv.gz`|
|                  |                                    |`PAMDA_2_norm_counts.csv`  |  
|                  |                                    |`PAMDA_3_rates.csv`        |
|                  |                                    | heatmaps                  |
 


Or the following functions can be used to perform these steps separately:

1. **fastq2count**: Convert fastq files to raw read counts for each sample, spacer, and PAM. Takes fastq files and CSV with sample barcodes as input. Outputs a compressed CSV file with raw read counts.
2. **rawcount2normcount**: Calculate normalized read counts from raw read counts. Takes a compressed CSV with raw counts as input. Outputs CSV file with normalized read counts.
3. **normcount2rate**: Calculate rate constants from normalized read counts. Takes CSV with normalized counts as input. Outputs CSV file with rate constants.
4. **rate2heatmap**: Plot heat maps representing the PAM profiles of each sample. Takes CSV with rate constants as input. Outputs heatmaps as PDF and CSV files.

|function          |input                      |output                      |
|:-----------------|:--------------------------|:---------------------------|
|fastq2count       |fastq files                |`PAMDA_1_raw_counts.csv.gz` |
|rawcount2normcount|`PAMDA_1_raw_counts.csv.gz`|`PAMDA_2_norm_counts.csv`   |
|normcount2rate    |`PAMDA_2_norm_counts.csv`  |`PAMDA_3_rates.csv`         |
|rate2heatmap      |`PAMDA_3_rates.csv`        |heatmaps                    |


**To run the HT-PAMDA analysis:**

Navigate to the `code` directory, open the `inputs.py` file, and enter the required parameters. Run the complete analysis:
```
python3 PAMDA_complete.py
```
Or run individual steps of the analysis:
```
python3 fastq2count.py
python3 rawcount2normcount.py
python3 normcount2rate.py
python3 rate2heatmap.py
```

The analysis can be run on the provided example data by leaving all parameters at their default values. On the example data, the complete analysis should take about a minute.

Additionally, barcodes for all samples from [Walton et al. (Science, 2020)](https://science.sciencemag.org/content/368/6488/290) are available (Table S7 - PAMDA data summary) and can be used to analyze HT-PAMDA data uploaded to the NCBI sequence read archive (SRA) under [BioProject ID: PRJNA605711](http://www.ncbi.nlm.nih.gov/bioproject/605711).

## Input parameters
### Required parameters

*These parameters are likely to change with any given implementation of the assay. Please enter values for all the following parameters. Leave default values to run an example dataset.*

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

### Additional parameters

*The following arguments have default values consistent with our standard HT-PAMDA implementation. They are defined here with the intention of providing flexibility for changes to assay design, analysis parameters, and figure generation. These values do not need to be changed to run the analysis for our standard HT-PAMDA protocol. Changes to assay design will likely require changes to these parameters.*


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
