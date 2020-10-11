# REQUIRED ARGUMENTS
RUN_NAME = 'example_PAMDA_data'
BARCODE_CSV = 'barcode_csv/example_PAMDA_barcodes.csv'
FASTQ_DIR = 'fastqs/example_PAMDA_data'
TIMEPOINT_FASTQ = {
                   'expRW086_pool_10_S10': 0,
                   'expRW086_pool_11_S11': 1, 
                   'expRW086_pool_12_S12': 2
                  }
PAM_ORIENTATION = 'three_prime'
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