# REQUIRED ARGUMENTS
RUN_NAME = 'example_PAM_library'
CONTROL_BARCODE_CSV = 'barcode_csv/example_PAM_library_barcodes.csv'
CONTROL_FASTQ_DIR = 'fastqs/example_PAM_library'
CONTROL_FASTQ = 'expRW086_pool_03_S3'
PAM_ORIENTATION = 'three_prime'
PAM_LENGTH = 4
PAM_START = 0

# ADDITIONAL ARGUMENTS (with default values entered)
MAX_PAM_LENGTH = 8
SPACERS = {'SPACER1':'GGGCACGGGCAGCTTGCCGG', 
           'SPACER2':'GTCGCCCTCGAACTTCACCT'}
P5_SAMPLE_BARCODE_START = 2
P7_SAMPLE_BARCODE_START = 2