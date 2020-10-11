import os
import sys
os.chdir('../')
sys.path.append('code')
import PAMDA
from library_QC_inputs import *

# calculate rate constants from normalized read counts
PAMDA.library_QC(RUN_NAME, 
                 CONTROL_BARCODE_CSV,
                 CONTROL_FASTQ_DIR,
                 CONTROL_FASTQ,
                 PAM_ORIENTATION,
                 PAM_LENGTH,
                 PAM_START,
                 MAX_PAM_LENGTH, 
                 SPACERS,
                 P5_SAMPLE_BARCODE_START,
                 P7_SAMPLE_BARCODE_START)