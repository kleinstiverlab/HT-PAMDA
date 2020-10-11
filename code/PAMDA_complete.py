import os
import sys
os.chdir('../')
sys.path.append('code')
import PAMDA
from inputs import *

PAMDA.PAMDA_complete(RUN_NAME,
                     BARCODE_CSV,
                     FASTQ_DIR,
                     TIMEPOINT_FASTQ,
                     PAM_ORIENTATION,
                     PAM_LENGTH,
                     PAM_START,
                     CONTROL_RAW_COUNT_CSV,
                     CONTROL_SAMPLE,
                     CONTROL_SAMPLE_TIMEPOINT_FASTQ = CONTROL_SAMPLE_TIMEPOINT_FASTQ,
                     TIMEPOINTS = TIMEPOINTS,
                     MAX_PAM_LENGTH = MAX_PAM_LENGTH,
                     SPACERS = SPACERS,
                     P5_SAMPLE_BARCODE_START = P5_SAMPLE_BARCODE_START,
                     P7_SAMPLE_BARCODE_START = P7_SAMPLE_BARCODE_START,
                     USE_TIMEPOINTS = USE_TIMEPOINTS,
                     TOP_N_NORMALIZE = TOP_N_NORMALIZE,
                     INIT_RATE_EST = INIT_RATE_EST,
                     READ_SUM_MIN = READ_SUM_MIN,
                     TPS_SUM_MIN = TPS_SUM_MIN,
                     PAM1_NT_RANK = PAM1_NT_RANK,
                     PAM2_NT_RANK = PAM2_NT_RANK,
                     PAM1_INDEX_RANK = PAM1_INDEX_RANK,
                     PAM2_INDEX_RANK = PAM2_INDEX_RANK,
                     AVERAGE_SPACER = AVERAGE_SPACER,
                     HEATMAP_FIXED_MIN = HEATMAP_FIXED_MIN,
                     HEATMAP_FIXED_MAX = HEATMAP_FIXED_MAX,
                     LOG_SCALE_HEATMAP = LOG_SCALE_HEATMAP)