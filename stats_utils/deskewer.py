"""
YOGO class deskewing using inverse confusion matrix

Based on cultured lab data
"""

import numpy as np

from stats_utils.constants import (
    RBC_CLASS_IDS,
    ASEXUAL_PARASITE_CLASS_IDS,
    YOGO_CMATRIX_MEAN_DIR,
    YOGO_INV_CMATRIX_STD_DIR,
)
from stats_utils.correct_counts import CountCorrector


class CountDeskewer(CountCorrector):
    def __init__(self):
        # Load confusion matrix data
        norm_cmatrix = np.load(YOGO_CMATRIX_MEAN_DIR)
        inv_cmatrix_std = np.load(YOGO_INV_CMATRIX_STD_DIR)

        # Compute inverse
        inv_cmatrix = np.linalg.inv(norm_cmatrix)

        super(CountDeskewer, self).__init__(
            inv_cmatrix,
            inv_cmatrix_std,
            RBC_CLASS_IDS,
            ASEXUAL_PARASITE_CLASS_IDS,
        )
