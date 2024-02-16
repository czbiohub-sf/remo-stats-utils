"""
YOGO class deskewing using inverse confusion matrix

Based on cultured lab data
"""

import numpy as np

from stats_utils.constants import (
    RBC_CLASS_IDS,
    ASEXUAL_PARASITE_CLASS_IDS,
    DATA_DIR,
    CMATRIX_MEAN_SUFFIX,
    INV_CMATRIX_STD_SUFFIX,
)
from stats_utils.correct_counts import CountCorrector


class CountDeskewer(CountCorrector):
    def __init__(self, model_name: str):
        """
        Initialize count deskewer

        Input(s)
        - model_name: 
            Include name and number (eg. "frightful-wendigo-1931")
        """
        cmatrix_mean_dir = str(DATA_DIR / model_name / (model_name + CMATRIX_MEAN_SUFFIX))
        inv_cmatrix_std_dir = str(DATA_DIR / model_name / (model_name + INV_CMATRIX_STD_SUFFIX))

        # Check that cmatrix data exists
        if not cmatrix_mean_dir.is_dir():
            raise FileNotFoundError(
                f"Could not find {model_name} confusion matrix mean file {cmatrix_mean_dir}"
            )
        if not inv_cmatrix_std_dir.is_dir():
            raise FileNotFoundError(
                f"Could not find {model_name} inverse confusion matrix std file {inv_cmatrix_std_dir}"
            )

        # Load confusion matrix data
        norm_cmatrix = np.load(cmatrix_mean_dir)
        inv_cmatrix_std = np.load(inv_cmatrix_std_dir)

        # Compute inverse
        inv_cmatrix = np.linalg.inv(norm_cmatrix)

        super(CountDeskewer, self).__init__(
            inv_cmatrix,
            inv_cmatrix_std,
            RBC_CLASS_IDS,
            ASEXUAL_PARASITE_CLASS_IDS,
        )
