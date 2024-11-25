"""
Correct for skew in YOGO classification, where skew is represented by the confusion matrix.

The correcting transformation matrix used by the base class CountCorrector is the inverse
confusion matrix. This matrix is computed from the average of k partitions of training data,
inspired by k-fold validation.

The error in each inverse confusion matrix term is individually computed as the standard deviation
across k confusion matrices.
"""

import numpy as np
import numpy.typing as npt

from typing import Tuple
from pathlib import Path

from stats_utils.constants import (
    RBC_CLASS_IDS,
    ASEXUAL_PARASITE_CLASS_IDS,
    DATA_DIR,
    CMATRIX_MEAN_SUFFIX,
    INV_CMATRIX_STD_SUFFIX,
)
from stats_utils.corrector import CountCorrector


class CountDeskewer(CountCorrector):
    def __init__(self, model_name: str):
        """
        Initialize count deskewer

        Input(s)
        - model_name:
            Include name and number (eg. "frightful-wendigo-1931")
        """
        cmatrix_mean_dir = str(
            DATA_DIR / model_name / (model_name + CMATRIX_MEAN_SUFFIX)
        )
        inv_cmatrix_std_dir = str(
            DATA_DIR / model_name / (model_name + INV_CMATRIX_STD_SUFFIX)
        )

        # Check that cmatrix data exists
        if not Path(cmatrix_mean_dir).is_file():
            raise FileNotFoundError(
                f"Could not find confusion matrix mean for {model_name} ({cmatrix_mean_dir})"
            )
        if not Path(inv_cmatrix_std_dir).is_file():
            raise FileNotFoundError(
                f"Could not find inverse confusion matrix std for {model_name} ({inv_cmatrix_std_dir})"
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

    def calc_parasitemia(self, counts: npt.NDArray) -> float:
        """
        Wrapper for base class method _calc_parasitemia()

        Input(s)
        - counts:
            Cell counts, formatted as 7x1 array with all YOGO classes
        """

        return self._calc_parasitemia(counts)

    def get_res_from_counts(
        self, raw_counts: npt.NDArray, units_ul_out: bool = False
    ) -> Tuple[float, float]:
        """
        Wrapper for base class method _get_res_from_counts()

        Input(s)
        - counts:
            Raw cell counts, formatted as 7x1 array with all YOGO classes
        - units_ul_out (optional):
            True to return parasitemia in parasitemia/uL
            False to return parasitemia in % (default)
        """

        return self._get_res_from_counts(raw_counts, units_ul_out=units_ul_out)
