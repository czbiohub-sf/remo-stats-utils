"""
YOGO class deskewing using inverse confusion matrix

Based on cultured lab data
"""
import numpy as np
import numpy.typing as npt

from typing import List, Tuple

from yogo.data import YOGO_CLASS_ORDERING
from stats_utils.constants import (
    YOGO_CLASS_IDX_MAP,
    RBC_CLASS_IDS,
    ASEXUAL_PARASITE_CLASS_IDS,
    YOGO_CMATRIX_MEAN_DIR,
    YOGO_INV_CMATRIX_STD_DIR,
    PERCENT_2_PARASITES_PER_UL,
)


class Deskewer:
    def __init__(self):
        self.matrix_dim = len(YOGO_CLASS_ORDERING)

        # Load confusion matrix data
        norm_cmatrix = np.load(YOGO_CMATRIX_MEAN_DIR)
        self.inv_cmatrix_std = np.load(YOGO_INV_CMATRIX_STD_DIR)

        # Compute inverse
        self.inv_cmatrix = np.linalg.inv(norm_cmatrix)

    def calc_res(
        self, raw_counts: npt.NDArray, units_ul: bool=False 
    ) -> Tuple[float, float, npt.NDArray]:
        """
        Return parasitemia, 95% confidence bound, and deskewed counts
        See remoscope manuscript for full derivation

        95% confidence interval can be defined as
            lower_bound = max(0, parasitemia - bound)
            upper_bound = min(1, parasitemia + bound)
        """
        # Deskew
        deskewed_counts = self.calc_deskewed_counts(raw_counts)
        
        # Get uncertainties
        count_vars = self.calc_class_count_vars(raw_counts, deskewed_counts)

        # Use rule of 3 if there are no parasites
        if parasitemia == 0:
            bound = 3 / deskewed_counts[YOGO_CLASS_IDX_MAP["healthy"]]
        else:
            bound = 1.69 * self.calc_parasitemia_rel_err(count_vars, deskewed_counts)

        if units_ul:
            return parasitemia * PERCENT_2_PARASITES_PER_UL, bound * PERCENT_2_PARASITES_PER_UL, 
        else:
            return parasitemia, bound, deskewed_counts
        
    def calc_deskewed_counts(self, raw_counts: npt.NDArray) -> npt.NDArray:
        """
        Deskew raw counts using inverse confusion matrix. Optional parameters

        Returns list of deskewed cell counts. that are whole number integers (ie. no negative vals)
        """
        deskewed_floats = np.matmul(raw_counts, self.inv_cmatrix)
        # Round all negative values to 0
        deskewed_floats[deskewed_floats < 0] = 0

        return deskewed_floats
    def calc_class_count_vars(
        self, raw_counts : npt.NDArray, deskewed_counts: npt.NDArray
    ) -> npt.NDArray:
        """
        Return absolute uncertainty of each class count based on deskewing and Poisson statistics
        See remoscope manuscript for full derivation
        """
        poisson_terms = self.calc_poisson_count_var_terms(raw_counts)
        deskew_terms = self.calc_deskew_count_var_terms(raw_counts)

        class_vars = poisson_terms + deskew_terms

        return class_vars

    def calc_poisson_count_var_terms(self, raw_counts: npt.NDArray) -> npt.NDArray:
        """
        Return absolute uncertainty term of each class count based on Poisson statistics
        See remoscope manuscript for full derivation
        """
        return np.matmul(raw_counts, np.square(self.inv_cmatrix))

    def calc_deskew_count_var_terms(self, raw_counts: npt.NDArray) -> npt.NDArray:
        """
        Return absolute uncertainty term of each class count based on deskewing
        See remoscope manuscript for full derivation
        """
        # Commented out for now because int overflow should not be an issue with bug fixes
        # # TODO does division by 0 cause error?
        # RBC_count = np.sum(raw_counts[RBC_CLASS_IDS])

        # # Use ratio of class relative to RBC count to avoid overflow
        # class_ratios = raw_counts / RBC_count
        # unscaled_err = np.matmul(np.square(class_ratios), np.square(self.inv_cmatrix_std))

        # return unscaled_err * RBC_count **2

        return np.matmul(np.square(raw_counts), np.square(self.inv_cmatrix_std))

    def calc_parasitemia(self, deskewed_counts: npt.NDArray) -> float:
        """
        Return total parasitemia count
        """
        parasites = np.sum(deskewed_counts[ASEXUAL_PARASITE_CLASS_IDS])
        RBCs = np.sum(deskewed_counts[RBC_CLASS_IDS])
        return 0 if RBCs == 0 else parasites / RBCs

    def calc_parasitemia_rel_err(self, raw_counts: npt.NDArray) -> float:
        """
        Return relative uncertainty of total parasitemia count
        See remoscope manuscript for full derivation
        """
        deskewed_counts = self.calc_deskewed_counts(raw_counts)
        count_vars = self.calc_class_count_vars(raw_counts, deskewed_counts)

        # Filter for parasite classes only
        parasite_count_vars = count_vars[ASEXUAL_PARASITE_CLASS_IDS]
        parasite_count = np.sum(deskewed_counts[ASEXUAL_PARASITE_CLASS_IDS])

        # Compute error
        return np.inf if parasite_count == 0 else np.sqrt(np.sum(parasite_count_vars)) / parasite_count