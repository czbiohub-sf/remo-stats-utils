"""
Base class for correcting YOGO class counts
"""
import numpy as np
import numpy.typing as npt

from typing import List, Tuple
from abc import ABC, abstractmethod

from yogo.data import YOGO_CLASS_ORDERING

class Corrector:
    def __init__(self):
        pass
    
    def correct_counts(self, counts: npt.NDArray, matrix: npt.NDArray):
        """
        Correct raw counts using inverse confusion matrix

        Returns list of corrected cell counts and rounds negative values to 0
        """
        corrected_counts = np.matmul(counts, matrix)

        # Round all negative values to 0
        corrected_counts[corrected_counts < 0] = 0

        return corrected_counts
    
    def calc_count_vars(
        self, raw_counts : npt.NDArray
    ) -> npt.NDArray:
        """
        Return absolute uncertainty of each class count based on deskewing and Poisson statistics

        See remoscope manuscript for full derivation
        """
        poisson_terms = self.calc_poisson_var_terms(raw_counts)
        deskew_terms = self.calc_deskew_var_terms(raw_counts)

        count_vars = poisson_terms + deskew_terms

        return count_vars

    def calc_poisson_var_terms(self, raw_counts: npt.NDArray, inv_cmatrix: npt.NDArray) -> npt.NDArray:
        """
        Return absolute uncertainty term of each class count based on Poisson statistics

        See remoscope manuscript for full derivation
        """
        return np.matmul(raw_counts, np.square(inv_cmatrix))

    def calc_deskew_var_terms(self, raw_counts: npt.NDArray, inv_cmatrix_std: npt.NDArray) -> npt.NDArray:
        """
        Return absolute uncertainty term of each class count based on correction

        See remoscope manuscript for full derivation
        """
        return np.matmul(np.square(raw_counts), np.square(inv_cmatrix_std))    

    def calc_parasitemia(self, corrected_counts: npt.NDArray, parasite_ids: List[int], rbc_ids: List[int]) -> float:
        """
        Return total parasitemia count
        """
        parasites = np.sum(corrected_counts[parasite_ids])
        rbcs = np.sum(corrected_counts[rbc_ids])
        return 0 if rbcs == 0 else parasites / rbcs

    def calc_parasitemia_rel_err(self, count_vars: npt.NDArray, parasite_ids: List[int], parasites: float) -> float:
        """
        Return relative uncertainty of total parasitemia count

        See remoscope manuscript for full derivation
        """
        parasite_count_vars = count_vars[parasite_ids]

        # Compute error
        return np.inf if parasites == 0 else np.sqrt(np.sum(parasite_count_vars)) / parasites