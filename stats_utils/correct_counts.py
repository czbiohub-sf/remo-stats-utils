"""
Base class for correcting YOGO class counts
"""
import numpy as np
import numpy.typing as npt

from typing import List, Tuple
from abc import ABC, abstractmethod

from yogo.data import YOGO_CLASS_ORDERING

class CountCorrector:
    def __init__(self, inv_cmatrix: npt.NDArray, inv_cmatrix_std: npt.NDArray, parasite_ids: List[int], rbc_ids: List[int]):
        self.inv_cmatrix = inv_cmatrix
        self.inv_cmatrix_std = inv_cmatrix_std
        self.parasite_ids = parasite_ids
        self.rbc_ids = rbc_ids
    
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

    def calc_poisson_var_terms(self, raw_counts: npt.NDArray) -> npt.NDArray:
        """
        Return absolute uncertainty term of each class count based on Poisson statistics

        See remoscope manuscript for full derivation
        """
        return np.matmul(raw_counts, np.square(self.inv_cmatrix))

    def calc_deskew_var_terms(self, raw_counts: npt.NDArray) -> npt.NDArray:
        """
        Return absolute uncertainty term of each class count based on correction

        See remoscope manuscript for full derivation
        """
        return np.matmul(np.square(raw_counts), np.square(self.inv_cmatrix_std))    

    def calc_parasitemia(self, corrected_counts: npt.NDArray) -> float:
        """
        Return total parasitemia count
        """
        parasites = np.sum(corrected_counts[self.parasite_ids])
        rbcs = np.sum(corrected_counts[self.rbc_ids])
        return 0 if rbcs == 0 else parasites / rbcs

    def calc_parasitemia_rel_err(self, count_vars: npt.NDArray, parasites: float) -> float:
        """
        Return relative uncertainty of total parasitemia count

        See remoscope manuscript for full derivation
        """
        parasite_count_vars = count_vars[self.parasite_ids]

        # Compute error
        return np.inf if parasites == 0 else np.sqrt(np.sum(parasite_count_vars)) / parasites