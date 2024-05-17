"""
Base class for correcting YOGO class counts
"""

import numpy as np
import numpy.typing as npt

from typing import List, Tuple, Union
from abc import ABC, abstractmethod

from stats_utils.constants import (
    YOGO_CLASS_ORDERING,
    YOGO_CLASS_IDX_MAP,
    PARASITES_P_UL_PER_PERCENT,
)


class CountCorrector:
    def __init__(
        self,
        inv_cmatrix: npt.NDArray,
        inv_cmatrix_std: npt.NDArray,
        rbc_ids: List[int],
        parasite_ids: List[int],
    ):
        self.inv_cmatrix = inv_cmatrix
        self.inv_cmatrix_std = inv_cmatrix_std

        self.rbc_ids = rbc_ids
        self.parasite_ids = parasite_ids

    def _correct_counts(self, raw_counts: npt.NDArray):
        """
        Correct raw counts using inverse confusion matrix

        Returns list of corrected cell counts and rounds negative values to 0
        """
        corrected_counts = np.matmul(raw_counts, self.inv_cmatrix)

        # Round all negative values to 0
        corrected_counts[corrected_counts < 0] = 0

        return corrected_counts

    def _calc_count_vars(self, raw_counts: npt.NDArray) -> npt.NDArray:
        """
        Return absolute uncertainty of each class count based on deskewing and Poisson statistics

        See remoscope manuscript for full derivation
        """
        poisson_terms = self._calc_poisson_var_terms(raw_counts)
        deskew_terms = self._calc_deskew_var_terms(raw_counts)

        count_vars = poisson_terms + deskew_terms

        return count_vars

    def _calc_poisson_var_terms(self, raw_counts: npt.NDArray) -> npt.NDArray:
        """
        Return absolute uncertainty term of each class count based on Poisson statistics

        See remoscope manuscript for full derivation
        """
        return np.matmul(raw_counts, np.square(self.inv_cmatrix))

    def _calc_deskew_var_terms(self, raw_counts: npt.NDArray) -> npt.NDArray:
        """
        Return absolute uncertainty term of each class count based on correction

        See remoscope manuscript for full derivation
        """
        return np.matmul(np.square(raw_counts), np.square(self.inv_cmatrix_std))

    def _calc_parasitemia(
        self, counts: npt.NDArray, parasites: Union[None, float] = None
    ) -> float:
        """
        Return total parasitemia count
        """
        rbcs = np.sum(counts[self.rbc_ids])
        if parasites is None:
            parasites = np.sum(counts[self.parasite_ids])

        return 0 if rbcs == 0 else parasites / rbcs

    def _calc_parasitemia_rel_err(
        self,
        corrected_counts: npt.NDArray,
        count_vars: npt.NDArray,
        parasites: Union[None, float] = None,
    ) -> float:
        """
        Return relative uncertainty of total parasitemia count

        See remoscope manuscript for full derivation
        """
        if parasites is None:
            parasites = np.sum(corrected_counts[self.parasite_ids])

        parasite_count_vars = count_vars[self.parasite_ids]

        # Compute error
        return (
            np.inf
            if parasites == 0
            else np.sqrt(np.sum(parasite_count_vars)) / parasites
        )

    def _get_res_from_counts(
        self, raw_counts: npt.NDArray, units_ul_out: bool = False
    ) -> Tuple[float, float]:
        """
        Return parasitemia and 95% confidence bound based on class counts

        See remoscope manuscript for full derivation
        """
        # Correct counts
        corrected_counts = self._correct_counts(raw_counts)
        parasites = np.sum(corrected_counts[self.parasite_ids])

        # Calc_parasitemia
        parasitemia = self._calc_parasitemia(corrected_counts, parasites=parasites)

        # Get uncertainties
        count_vars = self._calc_count_vars(raw_counts)

        # Use rule of 3 if there are no parasites
        if parasites == 0:
            abs_bound = 100 * 3 / corrected_counts[YOGO_CLASS_IDX_MAP["healthy"]]
        else:
            rel_bound = 1.69 * self._calc_parasitemia_rel_err(
                corrected_counts, count_vars, parasites=parasites
            )
            abs_bound = rel_bound * parasitemia

        if units_ul_out:
            return (
                parasitemia * PARASITES_P_UL_PER_PERCENT,
                abs_bound * PARASITES_P_UL_PER_PERCENT,
            )
        else:
            return parasitemia, abs_bound

    def _get_95_confidence_bound(self, parasitemia: float, bound: float) -> List[float]:
        """
        Return 95% confidence bound on parasitemia estimate
        """
        lower_bound = max(0, parasitemia - bound)
        upper_bound = parasitemia + bound

        return [lower_bound, upper_bound]
