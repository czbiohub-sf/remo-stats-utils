"""
Abstract base class for correcting YOGO class counts.

The class count correction is generalizable and is computed by multiplying a row vector of 
raw cell counts by a correcting transformation matrix. After correction, the parasitemia and
corresponding 95% confidence bound is computed.

The 95% confidence bound is computed by propogating errors in the linear expansion of this matrix
multiplication. The error in each raw count is based on Poisson and the error from each confusion
matrix term is provided as a matrix input argument. In case 0 parasites are counted, the error is computed
using the rule of three.
"""


import numpy as np
import numpy.typing as npt

from typing import List, Tuple, Union
from abc import ABC, abstractmethod

from stats_utils.constants import RBCS_P_UL


class CountCorrector(ABC):
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

    @abstractmethod
    def get_res_from_counts(self, raw_counts: npt.NDArray, units_ul_out: bool = False
    ) -> Tuple[float, float]:
        pass

    @abstractmethod
    def calc_parasitemia(self, counts: npt.NDArray) -> float:
        pass

    def _correct_counts(self, raw_counts: npt.NDArray):
        """
        Correct raw counts using transformation matrix

        Returns list of corrected cell counts and rounds negative values to 0
        """
        corrected_counts = np.matmul(raw_counts, self.inv_cmatrix)

        # Round all negative values to 0
        corrected_counts[corrected_counts < 0] = 0

        return corrected_counts

    def _calc_count_vars(self, raw_counts: npt.NDArray) -> npt.NDArray:
        """
        Return absolute uncertainty of each class count based on correction matrix error
        and Poisson statistics
        """
        poisson_terms = self._calc_poisson_var_terms(raw_counts)
        deskew_terms = self._calc_deskew_var_terms(raw_counts)

        count_vars = poisson_terms + deskew_terms

        return count_vars

    def _calc_poisson_var_terms(self, raw_counts: npt.NDArray) -> npt.NDArray:
        """
        Return absolute uncertainty of each class count based on Poisson statistics
        """
        return np.matmul(raw_counts, np.square(self.inv_cmatrix))

    def _calc_deskew_var_terms(self, raw_counts: npt.NDArray) -> npt.NDArray:
        """
        Return absolute uncertainty of each class count based on correction
        """
        return np.matmul(np.square(raw_counts), np.square(self.inv_cmatrix_std))

    def _calc_parasitemia(
        self, counts: npt.NDArray, parasites: Union[None, float] = None
    ) -> float:
        """
        Return total parasitemia as fractional percentage

        Input(s)
        - counts:
            Cell counts, formatted as a row vector
        - parasites (optional):
            Input parasite count if it has been previously computed
        """
        rbcs = np.sum(counts[self.rbc_ids])
        if parasites is None:
            parasites = np.sum(counts[self.parasite_ids])

        return 0 if rbcs == 0 else parasites / rbcs

    def _calc_parasites_abs_std(
        self,
        corrected_counts: npt.NDArray,
        count_vars: npt.NDArray,
        parasites: Union[None, float] = None,
    ) -> float:
        """
        Return absolute uncertainty of corrected parasitemia

        Input(s)
        - corrected_counts:
            Corrected cell counts, formatted as a row vector
        - count_vars:
            Absolute uncertainty of each class count, including Poisson error and
            confusion matrix error
        - parasites (optional):
            Input parasite count if it has been previously computed
        """
        if parasites is None:
            parasites = np.sum(corrected_counts[self.parasite_ids])

        parasite_count_vars = count_vars[self.parasite_ids]
        return np.sqrt(np.sum(parasite_count_vars))

    def _get_res_from_counts(
        self, raw_counts: npt.NDArray, units_ul_out: bool = False
    ) -> Tuple[float, float]:
        """
        Return parasitemia and 95% confidence bounds

        Input(s)
        - raw_counts:
            Raw cell counts, formatted as a row vector
        - units_ul_out (optional):
            True to return parasitemia in parasitemia/uL
            False to return parasitemia in % (default)
        """
        # Correct counts
        corrected_counts = self._correct_counts(raw_counts)
        parasites = np.sum(corrected_counts[self.parasite_ids])
        rbcs = np.sum(corrected_counts[self.rbc_ids])

        # Calc_parasitemia
        parasitemia = 100.0 * self._calc_parasitemia(
            corrected_counts, parasites=parasites
        )

        # Get uncertainties
        count_vars = self._calc_count_vars(raw_counts)

        # Use rule of 3 if there are no parasites
        if parasites == 0:
            parasites_95_conf_bounds = 3.0
            parasitemia_95_conf_bounds = (
                100 * parasites_95_conf_bounds / rbcs
            )  # unit: %
        else:
            parasites_95_conf_bounds = 1.69 * self._calc_parasites_abs_std(
                corrected_counts, count_vars, parasites=parasites
            )
            parasitemia_95_conf_bounds = (
                100 * parasites_95_conf_bounds / rbcs
            )  # unit: %

        if units_ul_out:
            return (
                parasitemia * RBCS_P_UL / 100.0,  # unit: parasitemia / uL
                parasitemia_95_conf_bounds
                * RBCS_P_UL
                / 100.0,  # unit: parasitemia / uL
            )
        else:
            return parasitemia, parasitemia_95_conf_bounds

    def _get_95_confidence_bound(self, parasitemia: float, bound: float) -> List[float]:
        """
        Return 95% confidence bounds on parasitemia estimate
        """
        lower_bound = max(0, parasitemia - bound)
        upper_bound = parasitemia + bound

        return [lower_bound, upper_bound]
