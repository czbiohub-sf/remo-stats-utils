"""
Class compensation using y = mx + b fit, where m and b are computed by comparing Remoscope
parasitemia estimates with clinical PCR values.

The correcting transformation matrix used by the base class CountCorrector is a matrix 
representation of y = mx + b. Before matrix multiplication, the 7x1 row vector of raw class
counts is simplified to a 2x1 vector [healthy, parasites], where parasites include all
asexual life stages.

The error in each transformation matrix term is based on the standard deviation of the m and b
values from the linear fit.
"""


import numpy as np
import numpy.typing as npt
import pandas as pd

from typing import Tuple, Union
from pathlib import Path

from stats_utils.constants import (
    DATA_DIR,
    CULTURED_COMPENSATION_SUFFIX1,
    CLINICAL_COMPENSATION_SUFFIX1,
    NO_HEATMAPS_SUFFIX2,
    W_HEATMAPS_SUFFIX2,
    RBCS_P_UL,
    YOGO_CLASS_IDX_MAP,
    ASEXUAL_PARASITE_CLASS_IDS,
)
from stats_utils.corrector import CountCorrector


class CountCompensator(CountCorrector):
    def __init__(
        self,
        model_name: str,
        conf_thresh: float,
        clinical: bool = True,
        heatmaps: bool = False,
        skip: bool = False,
    ):
        """
        Initialize count compensator

        Input(s)
        - model_name:
            Include name and number (eg. "frightful-wendigo-1931")
        - conf_thresh:
            The confidence threshold defines the required model confidence to accept a class call.
            Compensation using y = mx + b fit was computed for multiple thresholds, input the
            desired threshold to extract the corresponding m and b values
        - clinical (optional):
            True for clinical Uganda data (default)
            False for cultured lab data
        - heatmaps (optional):
            True for heatmap corrected data
            False otherwise (default)
        - skip (optional):
            True to skip compensation and compute parasitemia and error from raw data only.
            False to proceed with compensation (default)
        """

        self.conf_thresh = conf_thresh
        if skip:
            m = 1.0
            b = 0.0
            cov_m = 0.0
            cov_b = 0.0
        else:
            # Generate directory for compensation metrics csv
            if clinical:
                suffix1 = CLINICAL_COMPENSATION_SUFFIX1
            else:
                suffix1 = CULTURED_COMPENSATION_SUFFIX1
            if heatmaps:
                suffix2 = W_HEATMAPS_SUFFIX2
            else:
                suffix2 = NO_HEATMAPS_SUFFIX2
            compensation_csv_dir = str(
                DATA_DIR / model_name / (model_name + suffix1 + suffix2)
            )

            # Check that compensation metrics csv exists
            if not Path(compensation_csv_dir).is_file():
                raise FileNotFoundError(
                    f"Could not find compensation metrics for {model_name} ({compensation_csv_dir})"
                )

            m, b, cov_m, cov_b = self._get_fit_metrics(compensation_csv_dir)
            # ignore negative b
            if b < 0:
                b = 0
                cov_b = 0

        inv_cmatrix = self._get_matrix(m, b)
        inv_cmatrix_std = self._get_matrix_std(cov_m, cov_b)

        rbc_ids = [0, 1]
        parasite_ids = [1]

        super(CountCompensator, self).__init__(
            inv_cmatrix,
            inv_cmatrix_std,
            rbc_ids,
            parasite_ids,
        )

    def _get_fit_metrics(
        self, compensation_csv_dir: str
    ) -> Tuple[float, float, float, float]:
        """
        Extract fit metrics for a given confidence threshold from .csv

        Returns fit metrics as (m, b, cov_m, cov_b)
        """
        df = pd.read_csv(compensation_csv_dir, dtype=np.float64)
        row = df.loc[df["conf_val"] == self.conf_thresh]

        # Adjust b for parasitemia fractional percentage instead of parasites per uL
        fit_b = row["fit_b"].item() / RBCS_P_UL
        cov_b = row["cov_b"].item() / RBCS_P_UL

        return row["fit_m"].item(), fit_b, row["cov_m"].item(), cov_b

    def _get_matrix(self, m: float, b: float) -> npt.NDArray:
        """
        Return transformation matrix equivalent of y = mx + b
        """
        M12 = -b
        M22 = (1 / m) - b

        M11 = 1 - M12
        M21 = 1 - M22

        return np.array([[M11, M12], [M21, M22]])

    def _get_matrix_std(self, cov_m: float, cov_b: float) -> npt.NDArray:
        """
        Return standard deviation of transformation matrix
        """
        M11 = cov_b
        M12 = cov_b

        M21 = np.sqrt(cov_b**2 + cov_m**2)
        M22 = np.sqrt(cov_b**2 + cov_m**2)

        return np.array([[M11, M12], [M21, M22]])

    def _reformat_7x1_to_2x1(self, counts: npt.NDArray) -> npt.NDArray:
        """
        Reformats a 7x1 np array with all YOGO classes into a 2x1 np array
        with [healthy, parasites] only
        """
        healthy = counts[YOGO_CLASS_IDX_MAP["healthy"]]
        parasites = np.sum(counts[ASEXUAL_PARASITE_CLASS_IDS])

        return np.asarray([healthy, parasites])

    def calc_parasitemia(self, counts: npt.NDArray) -> float:
        """
        Wrapper for base class method _calc_parasitemia()

        Reformats 7x1 array into required 2x1 array before computing statistics

        Input(s)
        - counts:
            Cell counts, formatted as 7x1 array with all YOGO classes
        """

        reformatted_counts = self._reformat_7x1_to_2x1(counts)
        return self._calc_parasitemia(reformatted_counts)

    def get_res_from_counts(
        self, raw_counts: npt.NDArray, units_ul_out: bool = False
    ) -> Tuple[float, float]:
        """
        Wrapper for base class method _get_res_from_counts()

        Reformats 7x1 array into required 2x1 array before computing statistics

        Input(s)
        - counts:
            Raw cell counts, formatted as 7x1 array with all YOGO classes
        - units_ul_out (optional):
            True to return parasitemia in parasitemia/uL
            False to return parasitemia in % (default)
        """

        reformatted_counts = self._reformat_7x1_to_2x1(raw_counts)
        return self._get_res_from_counts(reformatted_counts, units_ul_out=units_ul_out)

    def get_95_bound_and_compensation_from_parasitemia(
        self,
        raw_parasitemia: float,
        rbcs: Union[int, float],
        units_ul_in: bool,
        units_ul_out: bool = False,
    ) -> Tuple[float, float]:
        """
        Return compensated parasitemia and 95% confidence bounds based on raw parasitemia and
        total rbcs

        Input(s)
        - raw_parasitemia:
            Uncorrected parasitemia estimate, given as fractional percentage
        - rbcs:
            RBC count, includes healthy cells and asexual parasite life stages
        - units_ul_in:
            True if raw_parasitemia is in parasites/uL
            False if raw_parasitemia is in fractional percentage
        - units_ul_out (optional):
            True to return parasitemia in parasitemia/uL
            False to return parasitemia in fractional percentage (default)
        """

        if raw_parasitemia < 0:
            raw_parasitemia = 0

        if units_ul_in:
            raw_parasitemia /= RBCS_P_UL

        # Compute counts based on parasitemia and rbcs
        parasites = raw_parasitemia * rbcs
        healthy = rbcs - parasites
        counts = [healthy, parasites]

        compensated_parasitemia, bound = self._get_res_from_counts(
            counts, units_ul_out=units_ul_out
        )

        return compensated_parasitemia, bound
