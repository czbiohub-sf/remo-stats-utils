"""
YOGO class compensation using y = mx + b fit

Based on clinical Uganda data
"""

import numpy as np
import numpy.typing as npt
import pandas as pd

from typing import List, Tuple, Union
from pathlib import Path

from stats_utils.constants import (
    DATA_DIR,
    CULTURED_COMPENSATION_SUFFIX1,
    CLINICAL_COMPENSATION_SUFFIX1,
    NO_HEATMAPS_SUFFIX2,
    W_HEATMAPS_SUFFIX2,
    CONFIDENCE_THRESHOLD,
    PARASITES_P_UL_PER_PERCENT,
)
from stats_utils.corrector import CountCorrector


class CountCompensator(CountCorrector):
    def __init__(self, model_name: str, clinical=True, heatmaps=False):
        """
        Initialize count compensator

        Input(s)
        - model_name: 
            Include name and number (eg. "frightful-wendigo-1931")
        - clinical: 
            True for clinical Uganda data
            False for cultured lab data
        - heatmaps: 
            True for heatmap nuked data
            False otherwise
        """
        
        # Generate directory for compensation metrics csv
        if clinical:
            suffix1 = CLINICAL_COMPENSATION_SUFFIX1
        else:
            suffix1 = CULTURED_COMPENSATION_SUFFIX1
        if heatmaps:
            suffix2 = W_HEATMAPS_SUFFIX2
        else:
            suffix2 = NO_HEATMAPS_SUFFIX2
        compensation_csv_dir = str(DATA_DIR / model_name / (model_name + suffix1 + suffix2))

        # Check that compensation metrics csv exists
        if not Path(compensation_csv_dir).is_file():
            raise FileNotFoundError(
                f"Could not find compensation metrics for {model_name} ({compensation_csv_dir})"
            )

        m, b, cov_m, cov_b = self.get_fit_metrics(compensation_csv_dir)
        inv_cmatrix = self.get_matrix(m, b)
        inv_cmatrix_std = self.get_matrix_std(cov_m, cov_b)

        rbc_ids = [0, 1]
        parasite_ids = [1]

        super(CountCompensator, self).__init__(
            inv_cmatrix,
            inv_cmatrix_std,
            rbc_ids,
            parasite_ids,
        )

    def get_fit_metrics(self, compensation_csv_dir: str) -> Tuple[float, float, float, float]:
        """
        Extract fit metrics from csv

        Returns fit metrics as (m, b, cov_m, cov_b)
        """

        df = pd.read_csv(compensation_csv_dir, dtype=np.float64)
        row = df.loc[df['conf_val'] == CONFIDENCE_THRESHOLD]

        # Adjust b for parasitemia % instead of parasites per uL
        fit_b = row["fit_b"] / PARASITES_P_UL_PER_PERCENT
        cov_b = row["cov_b"] / PARASITES_P_UL_PER_PERCENT

        return row["fit_m"], fit_b, row["cov_m"], cov_b

    def get_matrix(self, m: float, b: float) -> npt.NDArray:
        """
        Return inverse 2x2 confusion matrix

        See supplementary materials for derivation
        """
        M12 = -b
        M22 = (1 / m) - b

        M11 = 1 - M12
        M21 = 1 - M22

        return np.array([[M11, M12], [M21, M22]])

    def get_matrix_std(self, cov_m: float, cov_b: float) -> npt.NDArray:
        """
        Return standard deviations inverse 2x2 confusion matrix

        See supplementary materials for derivation
        """
        M11 = cov_b
        M12 = cov_b

        M21 = np.sqrt(cov_b**2 + cov_m**2)
        M22 = np.sqrt(cov_b**2 + cov_m**2)

        return np.array([[M11, M12], [M21, M22]])

    def get_95_bound_and_compensation_from_parasitemia(
        self,
        raw_parasitemia: float,
        rbcs: Union[int, float],
        units_ul_in: bool,
        units_ul_out: bool = False,
    ) -> Tuple[float, List[float]]:
        """
        Return compensated parasitemia and 95% confidence bound based on raw parasitemia and total rbcs

        See remoscope manuscript for full derivation

        Input(s)
        - raw_parasitemia:
            Uncorrected parasitemia estimate
        - rbcs:
            Total count of rbcs
        - units_ul_in:
            True if raw_parasitemia is in parasites/uL 
            False if raw_parasitemia is in %
        - units_ul_out:
            True to return parasitemia in parasitemia/uL
            False to return parasitemia in %
        """
        if units_ul_in:
            raw_parasitemia /= PARASITES_P_UL_PER_PERCENT
        parasitemia_fraction = raw_parasitemia / 100.0

        # Compute counts based on parasitemia and rbcs
        parasites = parasitemia_fraction * rbcs
        healthy = rbcs - parasites
        counts = [healthy, parasites]

        compensated_parasitemia, bound = self.get_res_from_counts(
            counts, units_ul_out=units_ul_out
        )

        return compensated_parasitemia, self.get_95_confidence_bound(
            compensated_parasitemia, bound
        )
