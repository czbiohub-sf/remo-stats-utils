"""
YOGO class compensation using y = mx + b fit

Based on clinical Uganda data
"""

import numpy as np
import numpy.typing as npt

from typing import Tuple, Union

from stats_utils.constants import (
    DATA_DIR,
    CULTURED_COMPENSATION_SUFFIX1,
    CLINICAL_COMPENSATION_SUFFIX1,
    NO_HEATMAPS_SUFFIX2,
    W_HEATMAPS_SUFFIX2,
    CONFIDENCE_THRESHOLD,
    PARASITES_P_UL_PER_PERCENT
)
from stats_utils.correct_counts import CountCorrector


class CountCompensator(CountCorrector):
    def __init__(self, model_name: str, clinical=True, heatmaps=False):
        
        # Generate directory for compensation metrics csv
        if clinical:
            suffix1 = CLINICAL_COMPENSATION_SUFFIX1
        else:
            suffix1 = CULTURED_COMPENSATION_SUFFIX1
        if heatmaps:
            suffix2 = W_HEATMAPS_SUFFIX2
        else:
            suffix2 = NO_HEATMAPS_SUFFIX2
        self.compensation_csv_dir = str(DATA_DIR / model_name / (model_name + suffix1 + suffix2))

        # Check that compensation metrics csv exists
        if not self.compensation_csv_dir.is_dir():
            raise FileNotFoundError(
                f"Could not find {model_name} compensation metrics file {self.compensation_csv_dir}"
            )

        m, b, cov_m, cov_b = self.get_fit_metrics()
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

    def get_fit_metrics(self) -> Tuple[float, float, float, float]:
        """
        Extract fit metrics from csv

        Returns fit metrics as (m, b, cov_m, cov_b)
        """        

        df = pd.read_csv(self.compensation_csv_dir, dtype=np.float64)
        row = df.loc[df['conf_val'] == CONFIDENCE_THRESHOLD]

        # Adjust b for parasitemia % instead of parasites per uL
        fit_b = row['fit_b'] / PARASITES_P_UL_PER_PERCENT
        cov_b = row['cov_b'] / PARASITES_P_UL_PER_PERCENT

        return row['fit_m'], fit_b, row['cov_m'], cov_b

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
        
        M21 = np.sqrt(cov_b ** 2 + cov_m ** 2)
        M22 = np.sqrt(cov_b ** 2 + cov_m ** 2)

        return np.array([[M11, M12], [M21, M22]])

    def get_95_bound_and_compensation_from_parasitemia(
        self, raw_parasitemia: float, rbcs: Union[int, float], units_ul_in: bool, units_ul_out: bool=False
    ) -> Tuple[float, List[float]]:
        """
        Return compensated parasitemia and 95% confidence bound based on raw parasitemia and total rbcs

        See remoscope manuscript for full derivation

        95% confidence interval can be defined as
            lower_bound = max(0, parasitemia - bound)
            upper_bound = min(1, parasitemia + bound)
        """
        
        if units_ul_in:
            raw_parasitemia /= PARASITES_P_UL_PER_PERCENT
        parasitemia_fraction = raw_parasitemia / 100.0
        
        # Compute counts based on parasitemia and rbcs
        parasites = parasitemia_fraction * rbcs
        healthy = rbcs - parasites
        counts = [healthy, parasites]

        compensated_parasitemia, bound = self.get_res_from_counts(counts, units_ul_out=units_ul_out)

        return compensated_parasitemia, self.get_95_confidence_bound(compensated_parasitemia, bound)
