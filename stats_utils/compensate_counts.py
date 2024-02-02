"""
YOGO class compensation using y = mx + b fit

Based on clinical Uganda data
"""

import numpy as np
import numpy.typing as npt

from typing import Tuple

from stats_utils.constants import (
    YOGO_COMPENSATION_CSV_DIR,
    CONFIDENCE_THRESHOLD,
)


class CountCompensator:
    def __init__(self):
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

        df = pd.read_csv(YOGO_COMPENSATION_CSV_DIR, dtype=np.float64)
        row = df.loc[df['conf_val'] == CONFIDENCE_THRESHOLD]

        return row['fit_m'], row['fit_b'], row['cov_m'], row['cov_b']

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