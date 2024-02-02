from yogo.data import YOGO_CLASS_ORDERING
from pathlib import Path

from typing import Dict, List


# Paths to data files
curr_dir = Path(__file__).parent.resolve()
YOGO_CMATRIX_MEAN_DIR = str(curr_dir / "data_files" / "frightful-wendigo-1931-cmatrix-mean.npy")
YOGO_INV_CMATRIX_STD_DIR = str(curr_dir / "data_files" / "frightful-wendigo-1931-inverse-cmatrix-std.npy")
YOGO_COMPENSATION_CSV_DIR= str(curr_dir / "data_files" / "frightful-wendigo-1931-with-heatmaps_titration_metrics.csv")

# Confidence threshold
CONFIDENCE_THRESHOLD = 0.90

# YOGO class IDs
YOGO_CLASS_IDX_MAP: Dict[str, int] = {k: idx for idx, k in enumerate(YOGO_CLASS_ORDERING)}
RBC_CLASS_IDS: List[int] = [
    YOGO_CLASS_IDX_MAP["healthy"],
    YOGO_CLASS_IDX_MAP["ring"],
    YOGO_CLASS_IDX_MAP["trophozoite"],
    YOGO_CLASS_IDX_MAP["schizont"],
]
ASEXUAL_PARASITE_CLASS_IDS: List[int] = [
    YOGO_CLASS_IDX_MAP["ring"],
    YOGO_CLASS_IDX_MAP["trophozoite"],
    YOGO_CLASS_IDX_MAP["schizont"],
]
PARASITE_CLASS_IDS: List[int] = [
    YOGO_CLASS_IDX_MAP["ring"],
    YOGO_CLASS_IDX_MAP["trophozoite"],
    YOGO_CLASS_IDX_MAP["schizont"],
    YOGO_CLASS_IDX_MAP["gametocyte"],
]

# Parasitemia unit conversion
PERCENT_2_PARASITES_PER_UL = 5E6 # parasites/uL = parasitemia % x 5E6