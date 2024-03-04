from typing import Dict, List

from pathlib import Path

curr_dir = Path(__file__).parent.resolve()

# Path to data files
DATA_DIR = curr_dir / "data_files"

# Suffixes for data files
CMATRIX_MEAN_SUFFIX = "-cmatrix-mean.npy"
INV_CMATRIX_STD_SUFFIX = "-inv-cmatrix-std.npy"

CLINICAL_COMPENSATION_SUFFIX1 = "-clinical-compensation"
CULTURED_COMPENSATION_SUFFIX1 = "-cultured-compensation"

W_HEATMAPS_SUFFIX2 = "-with-heatmaps.csv"
NO_HEATMAPS_SUFFIX2 = "-no-heatmaps.csv"

# Confidence threshold
CONFIDENCE_THRESHOLD = 0.90

# YOGO class IDs
YOGO_CLASS_ORDERING = [
    "healthy",
    "ring",
    "trophozoite",
    "schizont",
    "gametocyte",
    "wbc",
    "misc",
]
YOGO_CLASS_IDX_MAP: Dict[str, int] = {
    k: idx for idx, k in enumerate(YOGO_CLASS_ORDERING)
}
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
PARASITES_P_UL_PER_PERCENT = 5e6  # parasites/uL = parasitemia % x 5E6
