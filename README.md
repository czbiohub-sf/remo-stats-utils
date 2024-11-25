# Statistics utils

Generates parasitemia estimates and corresponding 95% confidence bounds for Remoscope. Intended for use with [YOGO](https://github.com/czbiohub-sf/yogo) classifications.

This repo includes the data files for deskewing and parasitemia compensation . See "Data files" below for file naming/organization convention.

## Local package installation
From top-level folder, run
```console
python3 -m pip install -e .
```

As an alternative to local installation, point to this repo and the desired release version in the `install_requires` argument of your codebase's `setup.py`.

## Usage
Use `CountCompensator` to correct parasitemia estimates and compute 95% confidence bounds according to a linear fit `y = mx + b`. Example instantiation (see documentation in `compensator.py` for input argument descriptions):
```console
# Instantiate compensator using y = mx + b fit with corresponding error
# Fit is based on 0.90 confidence thresholded frightful-wendigo model classification vs clinical PCR
corrector = CountCompensator("frightful-wendigo-1931", 0.90)

# Skip compensation, using computation of parasitemia and error from raw data only
# Ignores model and confidence threshold arguments
corrector = CountCompensator("frightful-wendigo-1931", 0.90, skip=True)
```

Use `CountDeskewer` to correct for skew in class counts according to the confusion matrix. Example instantiation:
```console
# Instantiate deskewer based on frightful-wendigo model confusion matrix
corrector = CountDeskewer("frightful-wendigo-1931")
```

Based on the input arguments, the fit values and confusion matrices are extracted from the appropriate .csv in `data_files`. See "Data files" for more details.

To compute parasitemia and corresponding 95% confidence bounds:
```console
# Example YOGO output
class_counts = np.array([
    100000, # healthy
    60, # ring
    40, # troph
    20, # schizont
    10, # gametocyte
    150, # WBC
    200, # misc
])

# Get parasitemia only
parasitemia = corrector.calc_parasitemia(class_counts)

# Get parasitemia and 95% confidence bounds as parasites/uL
parasitemia, conf_bounds = corrector.get_res_from_counts(class_counts, units_ul_out=True)

# Get parasitemia and 95% confidence bounds as percentage
parasitemia, conf_bounds = corrector.get_res_from_counts(class_counts, units_ul_out=False)
```

## Data files
`remo-stats-utils` requires the data files to be organized in a particular schema for dynamic loading. Dynamic loading is used to match the data with the YOGO model being run in [ulc-malaria-scope](https://github.com/czbiohub-sf/ulc-malaria-scope).

Let the model ID include the model name and number, separated by dashes:
* Group files by model ID in the subfolder ```data_files/<model ID>```
* Name data files ```<model ID><suffix>```, where the suffixes are defined in ```stats_utils/constants.py```
     * The suffix describes whether the data is generated from clinical vs cultured data and whether heatmap nuking was used in the data processing

For example, for the model ```frightful-wendigo-1981```, one may have the following file structure:
```
data_files/
├── frightful-wendigo-1931/
│   ├── frightful-wendigo-1931-cmatrix-mean.npy
│   ├── frightful-wendigo-1931-inv-cmatrix-std.npy
│   ├── frightful-wendigo-1931-cultured-compensation-no-heatmaps.csv
│   ├── frightful-wendigo-1931-cultured-compensation-with-heatmaps.csv
│   ├── frightful-wendigo-1931-clinical-compensation-no-heatmaps.csv
│   ├── frightful-wendigo-1931-clinical-compensation-with-heatmaps.csv
├── other-model-0000/
│   ... 
```

