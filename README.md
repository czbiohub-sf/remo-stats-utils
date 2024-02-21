# Statistics utils
**Live version: v0.0.1 (as of 2024-02-21)**

For analysis of remoscope results.

This repo includes the data files for deskewing and parasitemia compensation. See "Data files" below for file naming/organization convention.

## Package installation
From top-level folder, run
```console
python3 -m pip install -e .
```

## Data files
It's important to ensure files are named and organized appropriately so they can be dynamically loaded. Dynamic loading is used to match the data with the YOGO model being run in [ulc-malaria-scope](https://github.com/czbiohub-sf/ulc-malaria-scope).

Let the model ID include the model name and number, separated by dashes:
* Group files by model ID in the subfolder ```data_files/<model ID>```
* Name data files ```<model ID><suffix>```, where the suffixes are defined in ```stats_utils/constants.py```

For example, for model ID ```frightful-wendigo-1981```:
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

