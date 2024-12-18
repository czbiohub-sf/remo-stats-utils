# Statistics utils
This repo contains the statistical analysis tools used to compute uncertainty estimates for Remoscope data. The calculations are used both on the Remoscope instrument as well as Supplementary Note 2 in the Remoscope preprint for analytical expressions of the calculations: 

![image](https://github.com/user-attachments/assets/0cfd37a6-c543-4873-9a64-bd2e5f8506b7)

There are also additional tools, such as matrix deskewing analysis, which we explored but did not get implemented in the main project yet. 

## Local package installation
```console
python3 -m pip install .
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

