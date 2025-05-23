# FOCCUS
Softwares developped and used at Shom in the context of the Horizon Europe FOCCUS project (https://foccus-project.eu/).

# Installation

```
python setup.py install --user
```

# Convert GLO4ens native model outputs to CMEMS-type data structure

- Retrieve GLO4ens model output on the native NEMO grid, with one file per ensemble member, per day and per 3D variables. 2D variables, including SSH, were all combined in one file, the script ```extract_zos.sh``` is used to extract it (following the data structure below). 

- From original ```.tar``` files, the data was first organized follawing this convention:```/YOUR_LOCAL_DIR/YYYY/MM/memMMM/*.nc```

- ```convert_as_cmems.ipynb``` provide an use example to convert GLO4ens from native NEMO to CMEMS-type data format.
