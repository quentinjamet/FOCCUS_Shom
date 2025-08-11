# Description
Softwares developped and used at Shom in the context of the Horizon Europe FOCCUS project (https://foccus-project.eu/) to convert model output from native NEMO grid to a CMEMS-like format, including a capabilities to deal with ensemble simulation model outputs.


# Installation

```
python setup.py install --user
```

# Test case: the global $frac{1}{4}^{\circ}$ ensemble simulation GLO4ens

## Step-by-step procedure

- Retrieve GLO4ens model output on natice NEMO grid from Mercator Ocean International.

- Extract SSH field from 2D surface variables with ```extract_zos.sh```.

- Define various parameters in ```params.json``` file, including dates, ensemble and variable name information (see ```convert_to_cmems_type``` description for further informations).

- ```convert_GLO4ens.ipynb``` provides an use example of such convertion.
