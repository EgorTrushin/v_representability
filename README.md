### Violations of the $v$-representability condition underlying Kohn-Sham density-functional theory

Script [gen_inputs.py](https://github.com/EgorTrushin/vrepr-violation/blob/main/gen_inputs.py) creates all input files for Molpro, besides inputs for KS inversion using FCI target densities and other calculations which supplement FCI-based results. The script uses config files from the **configs** directory.

<details>
<summary>Example for generating inputs</summary>

```
./gen_inputs.py configs/config_EXX_AllSymCases.yaml
./gen_inputs.py configs/config_EXX_SpinSym_plot.yaml
./gen_inputs.py configs/config_CAHF.yaml
./gen_inputs.py configs/config_ROHF.yaml
./gen_inputs.py configs/config_KSINV_ROHF.yaml
./gen_inputs.py configs/config_KSINV_CAHF.yaml
```

</details>

Notebook [process_raw_outputs.ipynb](https://github.com/EgorTrushin/vrepr-violation/blob/main/process_raw_outputs.ipynb) processes Molpro output files. Processed data are stored in Python-friendly formats. Data inspection is also performed.

Notebook [vrepr-violation.ipynb](https://github.com/EgorTrushin/vrepr-violation/blob/main/vrepr-violation.ipynb) generates tables and figures, i.e., provides analytics. [vrepr-violation.html](https://github.com/EgorTrushin/vrepr-violation/blob/main/vrepr-violation.html) contains only output of [vrepr-violation.ipynb](https://github.com/EgorTrushin/vrepr-violation/blob/main/vrepr-violation.ipynb), i.e., cells with codes are absent.
