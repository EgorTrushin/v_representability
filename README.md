### Violations of the $v$-representability condition underlying Kohn-Sham density-functional theory

This repository contains the analysis of data and the preparation of tables and figures for the following publication:

- E. Trushin, J. Erhard and A. Görling. Violations of the v-representability condition underlying Kohn-Sham density-functional theory - Phys. Rev. A 110 L020802 (2024) https://doi.org/10.1103/PhysRevA.110.L020802
    - Click to see corresponding [Jupyter Notebook](https://github.com/EgorTrushin/supplementary-materials/vrepr_violation/vrepr-violation.ipynb)
    - Zenodo repository: [http://doi.org/10.5281/zenodo.8226492](http://doi.org/10.5281/zenodo.12772192)

---

Script [gen_inputs.py](https://github.com/EgorTrushin/v_representability/gen_inputs.py) creates all input files for Molpro, besides inputs for KS inversion using FCI target densities and other calculations which supplement FCI-based results. The script uses config files from the **configs** directory.

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

Notebook [process_raw_outputs.ipynb](https://github.com/EgorTrushin/v_representability/blob/main/process_raw_outputs.ipynb) processes Molpro output files. Processed data are stored in Python-friendly formats. Data inspection is also performed.

Notebook [vrepr-violation.ipynb](https://github.com/EgorTrushin/v_representability/vrepr-violation.ipynb) generates tables and figures, i.e., provides analytics. [vrepr-violation.html](https://github.com/EgorTrushin/v_representability/vrepr-violation.html) contains only output of [vrepr-violation.ipynb](https://github.com/EgorTrushin/v_representability/vrepr-violation.ipynb), i.e., cells with codes are absent.
