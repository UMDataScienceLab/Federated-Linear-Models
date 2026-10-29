# Federated-Linear-Models
Federated-Linear-Models

Federated hierarchical linear regression across devices that keep their data local.

CM_data.R: used to analyze the NASA data

variable_selection.stan: the federated EP models. Please install the Rstan package first.

The NASA dataset: https://ti.arc.nasa.gov/tech/dash/groups/pcoe/

The UCI dataset: https://archive.ics.uci.edu/ml/index.php

The windfarm data is confidential.

## Citation

[Paper](https://doi.org/10.1080/24725854.2022.2157912) (IISE Transactions, 2024). If you use this code, please cite:

```bibtex
@article{yue2024federated,
  title={Federated Data Analytics: A Study on Linear Models},
  author={Yue, Xubo and Al Kontar, Raed and Estrada G{\'o}mez, Ana Mar{\'i}a},
  journal={IISE Transactions},
  volume={56},
  number={1},
  pages={16--28},
  year={2024},
  doi={10.1080/24725854.2022.2157912}
}
```
