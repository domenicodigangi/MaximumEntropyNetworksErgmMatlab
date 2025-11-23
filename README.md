# Maximum Entropy Network Reconstruction for ERGM

[![MATLAB](https://img.shields.io/badge/MATLAB-R2014a+-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Research](https://img.shields.io/badge/Research-Published-blue.svg)](https://www.sciencedirect.com/science/article/pii/S0165188918301787)
[![Documentation](https://github.com/domenicodigangi/MaximumEntropyNetworksErgmMatlab/actions/workflows/documentation.yml/badge.svg)](https://github.com/domenicodigangi/MaximumEntropyNetworksErgmMatlab/actions/workflows/documentation.yml)
[![MATLAB Validation](https://github.com/domenicodigangi/MaximumEntropyNetworksErgmMatlab/actions/workflows/matlab-validation.yml/badge.svg)](https://github.com/domenicodigangi/MaximumEntropyNetworksErgmMatlab/actions/workflows/matlab-validation.yml)

> **Note**: This repository contains code developed during my Master's thesis research. It is maintained as an archive of published academic work and is no longer actively developed.

## Overview

This MATLAB package implements **Maximum Entropy Network Reconstruction** methods for **Exponential Random Graph Models (ERGM)**. It provides tools for estimating and sampling network ensembles when only partial information about the network structure is available.

The primary application is **assessing systemic risk in financial networks** through fire sales spillover analysis, particularly when complete network data is unavailable.

## Published Research

This code was developed as part of the following publication:

**Di Gangi, D., Lillo, F., & Pirino, D. (2018).** *Assessing systemic risk due to fire sales spillover through maximum entropy network reconstruction.* Journal of Economic Dynamics and Control, 94, 301-319.

📄 [Read the paper](https://www.sciencedirect.com/science/article/pii/S0165188918301787)

## Key Features

- **Multiple Network Models**: Implements various maximum entropy network ensembles
  - Binary Configuration Model (BIPCM)
  - Weighted Configuration Model (BIPWCM)
  - Enhanced Configuration Model (BIPECM)
  - CAPM-based models (MECAPM, EMECAPM)
  - Density-corrected variants

- **Systemic Risk Analysis**: Implementation of the Vulnerable Banks framework for fire sales spillover assessment

- **Flexible Estimation**: Analytical and numerical optimization methods for parameter estimation

- **Ensemble Sampling**: Generate synthetic networks from estimated ensembles

## Mathematical Background

Maximum entropy network reconstruction allows us to:

1. **Estimate unknown network structures** given partial information (e.g., only node strengths/degrees)
2. **Quantify uncertainty** in network-based risk measures
3. **Avoid biases** from arbitrary network assumptions

The method finds the least-biased probability distribution over all possible networks consistent with known constraints, based on the **maximum entropy principle**.

## Installation

### Requirements

- MATLAB R2014a or later
- Optimization Toolbox (required for most models)

### Setup

1. Clone this repository:
   ```bash
   git clone https://github.com/domenicodigangi/MaximumEntropyNetworksErgmMatlab.git
   cd MaximumEntropyNetworksErgmMatlab
   ```

2. Add to MATLAB path:
   ```matlab
   addpath(genpath('/path/to/MaximumEntropyNetworksErgmMatlab'))
   ```

## Usage

### Basic Example

```matlab
% Example: Estimate a bipartite weighted configuration model

% Input data: strength sequences for investors and assets
investor_strengths = [100; 150; 200; 120];  % Row sums
asset_strengths = [170; 200; 200];          % Column sums
in_data = {investor_strengths, asset_strengths};

% Estimate the model
model = Max_Entr_Nets('BIPWCM', in_data);

% Sample networks from the ensemble
n_samples = 100;
network_sample = model.sample(n_samples);

% Get expected adjacency matrix
expected_network = model.exp_matrix();
```

### Available Models

List all available models:
```matlab
Max_Entr_Nets('LIST')
```

List bipartite models only:
```matlab
Max_Entr_Nets('LIST-BIP')
```

### Model Descriptions

| Model | Description | Required Data |
|-------|-------------|---------------|
| **BIPCM** | Binary Configuration Model | Degree sequences |
| **BIPWCM** | Weighted Configuration Model | Strength sequences |
| **BIPECM** | Enhanced Configuration Model | Strength + degree sequences |
| **MECAPM** | Maximum Entropy CAPM | Strength sequences |
| **EMECAPM** | Enhanced MECAPM | Strength + degree sequences |
| **DCBIPWCM** | Density-Corrected WCM | Strength sequences + density |

See individual model files in `models/` for detailed descriptions.

### Systemic Risk Analysis

```matlab
% Compute systemic risk measures using the Vulnerable Banks framework

% Prepare data
X = [...];  % Adjacency matrix (investors × assets)
equity = [...];  % Investor equity values
shock = [0.01, 0.01, ...];  % Asset price shocks (e.g., 1% depreciation)

% Compute with known network
[AV, SYS, VUL] = Vulnerable_Banks('REAL', X, equity, shock);

% Estimate from partial information
in_data = {investor_caps, asset_caps};  % Only capitalization known
[AV_est, SYS_est, VUL_est] = Vulnerable_Banks('ESTIMATE', in_data, equity, shock);

% Output:
%   AV  - Aggregate Vulnerability (system-wide risk)
%   SYS - Systemicness (contribution to system risk)
%   VUL - Vulnerability (exposure to fire sales)
```

## Project Structure

```
.
├── Max_Entr_Nets.m          # Main function for network estimation
├── Vulnerable_Banks.m        # Systemic risk analysis
├── models/                   # Network model implementations
│   ├── BIPCM.m              # Binary Configuration Model
│   ├── BIPWCM.m             # Weighted Configuration Model
│   ├── BIPECM.m             # Enhanced Configuration Model
│   ├── MECAPM.m             # Maximum Entropy CAPM
│   ├── EMECAPM.m            # Enhanced MECAPM
│   ├── DCBIPWCM.m           # Density-Corrected models
│   └── DCMECAPM.m
├── useful_functions/         # Utility functions
│   ├── indata_from_matrix_Nets.m
│   └── list_to_mat_uni_und.m
└── readme.pdf               # Original technical documentation
```

## Technical Details

### Precision Control

Control the optimization precision:
```matlab
precision = 1e-3;  % Maximum relative error tolerance
model = Max_Entr_Nets('BIPWCM', in_data, precision);
```

### Custom Parameters

If parameters are pre-computed:
```matlab
known_parameters = [...];
model = Max_Entr_Nets('BIPWCM', in_data, precision, known_parameters);
```

## References

If you use this code in your research, please cite:

```bibtex
@article{DiGangi2018,
  title={Assessing systemic risk due to fire sales spillover through maximum entropy network reconstruction},
  author={Di Gangi, Domenico and Lillo, Fabrizio and Pirino, Davide},
  journal={Journal of Economic Dynamics and Control},
  volume={94},
  pages={301--319},
  year={2018},
  publisher={Elsevier},
  doi={10.1016/j.jedc.2018.07.001}
}
```

### Related Work

The Vulnerable Banks framework is based on:

- **Greenwood, R., Landier, A., & Thesmar, D. (2015).** *Vulnerable banks.* Journal of Financial Economics, 115(3), 471-485.

- **Duarte, F., & Eisenbach, T. M. (2015).** *Fire-sale spillovers and systemic risk.* FRB of New York Staff Report No. 645.

## Author

**Domenico Di Gangi**
Email: domenico.digangi@sns.it

Developed as part of Master's thesis research at Scuola Normale Superiore.

## Contributing

This is archived research code, but contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on:
- Bug fixes and documentation improvements
- Compatibility updates
- Performance optimizations
- How to report issues

## Citation

If you use this code in your research, please cite the original paper. A `CITATION.cff` file is provided for easy citation management with tools like Zotero, Mendeley, and GitHub.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Archive Notice

This repository is maintained as an **archive** of published research code. It represents the implementation used in the published paper and is provided for:

- **Reproducibility** of research results
- **Reference** for researchers in network reconstruction and systemic risk
- **Educational purposes**

For production use or extensions, consider implementing modern versions with updated dependencies and best practices.

## Acknowledgments

This research was conducted with the support of:
- Scuola Normale Superiore di Pisa
- IMT School for Advanced Studies Lucca

Special thanks to Prof. Fabrizio Lillo and Dr. Davide Pirino for their supervision and collaboration.
