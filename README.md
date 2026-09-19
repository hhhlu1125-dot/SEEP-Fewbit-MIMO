# SEEP-Fewbit-MIMO

MATLAB implementation of the **SEEP** and **ApproSEEP** detectors for MIMO systems with few-bit ADCs.

This repository is provided to facilitate reproducibility of the main BER results reported in:

> **SEEP: Sample Reconstruction-based Expectation Propagation Detection for MIMO Systems with Few-bit ADCs**

## Repository Structure

The MATLAB source code is provided in the `Matlab code/` directory.

```text
Matlab code/
├── reproduce_BER_results.m
├── detectors/
├── baselines/
├── utils/
└── analysis/
```

## Quick Start

Run

```matlab
reproduce_BER_results
```

to perform the BER simulations.

The detector type, antenna configuration, modulation order, ADC resolution, SNR range, and other simulation parameters can be specified in `reproduce_BER_results.m`.

## Main Algorithms

- `Matlab code/detectors/SEEP_detection.m`
- `Matlab code/detectors/ApproSEEP_detection.m`
