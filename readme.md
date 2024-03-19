# openmm-nf

![](https://img.shields.io/badge/current_version-0.1.34-blue)
![](https://github.com/stracquadaniolab/openmm-nf/workflows/build/badge.svg)
## Overview
A simple workflow to evaluate protein stability with respect to the wildtype

## Configuration

- param1: this is the parameter description (default: "hello")
- param2: this is the parameter description (default: "world")
- ...
- paramN: this is the parameter description (default: "flow")

## Running the workflow

### Install or update the workflow

```bash
nextflow pull stracquadaniolab/openmm-nf
```

### Run the analysis

```bash
nextflow run stracquadaniolab/openmm-nf
```

## Results

- `results/analysis.txt`: file with the analysis results.
- `results/tuning.txt`: file with the parameter tuning.
- ...

## Authors

- Josh David Littlefair
