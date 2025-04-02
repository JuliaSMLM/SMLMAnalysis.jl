# SMLMAnalysis

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSMLM.github.io/SMLMAnalysis.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSMLM.github.io/SMLMAnalysis.jl/dev/)
[![Build Status](https://github.com/JuliaSMLM/SMLMAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSMLM/SMLMAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSMLM/SMLMAnalysis.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSMLM/SMLMAnalysis.jl)

## Overview

SMLMAnalysis is a high-level umbrella package that organizes and provides a unified interface to the JuliaSMLM ecosystem for Single Molecule Localization Microscopy (SMLM) data analysis. This package streamlines common SMLM workflows by connecting specialized packages into integrated pipelines.

## Package Architecture

SMLMAnalysis integrates several specialized SMLM packages to provide complete analysis workflows:

```
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│  MicroscopePSFs │    │    SMLMSim      │    │PSFConvolutions  │
│  (PSF Models)   │    │  (Simulations)  │    │ (PSF Operations)│
└────────┬────────┘    └────────┬────────┘    └────────┬────────┘
         │                      │                      │
         └──────────────────────┼──────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────┐
│             SMLMDeepFit                 │
│     (Deep Learning Localization)        │
└────────────────────┬───────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────┐
│          SMLMDriftCorrection            │
│      (Spatiotemporal Stabilization)     │
└────────────────────┬───────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────┐
│              SMLMBaGoL                  │
│      (Bayesian Grouping of Localizations)│
└───────────┬────────────────────┬────────┘
            │                    │
            ▼                    ▼
┌────────────────┐      ┌────────────────┐
│   SMLMVis      │      │  SMLMMetrics   │
│(Visualization) │      │   (Analysis)   │
└────────────────┘      └────────────────┘


┌─────────────┐     ┌─────────────┐
│    Boxer    │     │  GaussMLE   │
│ (Detection) │────▶│ (Fitting)   │───┐
└─────────────┘     └─────────────┘   │
                                       │
                                       ▼
                    ┌─────────────────────────────┐
                    │    SMLMDriftCorrection     │
                    └───────────────┬─────────────┘
                                    │
                                    ▼
                    ┌─────────────────────────────┐
                    │        SMLMBaGoL           │
                    └──────────┬─────────┬────────┘
                               │         │
                               ▼         ▼
                    ┌────────────────┐  ┌────────────────┐
                    │   SMLMVis      │  │  SMLMMetrics   │
                    │(Visualization) │  │   (Analysis)   │
                    └────────────────┘  └────────────────┘
```

## Key Components

### Workflow 1: Deep Learning Pipeline

- **Input Preparation**: 
  - `MicroscopePSFs.jl` - Provides accurate PSF models for training and simulation
  - `SMLMSim.jl` - Generates realistic SMLM data with ground truth positions
  - `PSFConvolutions.jl` - Provides efficient convolution operations for PSF modeling

- **Localization**:
  - `SMLMDeepFit.jl` - Performs deep learning-based localization of single molecules

- **Post-Processing**:
  - `SMLMDriftCorrection.jl` - Corrects for sample drift and stage movement
  - `SMLMBaGoL.jl` - Performs Bayesian grouping of localizations for super-resolution reconstruction

- **Output & Analysis**:
  - `SMLMVis.jl` - Visualizes localization data and super-resolution reconstructions
  - `SMLMMetrics.jl` - Provides quantitative metrics for evaluating localization and reconstruction quality

### Workflow 2: Traditional Analysis Pipeline

- **Detection & Fitting**:
  - `Boxer.jl` - Detects candidate single-molecule spots in raw images
  - `GaussMLE.jl` - Performs Maximum Likelihood Estimation fitting with Gaussian PSF model

- **Post-Processing**:
  - `SMLMDriftCorrection.jl` - Corrects for sample drift and stage movement
  - `SMLMBaGoL.jl` - Performs Bayesian grouping of localizations for super-resolution reconstruction

- **Output & Analysis**:
  - `SMLMVis.jl` - Visualizes localization data and super-resolution reconstructions
  - `SMLMMetrics.jl` - Provides quantitative metrics for evaluating localization and reconstruction quality

## Installation

```julia
using Pkg
Pkg.add("SMLMAnalysis")
```

This will automatically install all the required components of the JuliaSMLM ecosystem.

## Basic Usage

SMLMAnalysis provides a simplified interface to common SMLM workflows. The package offers:

- Configuration interfaces for setting up analysis parameters
- Data loading and conversion utilities
- Pipeline management for both deep learning and traditional analysis workflows
- Unified data structures for storing and passing results between components
- Automated handling of file and data formats

## Contributing

Contributions to SMLMAnalysis are welcome! Please feel free to:

1. Report bugs and request features via GitHub issues
2. Submit pull requests with bug fixes or enhancements
3. Improve documentation or examples

## License

This project is licensed under the MIT License - see the LICENSE file for details.