# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

SMLMAnalysis.jl is the high-level integration package for the JuliaSMLM ecosystem. It serves as a wrapper/orchestrator that brings together all the individual SMLM (Single Molecule Localization Microscopy) analysis packages into coherent workflows. While currently a stub, this package is intended to be the primary entry point for users to access the complete SMLM analysis pipeline.

## JuliaSMLM Ecosystem

The JuliaSMLM organization provides a comprehensive suite of packages for SMLM data analysis:

### Core Data Infrastructure
- **SMLMData.jl**: Core data types and utilities for SMLM coordinate data. Defines fundamental types like emitters, cameras, and SMLD containers. All other packages build on these types.

### Detection & Fitting
- **SMLMBoxer.jl**: Fast particle detection in multidimensional image stacks using Difference of Gaussians
- **GaussMLE.jl**: Maximum likelihood Gaussian blob fitting with CPU/GPU support
- **SMLMDeepFit.jl**: Deep learning-based high-density emitter fitting using U-Net architectures

### PSF Modeling
- **MicroscopePSFs.jl**: Microscope point spread function calculations including scalar and vector models

### Post-Processing
- **SMLMFrameConnection.jl**: Connects repeated localizations across frames into single higher-precision localizations

### Analysis & Metrics
- **SMLMMetrics.jl**: Performance metrics for SMLM including Jaccard Index, RMSE, and efficiency metrics

### Simulation
- **SMLMSim.jl**: Comprehensive simulation of SMLM datasets including photophysics, diffusion, and camera effects

### Visualization
- **SMLMVis.jl**: Visualization tools for SMLM data and analysis results

### Infrastructure
- **ModelContextProtocol.jl**: MCP server implementation for tool integration

## Development Commands

### Testing
```bash
# Run full test suite
julia --project=. -e 'using Pkg; Pkg.test()'

# Run tests directly
julia --project=. test/runtests.jl
```

### Documentation
```bash
# Build documentation
julia --project=docs docs/make.jl

# Serve documentation locally (if LiveServer is available)
julia --project=docs -e 'using LiveServer; serve(dir="docs/build")'
```

### Package Management
```bash
# Install dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Update dependencies
julia --project=. -e 'using Pkg; Pkg.update()'

# Add a JuliaSMLM package as dependency
julia --project=. -e 'using Pkg; Pkg.add("SMLMData")'
```

### Development Setup
```bash
# Develop local packages simultaneously
julia --project=. -e 'using Pkg; Pkg.develop(path="../SMLMData")'
```

## Architecture

### Package Integration Strategy

SMLMAnalysis.jl will integrate the ecosystem packages to provide:

1. **Unified API**: Single import for common workflows
2. **Pipeline Management**: Orchestrate multi-step analysis pipelines
3. **Data Flow**: Handle data transformation between package formats
4. **Configuration**: Centralized parameter management
5. **Workflow Templates**: Pre-configured analysis pipelines

### Typical Analysis Pipelines

#### Standard Localization Pipeline
```
Raw Images → SMLMBoxer (detection) → GaussMLE (fitting) → 
SMLMFrameConnection (connection) → SMLMMetrics (validation)
```

#### Deep Learning Pipeline
```
Raw Images → SMLMDeepFit (detection + fitting) → 
SMLMFrameConnection (connection) → SMLMVis (visualization)
```

#### Simulation-Based Validation
```
SMLMSim (generate data) → Analysis Pipeline → 
SMLMMetrics (compare to ground truth)
```

### Data Flow
- All packages use SMLMData.jl types for interoperability
- Coordinate system: microns for spatial coordinates
- Standard containers: BasicSMLD, SmiteSMLD for compatibility

## Implementation Roadmap

As SMLMAnalysis.jl develops, it should:

1. **Re-export Core Types**: Make SMLMData types available directly
2. **Provide Workflow Functions**: High-level functions for common pipelines
3. **Handle Package Dependencies**: Manage optional dependencies gracefully
4. **Offer Configuration Management**: Centralized parameter handling
5. **Include Example Notebooks**: Demonstrate complete workflows

## Key Design Principles

1. **Composability**: Each package works independently but integrates seamlessly
2. **Performance**: Leverage Julia's speed with GPU support where available
3. **Interoperability**: Maintain compatibility with MATLAB SMITE toolbox formats
4. **Ease of Use**: High-level API hiding complexity while allowing low-level access

## GitHub Resources

- **Organization**: https://github.com/JuliaSMLM
- **Main Repository**: https://github.com/JuliaSMLM/SMLMAnalysis.jl
- **Related**: LidkeLab organization for MATLAB tools and broader microscopy software

## Usage Examples

Future implementation should support workflows like:

```julia
using SMLMAnalysis

# Load and process data
data = load_smlm_data("experiment.h5")
detected = detect_particles(data, method=:boxer)
fitted = fit_emitters(detected, method=:gaussmle)
connected = connect_frames(fitted)
results = calculate_metrics(connected, ground_truth)

# Or use a pre-configured pipeline
results = standard_smlm_pipeline(data)
```

## Testing Strategy

- Unit tests for each integration point
- Integration tests for complete pipelines
- Performance benchmarks comparing different methods
- Validation against established MATLAB implementations