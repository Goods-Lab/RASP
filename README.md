# Randomized Spatial PCA (RASP)


<img src="RASP/figures/figure_1_A.pdf">

## Description

Here, we present Randomized Spatial PCA (RASP), a novel spatially aware dimensionality reduction method for spatial transcriptomics (ST) data. 
RASP is designed to be orders-of-magnitude faster than existing techniques, scale to ST data with hundreds of thousands of locations, support the 
flexible integration of non-transcriptomic covariates, and enable the reconstruction of de-noised and spatially smoothed expression values for individual genes. 
To achieve these goals, RASP uses a randomized two-stage principal component analysis (PCA) framework that leverages sparse matrix operations and configurable spatial smoothing.

## Features

- **High-Speed Performance**: RASP is optimized for fast processing of large spatial transcriptomics datasets.
- **Flexible Integration**: Seamlessly integrates non-transcriptomic covariates into the analysis.
- **Spatially Smoothed Values**: Produces reconstructed expression values that account for spatial context.
- **User-Friendly**: Designed to be accessible for researchers in spatial biology.

## Requirements

To run RASP, you need the following Python packages:

- numpy
- pandas
- scanpy
- squidpy
- matplotlib
- scipy
- sklearn
- rpy2
- mclustpy
- igraph

## Citation
If you use RASP in your research, please cite the following preprint: 
