PhyloGeographeR
===========

A set of R scripts with utility functions for Bayesian phylogeography with [BEAST](https://github.com/beast-dev/beast-mcmc).

The plan is to have three basic modules 

- elicitatoR : this module helps create design matrices for generalised linear models (GLM) and set up marginal likelihood analyses for different predictors;

- annotatorR : this module contains functions to summarise and visualise the output from BEAST;

- simulatoR : this module will simulate traits down a tree using CTMC-based phylogeographic models;

The objective here is to combine the capabilities from geiger, phylosim and spdep/maptools to incorporate geographic information, such as distances, availability of roads and neighborhood structure into the formulation of CMTC-based discrete trait simulation. This simulation tool shall be of great value in generating synthetic data to evaluate phylogeography reconstruction methods/software.
