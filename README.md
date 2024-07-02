# Julia-based FEM Implementation for Axisymmetrical Continuum Elastic Modeling of Bilayers

This repository contains a Julia-based implementation of a finite element method (FEM) for the axisymmetrical continuum elastic modeling of bilayers. Building on the work of Rolf Ryham (2016), this implementation introduces volume and compositional constraints for each leaflet and encapsulated water using a modified augmented Lagrangian method. The modules are designed with a geometrical approach, facilitating easy setup of initial structures and code modification. Further details can be found in the Supplementary Information of my forthcoming paper.

## Background

Bilayers are a crucial class of biological membranes involved in various biological processes. Continuum elastic modeling of bilayers is a powerful tool for understanding their mechanical properties and the effects of external stimuli. The FEM implementation presented here extends Ryham's axisymmetrical model by incorporating volume constraints for each leaflet and encapsulated water through a modified augmented Lagrangian method.

## Features

This FEM implementation includes the following features:

- **Axisymmetrical continuum elastic modeling of bilayers**
- **Composition fields over the membrane**: Both predefined domains (non-varying) and general domains with mixing free energy
- **Volume and compositional constraints** for leaflets and encapsulated water via a modified augmented Lagrangian method
- **Geometrical approach**: Facilitates easy setup of initial structures and general code modifications
- **Remeshing scheme**: Provides resolution control
- **Free energy path**: Implemented using the modified Nudged Elastic Band (mNEB) method
- **State management**: Save and load states

## Usage

Example scripts used to create figures for the paper are provided. The general workflow involves creating a file with all necessary predefined parameters (as shown in most example Julia scripts), including `ElasticTool.jl`, building the required elements (coordinates, volume, and remesh patches) for the mesh (see `Builder.jl` for an example), and running the optimization.

## File Descriptions

- **ElasticTool.jl**: Loads all the required modules
- **Builder.jl**: Provides functions for constructing the initial configuration with the geometrical elements used in the FEM scheme
- **Diff.jl**: Handles gradient calculations
- **Energetics.jl**: Contains functionalities related to energy calculations over the mesh
- **Geometry.jl**: Defines the geometrical structures that form the model
- **Remesh.jl**: Provides remeshing functionality
- **SaveStates.jl**: Manages saving, loading, and modifying configurations

## Citation

When using this code or some modification of it. please cite the paper (arxiv if needed ***):

doi ***

As the model is based on Ryham's model, please cite his papers as well:

doi ***
