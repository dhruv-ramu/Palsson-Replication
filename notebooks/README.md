# Computational Biology Analysis Notebooks

This directory contains Jupyter notebooks for metabolic modeling and analysis, organized into logical sections for clarity and progression.

## Directory Structure

### 01_Basic_Analysis/
Fundamental metabolic modeling analyses:
- `01_Phenotype_Phase_Plane.ipynb`: PhPP analysis for substrate uptake vs growth
- `02_Gene_Essentiality.ipynb`: Gene knockout analysis and validation
- `03_Flux_Analysis.ipynb`: Basic flux analysis and visualization

### 02_Advanced_Analysis/
More sophisticated metabolic analyses:
- `01_Dynamic_FBA.ipynb`: Dynamic FBA for diauxic growth
- `02_13C_MFA_Integration.ipynb`: Integration with 13C-MFA data
- `03_Proteome_Constraints.ipynb`: Proteome-constrained analysis

### 03_Model_Improvements/
Model enhancement and customization:
- `01_Custom_Model_Development.ipynb`: Building custom metabolic models
- `02_Model_Validation.ipynb`: Validation against experimental data
- `03_Advanced_Constraints.ipynb`: Implementing advanced constraints

### 04_Utilities/
Helper notebooks and utilities:
- `01_Model_Download.ipynb`: BiGG model download and setup
- `02_Data_Processing.ipynb`: Data processing utilities
- `03_Visualization.ipynb`: Reusable visualization functions

## Usage

1. Start with notebooks in `01_Basic_Analysis/` to understand fundamental concepts
2. Progress to `02_Advanced_Analysis/` for more complex analyses
3. Use `03_Model_Improvements/` to enhance model accuracy
4. Reference `04_Utilities/` as needed for common tasks

## Dependencies

Required Python packages:
- cobra
- numpy
- pandas
- matplotlib
- seaborn
- scipy

## Data

- Model files are stored in `bigg_models/`
- Results are saved to `results/`
- Each analysis saves its outputs in a dedicated subfolder

## References

- Edwards et al. (2001) - Phenotype predictions
- Orth et al. (2011) - Gene essentiality
- Mahadevan et al. (2002) - Dynamic FBA
- Palsson Lab publications and methods