# Scripts Directory

This directory contains all Python scripts for metabolic modeling and phenotype phase plane analysis.

## Files

### `edwards_2001_reproduction.py`
- **Purpose**: Reproduces the phenotype phase plane analysis from Edwards et al. 2001
- **Function**: Analyzes acetate vs oxygen and succinate vs oxygen uptake relationships
- **Output**: Generates CSV data files and PNG visualizations in `../results/edwards_2001/`

### `simple_phpp_analysis.py`
- **Purpose**: Basic phenotype phase plane analysis for glucose and acetate
- **Function**: Demonstrates working PhPP analysis with common substrates
- **Output**: Generates CSV data files and PNG visualizations in `../results/simple_phpp/`

### `bigg_download.py`
- **Purpose**: Downloads metabolic models from the BiGG database
- **Function**: Retrieves E. coli iJO1366 model in SBML format
- **Output**: Saves models to `../bigg_models/`

## Usage

All scripts can be run from the project root directory:

```bash
python scripts/script_name.py
```

## Dependencies

- COBRApy
- NumPy
- Pandas
- Matplotlib
- SciPy

See `../environment.yml` for complete dependency list. 