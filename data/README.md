# Data Directory

This directory contains input data, notebooks, and reference materials for the metabolic modeling analyses.

## Files

### `GlucoseMinimalMedium.ipynb`
- **Purpose**: Jupyter notebook for glucose minimal medium analysis
- **Content**: Initial exploration of E. coli growth on glucose
- **Status**: Reference material for basic metabolic modeling

## Related Data

### Metabolic Models
Metabolic models are stored in the `../bigg_models/` directory:
- **iJO1366.xml**: E. coli metabolic model in SBML format
- **iJO1366.json**: E. coli metabolic model in JSON format
- **Compressed versions**: .gz files for efficient storage

### Model Sources
Models are downloaded from the BiGG database using `../scripts/bigg_download.py`

## Usage

This directory serves as a repository for:
- Input data files
- Jupyter notebooks for exploration
- Reference materials
- Documentation of data sources

## Data Sources

- **BiGG Database**: Metabolic models (bigg.ucsd.edu)
- **Edwards et al. 2001**: Experimental data for comparison
- **COBRApy**: Software framework for constraint-based modeling 