# Phenotype Phase Plane Analysis

This repository contains a Python script for performing phenotype phase plane analysis on metabolic models, specifically examining the relationship between substrate uptake and oxygen consumption in E. coli iJO1366.

## Overview

Phenotype phase plane analysis is a powerful tool in metabolic engineering that helps understand how the growth rate of an organism changes as a function of two different substrate uptake rates. This analysis can reveal optimal growth conditions and metabolic trade-offs.

## Files

- `simple_phpp_analysis.py` - Main analysis script
- `bigg_models/iJO1366.xml` - E. coli metabolic model (SBML format)
- `glucose_simple_phpp.csv` - Results from glucose analysis (generated)
- `acetate_simple_phpp.csv` - Results from acetate analysis (generated)
- `glucose_simple_phpp.png` - Visualization plot for glucose (generated)
- `acetate_simple_phpp.png` - Visualization plot for acetate (generated)

## Requirements

Install the required Python packages:

```bash
pip install cobra numpy pandas matplotlib
```

## Usage

Run the analysis script:

```bash
python simple_phpp_analysis.py
```

## What the Script Does

1. **Loads the metabolic model** from SBML file
2. **Tests model growth capabilities** with basic constraints
3. **Performs phenotype phase plane analysis** for glucose vs oxygen uptake
4. **Performs phenotype phase plane analysis** for acetate vs oxygen uptake
5. **Creates visualizations** showing growth rate as a function of substrate uptake
6. **Saves results** to CSV and PNG files

## Output

The script generates:

- **Console output**: Summary statistics and growth information
- **CSV files**: Raw data with substrate uptake, oxygen uptake, and growth rate
- **PNG files**: Visualization of the phenotype phase planes

## Key Results

The analysis reveals:
- **Growth capabilities**: How well the model can grow on different substrates
- **Optimal conditions**: The substrate and oxygen uptake rates that maximize growth
- **Metabolic trade-offs**: How growth rate changes with different substrate combinations

## Technical Details

The script uses:
- **COBRApy**: For metabolic model manipulation and optimization
- **NumPy**: For numerical computations
- **Pandas**: For data handling and analysis
- **Matplotlib**: For visualization

## Example Results

The analysis shows that the iJO1366 model can achieve:
- **Glucose growth**: Up to 0.899 1/h growth rate
- **Acetate growth**: Up to 0.899 1/h growth rate
- **176 data points** analyzed for each substrate
- **Clear phenotype phase planes** showing optimal growth regions
