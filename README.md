# Phenotype Phase Plane Analysis

This repository contains a Python script for performing phenotype phase plane analysis on metabolic models, specifically examining the relationship between acetate uptake and oxygen consumption in E. coli iJO1366.

## Overview

Phenotype phase plane analysis is a powerful tool in metabolic engineering that helps understand how the growth rate of an organism changes as a function of two different substrate uptake rates. This analysis can reveal optimal growth conditions and metabolic trade-offs.

## Files

- `phenotype_phase_plane_analysis.py` - Main analysis script
- `bigg_models/iJO1366.xml` - E. coli metabolic model (SBML format)
- `phpp_ac_o2.csv` - Results from the analysis (generated)
- `phenotype_phase_plane.png` - Visualization plot (generated)

## Requirements

Install the required Python packages:

```bash
pip install cobra numpy pandas matplotlib
```

## Usage

Run the analysis script:

```bash
python phenotype_phase_plane_analysis.py
```

## What the Script Does

1. **Loads the metabolic model** from SBML file
2. **Sets up exchange reactions** for acetate, succinate, and oxygen
3. **Fixes other carbon sources** to zero to focus on acetate metabolism
4. **Computes phenotype phase plane** by systematically varying acetate and oxygen uptake rates
5. **Finds the Line of Optimality** (LO) - the relationship between optimal acetate and oxygen uptake
6. **Creates visualizations** showing growth rate as a function of substrate uptake
7. **Saves results** to CSV and PNG files

## Output

The script generates:

- **Console output**: Summary statistics and best growth conditions
- **CSV file** (`phpp_ac_o2.csv`): Raw data with acetate uptake, oxygen uptake, and growth rate
- **PNG file** (`phenotype_phase_plane.png`): Visualization of the phenotype phase plane

## Key Results

The analysis reveals:
- **Line of Optimality slope**: How much oxygen is needed per unit of acetate for optimal growth
- **Best growth conditions**: The acetate and oxygen uptake rates that maximize growth
- **Metabolic trade-offs**: How growth rate changes with different substrate combinations

## Customization

You can modify the script to:
- Change the substrate ranges (modify `AC_LIMITS` and `O2_LIMITS`)
- Analyze different substrates (modify the exchange reaction IDs)
- Adjust the optimality threshold (modify the `threshold` parameter in `find_line_of_optimality`)

## Troubleshooting

If you encounter issues:

1. **Model loading errors**: Ensure the SBML file path is correct
2. **Zero growth rates**: Check if the model constraints are too restrictive
3. **Memory issues**: Reduce the number of points in the analysis by modifying the limits

## Technical Details

The script uses:
- **COBRApy**: For metabolic model manipulation and optimization
- **NumPy**: For numerical computations
- **Pandas**: For data handling and analysis
- **Matplotlib**: For visualization

## References

- Edwards, J. S., & Palsson, B. O. (2000). The Escherichia coli MG1655 in silico metabolic genotype: its definition, characteristics, and capabilities. *Proceedings of the National Academy of Sciences*, 97(10), 5528-5533.
- Feist, A. M., et al. (2007). A genome-scale metabolic reconstruction for Escherichia coli K-12 MG1655 that accounts for 1260 ORFs and thermodynamic information. *Molecular Systems Biology*, 3(1), 121.
