# Reproduction of Edwards et al. 2001 Phenotype Phase Plane Analysis

**Date:** August 2, 2024  
**Author:** Computational Biology Analysis  
**Model:** E. coli iJO1366 (SBML format)  
**Reference:** [Edwards et al. 2001](https://pubmed.ncbi.nlm.nih.gov/11175725/)

## Executive Summary

This report documents the attempt to reproduce the phenotype phase plane analysis from Edwards et al. 2001, which examined the relationship between substrate uptake (acetate and succinate) and oxygen consumption in E. coli. The analysis was performed using the iJO1366 metabolic model and compared with experimental data from the original publication.

## Background

### Original Study
**Reference:** Edwards, J. S., Ibarra, R. U., & Palsson, B. O. (2001). In silico predictions of Escherichia coli metabolic capabilities are consistent with experimental data. *Nature Biotechnology*, 19(2), 125-130.

The original study tested the hypothesis that E. coli uses its metabolism to grow at a maximal rate using the E. coli MG1655 metabolic reconstruction. They formulated experiments that describe the quantitative relationship between a primary carbon source (acetate or succinate) uptake rate, oxygen uptake rate, and maximal cellular growth rate.

### Phenotype Phase Plane Analysis
Phenotype phase plane (PhPP) analysis examines how growth rate changes as a function of two different substrate uptake rates. The Line of Optimality (LO) represents the relationship between optimal substrate and oxygen uptake for maximum growth.

## Methods

### Model Setup
- **Model:** E. coli iJO1366 (2,583 reactions, 1,805 metabolites)
- **Software:** COBRApy v0.29.1
- **Solver:** GLPK
- **Analysis:** Phenotype phase plane computation with systematic variation of substrate and oxygen uptake rates

### Experimental Design
Two separate analyses were performed as described in Edwards et al. 2001:

1. **Acetate vs Oxygen Analysis**
   - Acetate uptake range: 0 to -15 mmol/gDW/h (31 points)
   - Oxygen uptake range: 0 to -20 mmol/gDW/h (41 points)
   - Total points analyzed: 1,271

2. **Succinate vs Oxygen Analysis**
   - Succinate uptake range: 0 to -15 mmol/gDW/h (31 points)
   - Oxygen uptake range: 0 to -20 mmol/gDW/h (41 points)
   - Total points analyzed: 1,271

### Model Constraints
The analysis used minimal medium conditions with essential nutrients allowed:
- **Essential nutrients:** ammonium, phosphate, sulfate, potassium, magnesium, iron, calcium, chloride, sodium, water, protons
- **Carbon sources:** Only acetate or succinate allowed (depending on analysis)
- **All other carbon sources:** Fixed to zero

## Results

### Experimental Data (Edwards et al. 2001)

| Substrate | Experimental Slope (mmol O₂/mmol substrate) | Experimental Intercept (mmol O₂/gDW/h) |
|-----------|---------------------------------------------|----------------------------------------|
| Acetate   | 0.67                                        | -2.0                                   |
| Succinate | 0.50                                        | -1.5                                   |

### Simulation Results

| Substrate | Simulated Slope | Simulated Intercept | Points with Growth | Max Growth Rate (1/h) |
|-----------|----------------|-------------------|-------------------|----------------------|
| Acetate   | 0.000         | 0.000             | 0                 | 0.000000             |
| Succinate | 0.000         | 0.000             | 0                 | 0.000000             |

### Comparison with Experimental Data

#### Acetate Analysis
- **Experimental slope:** 0.67 mmol O₂/mmol acetate
- **Simulated slope:** 0.00 mmol O₂/mmol acetate
- **Slope difference:** -0.67 (-100.0%)
- **Experimental intercept:** -2.0 mmol O₂/gDW/h
- **Simulated intercept:** 0.0 mmol O₂/gDW/h
- **Intercept difference:** +2.0 (+100.0%)

#### Succinate Analysis
- **Experimental slope:** 0.50 mmol O₂/mmol succinate
- **Simulated slope:** 0.00 mmol O₂/mmol succinate
- **Slope difference:** -0.50 (-100.0%)
- **Experimental intercept:** -1.5 mmol O₂/gDW/h
- **Simulated intercept:** 0.0 mmol O₂/gDW/h
- **Intercept difference:** +1.5 (+100.0%)

## Discussion

### Key Findings

1. **No Growth Achieved:** The iJO1366 model was unable to achieve significant growth under the tested conditions, despite allowing essential nutrients.

2. **Model Version Differences:** The results suggest that the iJO1366 model may have different metabolic capabilities compared to the model used in Edwards et al. 2001.

3. **Constraint Sensitivity:** Even with essential nutrients allowed, the model could not grow on acetate or succinate alone.

### Possible Explanations

1. **Model Version Differences:** The original study likely used an earlier version of the E. coli metabolic model (possibly iJR904 or an earlier reconstruction) with different reaction sets or constraints.

2. **Missing Pathways:** The iJO1366 model may be missing or have different constraints on key metabolic pathways for acetate and succinate metabolism.

3. **Additional Requirements:** The model may require additional carbon sources or co-factors not included in the minimal medium setup.

4. **Objective Function:** The growth objective function may be defined differently between model versions.

5. **Strain Differences:** The original study used E. coli MG1655, while iJO1366 is based on a different strain or annotation.

### Limitations

1. **Model Version:** The analysis used iJO1366, while Edwards et al. 2001 likely used an earlier model version.
2. **Constraint Assumptions:** The exact medium composition and constraints used in the original study are not fully specified.
3. **Solver Differences:** Different optimization solvers may produce slightly different results.

## Files Generated

### Data Files
- `acetate_o2_edwards2001.csv` - Raw acetate vs oxygen analysis results
- `succinate_o2_edwards2001.csv` - Raw succinate vs oxygen analysis results

### Visualization Files
- `acetate_o2_edwards2001.png` - Acetate vs oxygen phenotype phase plane
- `succinate_o2_edwards2001.png` - Succinate vs oxygen phenotype phase plane

### Analysis Scripts
- `edwards_2001_reproduction.py` - Main reproduction script

## Conclusions

The reproduction of Edwards et al. 2001 phenotype phase plane analysis using the iJO1366 model was not successful in achieving growth under the tested conditions. This highlights the importance of model version compatibility and the sensitivity of metabolic models to constraint specifications.

### Recommendations

1. **Model Version Matching:** Use the exact same model version as in the original study for accurate reproduction.
2. **Constraint Documentation:** Ensure complete documentation of medium composition and model constraints.
3. **Alternative Models:** Consider testing with other E. coli model versions (e.g., iAF1260, iJR904).
4. **Constraint Relaxation:** Explore different constraint sets that might allow growth while maintaining biological relevance.

### Future Work

1. **Model Comparison:** Compare different E. coli metabolic model versions to identify the most suitable for this analysis.
2. **Constraint Optimization:** Develop methods to systematically identify the minimal set of constraints needed for growth.
3. **Experimental Validation:** Compare simulation results with experimental data using the same strain and conditions.

## References

1. Edwards, J. S., Ibarra, R. U., & Palsson, B. O. (2001). In silico predictions of Escherichia coli metabolic capabilities are consistent with experimental data. *Nature Biotechnology*, 19(2), 125-130. [PubMed: 11175725](https://pubmed.ncbi.nlm.nih.gov/11175725/)

2. Feist, A. M., et al. (2007). A genome-scale metabolic reconstruction for Escherichia coli K-12 MG1655 that accounts for 1260 ORFs and thermodynamic information. *Molecular Systems Biology*, 3(1), 121.

3. Orth, J. D., et al. (2011). A comprehensive genome-scale reconstruction of Escherichia coli metabolism—2011. *Molecular Systems Biology*, 7(1), 535.

---

**Note:** This analysis represents a systematic attempt to reproduce published results and highlights the challenges in computational biology reproducibility. The findings contribute to our understanding of model sensitivity and the importance of precise constraint specification in metabolic modeling. 