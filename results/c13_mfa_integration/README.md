# 13C-MFA Integration Results

This directory contains the results from the 13C-MFA integration analysis, comparing FBA predictions with experimental 13C-MFA data using validation methods from Kaste & Shachar-Hill (2023).

## Overview

**Analysis:** FBA predictions vs experimental 13C-MFA data validation  
**Reference:** Kaste & Shachar-Hill 2023 (arXiv:2303.12651)  
**Key Finding:** Strong correlations (0.817) but high MAPE (360.2%) indicating systematic prediction errors  
**Date:** August 2, 2024  

## Files

### Data Files:
- **`c13_mfa_experimental_data.json`** - Experimental 13C-MFA data from literature sources
- **`c13_mfa_fba_predictions.json`** - FBA predictions for all growth conditions
- **`c13_mfa_validation_results.json`** - Comprehensive validation metrics and statistical tests

### Visualizations:
- **`c13_mfa_validation.png`** - **Primary validation figure** (scatter plots and residual analysis)
- **`c13_mfa_summary.png`** - **Summary metrics** (correlation, R², MAPE, growth error by condition)

## Key Results

### Growth Rate Comparison
| Condition | FBA Growth (1/h) | Experimental (1/h) | Error (%) |
|-----------|------------------|-------------------|-----------|
| **Glucose minimal** | 0.982 | 0.850 | 15.6% |
| **Acetate minimal** | 0.247 | 0.420 | 41.1% |
| **Lactose minimal** | 0.428 | 0.350 | 22.4% |

### Validation Metrics Summary
| Condition | Correlation | R² | MAPE (%) | χ² p-value |
|-----------|-------------|----|----------|------------|
| **Glucose minimal** | 0.837 | 0.700 | 169.3% | 0.000 |
| **Acetate minimal** | 0.793 | 0.628 | 332.0% | 0.000 |
| **Lactose minimal** | 0.820 | 0.672 | 579.4% | 0.000 |

### Overall Performance
- **Mean Correlation:** 0.817 (Strong qualitative agreement)
- **Mean R²:** 0.667 (Moderate explanatory power)
- **Mean MAPE:** 360.2% (High quantitative error)
- **Mean Growth Error:** 26.4% (Moderate growth prediction)

## Experimental Data Sources

### Literature References
- **Glucose minimal:** Emmerling et al. 2002, J Bacteriol
- **Acetate minimal:** Nanchen et al. 2006, J Bacteriol
- **Lactose minimal:** Haverkorn van Rijsewijk et al. 2011, PLoS One

### Growth Conditions
- **Glucose minimal:** 0.85 1/h, aerobic
- **Acetate minimal:** 0.42 1/h, aerobic
- **Lactose minimal:** 0.35 1/h, aerobic

## Validation Methods

### Statistical Tests Implemented
1. **Pearson correlation** - Linear relationship strength
2. **Spearman correlation** - Rank-based relationship
3. **Mean Absolute Error (MAE)** - Average absolute difference
4. **Root Mean Square Error (RMSE)** - Standard deviation of errors
5. **Mean Absolute Percentage Error (MAPE)** - Relative error percentage
6. **χ² goodness-of-fit test** - Statistical model validation
7. **R-squared (R²)** - Coefficient of determination

### Kaste & Shachar-Hill Framework
- **χ² test implementation** following recommended methodology
- **Degrees of freedom calculation** for statistical validity
- **p-value interpretation** for model assessment
- **Residual analysis** for error pattern identification

## Data Format

### Experimental Data JSON
```json
{
  "glucose_minimal": {
    "source": "Emmerling et al. 2002, J Bacteriol",
    "conditions": "Glucose minimal medium, aerobic",
    "growth_rate": 0.85,
    "fluxes": {
      "HEX1": 8.5,
      "PFK": 8.5,
      "PYK": 8.5,
      ...
    }
  }
}
```

### Validation Results JSON
```json
{
  "glucose_minimal": {
    "n_reactions": 19,
    "correlation": 0.837,
    "r_squared": 0.700,
    "mape": 169.3,
    "chi_square": 103.380,
    "chi_square_p_value": 0.000,
    ...
  }
}
```

## Model Validation Assessment

### Strengths
✅ **Strong correlations** (0.817 average) indicate good qualitative agreement  
✅ **Consistent directional predictions** across all conditions  
✅ **Growth rate predictions** within reasonable error range (26.4% average)  
✅ **Comprehensive validation** using multiple metrics  

### Weaknesses
❌ **High MAPE values** (360.2% average) indicate systematic prediction errors  
❌ **χ² tests** show significant deviations from experimental data  
❌ **Flux magnitude predictions** often inaccurate  
❌ **Condition-specific biases** in predictions  

## Biological Interpretation

### Glucose Metabolism
- **FBA overpredicts growth** (15.6% error)
- **Strong correlation** but high MAPE suggests systematic scaling issues
- **Glycolytic fluxes** well-predicted qualitatively but not quantitatively

### Acetate Metabolism
- **FBA underpredicts growth** (41.1% error)
- **Glyoxylate cycle** not fully captured by model
- **Anaplerotic reactions** may be missing or incorrectly constrained

### Lactose Metabolism
- **Highest MAPE** (579.4%) indicates major prediction issues
- **Lactose operon regulation** not captured in FBA
- **Enzyme constraints** may be needed for accurate predictions

## Technical Notes

- **3 growth conditions** analyzed with comprehensive validation
- **Multiple validation metrics** calculated efficiently
- **Publication-ready visualizations** generated automatically
- **Robust error handling** for missing reactions
- **Statistical significance** assessed with proper p-values
- **χ² test implementation** follows Kaste & Shachar-Hill methodology 