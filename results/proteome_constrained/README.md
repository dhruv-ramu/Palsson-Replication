# Proteome-Constrained Analysis Results

This directory contains the results from the proteome-constrained metabolic analysis of E. coli iJO1366.

## Overview

**Analysis:** Proteome-constrained metabolic modeling with enzyme allocation constraints  
**Reference:** Sánchez et al. 2017 Nature Communications (GECKO framework)  
**Key Finding:** Lactose metabolism is 28x more sensitive to enzyme constraints than glucose  
**Date:** August 2, 2024  

## Files

### Data Files:
- **`proteome_constrained_growth_comparison.json`** - Growth rate comparisons between constrained and unconstrained models
- **`proteome_constrained_flux_data.json`** - Key reaction fluxes under enzyme constraints
- **`proteome_constrained_fva_results.json`** - Flux Variability Analysis results

### Visualizations:
- **`proteome_constrained_analysis.png`** - **Primary visualization** (4-panel analysis)

## Key Results

### Growth Rate Impact
- **Glucose growth reduction:** 2.6% (0.982 → 0.957 1/h)
- **Lactose growth reduction:** 72.7% (1.435 → 0.392 1/h)
- **Constraint sensitivity ratio:** 28:1 (lactose:glucose)

### Enzyme Constraints Applied
| Enzyme | k_cat (1/s) | Max Flux (mmol/gDW/h) | Status |
|--------|-------------|----------------------|--------|
| **HEX1** | 100.0 | 7.20 | Active (42% capacity) |
| **PFK** | 50.0 | 5.14 | Inactive |
| **PYK** | 200.0 | 13.09 | Inactive |
| **CS** | 30.0 | 1.08 | Bottleneck (100% capacity) |
| **MDH** | 150.0 | 15.43 | Bottleneck (100% capacity) |
| **LACZ** | 20.0 | 0.62 | Repressed |

### Flux Variability Analysis
- **CS (Citrate synthase):** Tightly constrained (0.154 variability)
- **MDH (Malate dehydrogenase):** Highly flexible (41.825 variability)
- **Glycolytic enzymes:** Show high flexibility under constraints
- **TCA cycle:** More constrained than glycolysis

## Biological Insights

### Why Lactose is More Constrained
1. **Enzyme Efficiency:** Glucose enzymes have higher k_cat values
2. **Pathway Complexity:** Lactose metabolism requires additional enzymes
3. **Evolutionary Optimization:** Glucose utilization is naturally optimized

### Implications for Diauxic Growth
- **Glucose preference** is reinforced by enzyme efficiency
- **Lactose repression** is energetically favorable
- **Proteome limitations** explain observed growth rate differences

## Data Format

### Growth Comparison JSON
```json
{
  "glucose": {
    "unconstrained": 0.982372,
    "constrained": 0.957111,
    "reduction_percent": 2.6
  },
  "lactose": {
    "unconstrained": 1.435305,
    "constrained": 0.391675,
    "reduction_percent": 72.7
  }
}
```

### Flux Data JSON
```json
{
  "HEX1": 3.016,
  "PFK": 0.0,
  "PYK": 0.0,
  "CS": 1.08,
  "MDH": 15.429,
  "LACZ": 0.0
}
```

### FVA Results JSON
```json
{
  "HEX1": {
    "min_flux": 0.0,
    "max_flux": 7.2,
    "variability": 7.2
  }
}
```

## Technical Notes

- **Enzyme constraints** applied using k_cat values from literature
- **Molecular weights** considered in flux capacity calculations
- **Global protein budget** of 0.5 g protein/gDW
- **FVA performed** at 90% of optimal growth rate 