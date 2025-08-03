# Omics Data Directory

This directory contains all omics datasets for multi-omics integration with the metabolic network.

## Overview

**Purpose:** Centralized repository for all omics data (transcriptomics, proteomics, metabolomics, fluxomics)  
**Integration Target:** E. coli iJO1366 metabolic network  
**Data Sources:** Literature, experimental results, and computational predictions  

## Available Datasets

### 1. Fluxomics (13C-MFA)
- **Source:** `../results/c13_mfa_integration/`
- **Data:** Experimental 13C-MFA flux measurements
- **Conditions:** Glucose, acetate, lactose minimal media
- **Format:** JSON with reaction fluxes
- **Reference:** Emmerling et al. 2002, Nanchen et al. 2006, Haverkorn van Rijsewijk et al. 2011

### 2. Proteomics (Enzyme Constraints)
- **Source:** `../results/proteome_constrained/`
- **Data:** Enzyme abundance and k_cat constraints
- **Conditions:** Glucose and lactose metabolism
- **Format:** JSON with enzyme flux capacities
- **Reference:** Sánchez et al. 2017 Nature Communications

### 3. Genomics (Gene Essentiality)
- **Source:** `../results/gene_essentiality/`
- **Data:** Gene essentiality predictions vs experimental data
- **Conditions:** Single-gene knockouts
- **Format:** CSV with gene essentiality scores
- **Reference:** Orth et al. 2011

### 4. Transcriptomics (Gene Expression)
- **Status:** Need to source from literature
- **Target:** E. coli gene expression data under different growth conditions
- **Format:** CSV/JSON with gene expression levels
- **Priority:** High - needed for complete multi-omics integration

### 5. Metabolomics (Metabolite Concentrations)
- **Status:** Need to source from literature
- **Target:** E. coli metabolite concentration data
- **Format:** CSV/JSON with metabolite levels
- **Priority:** High - needed for complete multi-omics integration

## Data Integration Plan

### Phase 1: Data Collection (Current)
- [x] Identify existing omics data sources
- [x] Create centralized omics directory
- [ ] Source transcriptomics data from literature
- [ ] Source metabolomics data from literature
- [ ] Standardize all data formats

### Phase 2: Data Standardization
- [ ] Convert all datasets to consistent format
- [ ] Create mapping files for network integration
- [ ] Validate data quality and completeness
- [ ] Document data preprocessing steps

### Phase 3: Network Integration
- [ ] Map omics data to metabolic network nodes
- [ ] Update graph construction pipeline
- [ ] Add omics features to node attributes
- [ ] Validate integration results

## File Structure

```
data/omics/
├── README.md                    # This file
├── fluxomics/                   # 13C-MFA flux data
│   ├── experimental_data.json
│   ├── fba_predictions.json
│   └── validation_results.json
├── proteomics/                  # Enzyme constraint data
│   ├── enzyme_constraints.json
│   ├── flux_data.json
│   └── fva_results.json
├── genomics/                    # Gene essentiality data
│   ├── gene_essentiality.csv
│   ├── comparison_with_orth.csv
│   └── metrics.json
├── transcriptomics/             # Gene expression data (to be sourced)
│   └── README.md
├── metabolomics/                # Metabolite concentration data (to be sourced)
│   └── README.md
└── integration/                 # Integration scripts and mapping files
    ├── data_mapper.py
    ├── format_converter.py
    └── network_integrator.py
```

## Data Sources to Source

### Transcriptomics Data
1. **E. coli gene expression under different carbon sources**
   - Glucose vs acetate vs lactose
   - Aerobic vs anaerobic conditions
   - Growth phase-dependent expression

2. **Potential Sources:**
   - GEO database (Gene Expression Omnibus)
   - ArrayExpress database
   - Literature: Covert et al. 2004, Ishii et al. 2007
   - E. coli gene expression compendiums

### Metabolomics Data
1. **E. coli metabolite concentrations**
   - Central carbon metabolism intermediates
   - Amino acid pools
   - Nucleotide pools
   - Co-factor levels

2. **Potential Sources:**
   - MetaboLights database
   - Literature: Bennett et al. 2009, Ishii et al. 2007
   - E. coli metabolomics studies

## Data Format Standards

### JSON Format for Omics Data
```json
{
  "dataset_info": {
    "name": "dataset_name",
    "source": "literature_reference",
    "conditions": "growth_conditions",
    "date": "YYYY-MM-DD"
  },
  "data": {
    "entity_id": {
      "value": numeric_value,
      "unit": "unit_string",
      "condition": "condition_string"
    }
  }
}
```

### CSV Format for Omics Data
```csv
entity_id,value,unit,condition,replicate
gene_001,1.5,log2_fold_change,glucose_minimal,1
gene_002,0.8,log2_fold_change,glucose_minimal,1
```

## Integration Mapping

### Network Node Types
1. **Metabolites** → Metabolomics data
2. **Reactions** → Fluxomics data, Proteomics data
3. **Genes** → Transcriptomics data, Genomics data

### Mapping Strategy
- **Metabolite mapping:** Use BiGG IDs (e.g., "glc__D_c")
- **Reaction mapping:** Use BiGG IDs (e.g., "HEX1")
- **Gene mapping:** Use BiGG gene IDs (e.g., "b0002")

## Quality Control

### Data Validation Checklist
- [ ] Data completeness (no missing values)
- [ ] Data consistency (same units, conditions)
- [ ] Biological plausibility (reasonable ranges)
- [ ] Source reliability (peer-reviewed literature)
- [ ] Mapping accuracy (correct entity IDs)

### Integration Validation
- [ ] Feature correlation analysis
- [ ] Network topology preservation
- [ ] Biological pathway enrichment
- [ ] Cross-validation with known relationships

## Next Steps

1. **Source transcriptomics data** from GEO/ArrayExpress databases
2. **Source metabolomics data** from MetaboLights/literature
3. **Create data standardization scripts**
4. **Implement network integration pipeline**
5. **Validate multi-omics integration results**

## References

- **13C-MFA:** Emmerling et al. 2002, Nanchen et al. 2006, Haverkorn van Rijsewijk et al. 2011
- **Proteomics:** Sánchez et al. 2017 Nature Communications
- **Genomics:** Orth et al. 2011
- **Transcriptomics:** Covert et al. 2004, Ishii et al. 2007
- **Metabolomics:** Bennett et al. 2009, Ishii et al. 2007 