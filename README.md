# E. coli Metabolic Modeling: Phenotype Phase Plane Analysis

This repository contains computational biology analyses focusing on phenotype phase plane (PhPP) analysis of E. coli metabolism, including reproduction of published results from Edwards et al. 2001.

## 📁 Project Structure

```
Palsson/
├── scripts/                    # Python analysis scripts
│   ├── edwards_2001_reproduction.py
│   ├── simple_phpp_analysis.py
│   └── bigg_download.py
├── data/                       # Input data and notebooks
│   └── GlucoseMinimalMedium.ipynb
├── results/                    # Analysis results and visualizations
│   ├── edwards_2001/          # Edwards et al. 2001 reproduction results
│   └── simple_phpp/           # Basic PhPP analysis results
├── reports/                    # Comprehensive analysis reports
│   └── Edwards_2001_Reproduction_Report.md
├── docs/                       # Additional documentation
├── bigg_models/               # Metabolic model files
│   ├── iJO1366.xml
│   └── iJO1366.json
├── environment.yml            # Conda environment specification
└── README.md                  # This file
```

## 🚀 Quick Start

### 1. Environment Setup
```bash
# Create conda environment
conda env create -f environment.yml
conda activate palsson

# Or install dependencies manually
pip install cobra numpy pandas matplotlib scipy
```

### 2. Download Metabolic Model
```bash
python scripts/bigg_download.py
```

### 3. Run Basic Analysis
```bash
python scripts/simple_phpp_analysis.py
```

### 4. Run Edwards et al. 2001 Reproduction
```bash
python scripts/edwards_2001_reproduction.py
```

## 📊 Analyses Included

### Simple Phenotype Phase Plane Analysis
- **Purpose**: Demonstrate working PhPP analysis
- **Substrates**: Glucose and acetate
- **Output**: Growth rate vs substrate and oxygen uptake
- **Results**: `results/simple_phpp/`

### Edwards et al. 2001 Reproduction
- **Purpose**: Reproduce published phenotype phase plane analysis
- **Reference**: [Edwards et al. 2001](https://pubmed.ncbi.nlm.nih.gov/11175725/)
- **Substrates**: Acetate and succinate vs oxygen
- **Comparison**: In silico vs experimental slopes
- **Results**: `results/edwards_2001/`
- **Report**: `reports/Edwards_2001_Reproduction_Report.md`

## 📈 Key Findings

### Simple PhPP Analysis
- ✅ **Successful growth** on glucose and acetate
- ✅ **Clear phenotype phase planes** with growth gradients
- ✅ **Working demonstration** of PhPP methodology

### Edwards et al. 2001 Reproduction
- ⚠️ **No growth achieved** with iJO1366 model under tested conditions
- 📊 **Model version differences** identified as key factor
- 🔍 **Comprehensive analysis** of reproducibility challenges
- 📝 **Detailed documentation** of findings and limitations

## 📁 Directory Details

### `scripts/`
Contains all Python analysis scripts:
- `edwards_2001_reproduction.py` - Main reproduction script
- `simple_phpp_analysis.py` - Basic PhPP analysis
- `bigg_download.py` - Model download utility

### `results/`
Organized by analysis type:
- `edwards_2001/` - Edwards et al. 2001 reproduction results
- `simple_phpp/` - Basic PhPP analysis results

### `reports/`
Comprehensive analysis reports in Markdown format.

### `data/`
Input data, notebooks, and reference materials.

## 🔬 Scientific Context

This project explores:
- **Phenotype Phase Plane Analysis**: How growth rate changes with substrate uptake
- **Metabolic Modeling**: Constraint-based reconstruction and analysis
- **Reproducibility**: Challenges in computational biology
- **Model Validation**: Comparison with experimental data

## 📚 References

1. Edwards, J. S., Ibarra, R. U., & Palsson, B. O. (2001). In silico predictions of Escherichia coli metabolic capabilities are consistent with experimental data. *Nature Biotechnology*, 19(2), 125-130.

2. Feist, A. M., et al. (2007). A genome-scale metabolic reconstruction for Escherichia coli K-12 MG1655. *Molecular Systems Biology*, 3(1), 121.

3. Orth, J. D., et al. (2011). A comprehensive genome-scale reconstruction of Escherichia coli metabolism. *Molecular Systems Biology*, 7(1), 535.

## 🤝 Contributing

This is a research repository. For questions or contributions, please refer to the documentation in the `docs/` directory.

## 📄 License

This project is for research purposes. Please cite the original papers when using these analyses.
