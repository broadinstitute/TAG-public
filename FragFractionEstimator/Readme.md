# cfTF: Cell-Free Tumor Fraction Estimation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

> A next-generation tool for estimating tumor fraction from low-pass whole genome sequencing of cell-free DNA, optimized for clinical use and pediatric oncology.

---

## Overview

**cfTF** is a computational tool designed to accurately estimate the fraction of tumor-derived DNA in cell-free DNA (cfDNA) samples from patients with cancer. Unlike existing tools, cfTF uses a multi-modal approach combining fragmentomic features, copy number analysis, and rigorous statistical modeling to provide highly accurate tumor fraction estimates with proper uncertainty quantification.

### Key Features

- **Multi-modal feature integration**: Fragment sizes, end motifs, nucleosome spacing (FFT-based), and copy number alterations
- **Interpretability-first design**: Classical statistics → gradient boosting → transformers (only if needed)
- **Rigorous uncertainty quantification**: Proper confidence intervals using multiple methods (bootstrap, conformal prediction, Bayesian inference)
- **Clinical-grade QC**: Comprehensive, interpretable quality control metrics suitable for clinical deployment
- **Longitudinal modeling**: Joint analysis of multiple timepoints to track tumor dynamics
- **Flexible input modes**: Works with cfDNA alone, or can leverage matched tumor/normal tissue when available
- **Optimized for 3× coverage**: Designed to work well at moderate coverage depths (0.5-5×)
- **Pediatric-focused**: Tailored for pediatric cancer samples with age-appropriate normal references

---

## Motivation

### Why Build a New Tool?

Existing tools for tumor fraction estimation have significant limitations:

**ichorCNA:**
- ✅ Fast, open-source, works at ultra-low coverage (0.1×)
- ❌ Requires detectable copy number alterations (fails on diploid tumors)
- ❌ Has issues at higher coverage depths (>1×) due to polyclonal event filtering
- ❌ Limited sensitivity at low tumor fractions (<3-5%)
- ❌ No confidence intervals on predictions
- ❌ Minimal QC metrics

**DELFI-TF:**
- ✅ Mutation-independent, works across tumor types
- ✅ Uses fragmentomic features (fragment sizes, chromatin patterns)
- ❌ Proprietary/commercial (not open-source)
- ❌ Requires higher coverage (~6×) for optimal performance
- ❌ Black-box machine learning (limited interpretability)
- ❌ Not optimized for pediatric cancers

### Our Approach

cfTF addresses these limitations by:

1. **Combining the best of both approaches**: CNA-based features (like ichorCNA) + fragmentomic features (like DELFI-TF)
2. **Working efficiently at 3× coverage**: Our operational sequencing depth
3. **Prioritizing interpretability**: Understand *why* the model makes each prediction
4. **Providing rigorous uncertainty estimates**: Know when predictions are reliable
5. **Building for clinical deployment**: Comprehensive QC, clear documentation, validated performance
6. **Optimizing for pediatric oncology**: Age-appropriate references, pediatric tumor signatures

---

## Scientific Approach

### Feature Engineering

cfTF extracts multiple complementary signals from low-pass WGS data:

#### 1. **Fragment Size Features**
- Overall size distributions (100-220 bp range)
- Short/long fragment ratios at multiple thresholds
- Modal positions (nucleosome-associated peaks at ~167 bp)
- Size entropy and diversity metrics
- Window-specific size profiles (5 Mb bins)

#### 2. **Fragment End Motifs**
- 4-mer and 5-mer frequencies at fragment termini
- Strand-specific end motif patterns
- End motif diversity (entropy)
- Deviation from healthy cfDNA end motif profiles
- Nuclease-specific signatures (DNASE1L3 cleavage patterns)

#### 3. **Nucleosome Spacing (FFT-Based)**
- Fast Fourier Transform analysis of coverage profiles
- Detection of ~167 bp nucleosomal periodicity
- Amplitude and phase measurements
- Comparison to normal tissue nucleosome positioning
- Aberrant spacing indicative of disrupted tumor chromatin

#### 4. **Copy Number Alterations**
- Arm-level z-scores and segmental CNAs
- CNA burden (fraction of genome altered)
- Number and amplitude of breakpoints
- Focal vs. broad alteration patterns
- Integration with known pediatric tumor CNA signatures

#### 5. **Coverage-Based Features**
- GC-corrected and mappability-corrected coverage
- Bin-to-bin variance and regional anomalies
- Coverage correlation with reference samples

### Statistical Modeling Strategy

We follow a **progressive complexity approach**:

```
Phase 1: Classical Statistics (Baseline)
├── Linear regression with informative features
├── Beta regression (TF naturally bounded 0-1)
└── Generalized Additive Models (GAMs) for non-linear relationships

Phase 2: Tree-Based Methods
├── Gradient Boosting (XGBoost/LightGBM) with SHAP interpretability
└── Quantile Regression Forests (direct uncertainty estimation)

Phase 3: Deep Learning (if needed)
└── Transformer architecture (only for final performance gains)
```

**Philosophy**: Use the simplest model that achieves target performance. Only increase complexity when necessary and quantifiable gains are demonstrated.

### Uncertainty Quantification

cfTF provides **multiple levels of uncertainty**:

1. **Aleatoric uncertainty**: Inherent data noise (from quantile regression or neural net variance)
2. **Epistemic uncertainty**: Model uncertainty (from bootstrap ensembles or Monte Carlo dropout)
3. **Feature-based reliability**: Quality scores based on input characteristics
4. **Conformal prediction intervals**: Distribution-free 95% confidence intervals
5. **Sample-specific flags**: Warnings when predictions may be unreliable

---

## Use Cases

### Primary Applications

1. **Treatment Response Monitoring**
   - Track tumor burden during therapy
   - Detect early progression or resistance
   - Guide treatment decisions

2. **Minimal Residual Disease (MRD) Detection**
   - Post-treatment surveillance
   - Early relapse detection
   - Risk stratification

3. **Longitudinal Disease Tracking**
   - Joint modeling of multiple timepoints
   - Trajectory analysis (response, plateau, progression)
   - Change point detection

### Optional Enhanced Modes

When additional data is available, cfTF can leverage:

- **Tumor-informed mode**: Use matched tumor biopsy CNAs for improved accuracy
- **Matched normal mode**: Patient-specific baseline from normal tissue
- **Tissue-of-origin informed**: Cancer type-specific feature weighting
- **Temporal modeling**: Borrow strength across timepoints for improved low-TF detection

---

## Technical Specifications

### Input Requirements

**Minimum:**
- BAM/CRAM file from low-pass WGS of cfDNA
- Sequencing coverage: 0.5-5× (optimized for ~3×)
- Reference genome: hg38 (hg19 supported)

**Optional:**
- Matched tumor biopsy (BAM/CRAM or VCF with CNAs)
- Matched normal tissue (BAM/CRAM)
- Patient metadata: age, cancer type, treatment timeline
- Clinical information: tissue of origin, disease stage

### Output

```
tumor_fraction_estimate: 0.087
confidence_interval_95: [0.052, 0.134]
quality_score: 0.89
prediction_reliable: true

uncertainty_components:
  aleatoric: 0.021
  epistemic: 0.015
  
feature_contributions:
  copy_number_signal: 0.45
  fragment_size_signal: 0.32
  end_motif_signal: 0.15
  nucleosome_signal: 0.08

quality_control:
  sequencing_pass: true
  coverage_mean: 3.2
  fragment_size_peak: 167
  model_confidence: high
  warnings: []
```

### Performance Targets

- **Sensitivity at 3% TF**: >95%
- **Sensitivity at 1% TF**: >80% (stretch goal)
- **Specificity**: >98% (healthy donors)
- **Uncertainty calibration**: 95% CI should contain true TF 95% of the time
- **Runtime**: <10 minutes per sample on standard hardware
- **Reproducibility**: CV <5% on technical replicates

---

## Design Principles

### 1. **Interpretability Over Black-Box Performance**
Every prediction should be explainable. Feature importance, contribution plots, and decision pathways must be transparent.

### 2. **Uncertainty is Not Optional**
Never provide a point estimate without confidence intervals. Flag predictions when uncertainty is high.

### 3. **Quality Control First**
Catch bad samples before they produce misleading results. Comprehensive QC metrics at every step.

### 4. **Clinical Validation Standards**
Follow CLIA/CAP guidelines for analytical validation. Document performance characteristics rigorously.

### 5. **Pediatric-Specific**
Age-appropriate normal references. Pediatric tumor biology (neuroblastoma, Wilms, medulloblastoma, etc.) built into model.

### 6. **Open Science**
Open-source code, pre-trained models, documented methodologies. Enable reproducibility and community improvement.

---

## Comparison to Existing Tools

| Feature | ichorCNA | DELFI-TF | cfTF (this tool) |
|---------|----------|----------|------------------|
| **Coverage required** | 0.1× | 6× | 0.5-5× (optimal: 3×) |
| **Requires CNAs** | Yes | No | No (but uses when present) |
| **Fragmentomic features** | No | Yes | Yes (multiple types) |
| **Confidence intervals** | No | No | Yes (rigorous) |
| **Open source** | Yes | No | Yes |
| **Interpretability** | Moderate | Low | High |
| **Clinical QC metrics** | Minimal | Unknown | Comprehensive |
| **Pediatric-optimized** | No | No | Yes |
| **Longitudinal modeling** | No | No | Yes |
| **Tumor-informed mode** | No | No | Yes |

---

## Project Status

🚧 **Under Active Development** 🚧

### Roadmap

- [x] Project design and requirements gathering
- [ ] **Phase 1**: Foundation (Months 1-2)
  - [ ] Data pipeline (BAM → features)
  - [ ] Feature extraction modules
  - [ ] QC metrics calculator
  - [ ] Ground truth dataset curation
  
- [ ] **Phase 2**: Classical Models (Months 3-4)
  - [ ] GLM/GAM implementation
  - [ ] Uncertainty quantification
  - [ ] Validation framework
  - [ ] Benchmark vs. ichorCNA
  
- [ ] **Phase 3**: ML Enhancement (Months 5-6)
  - [ ] Gradient boosting + SHAP
  - [ ] Quantile regression forests
  - [ ] Ensemble methods
  
- [ ] **Phase 4**: Advanced Features (Months 7-8)
  - [ ] Longitudinal modeling
  - [ ] Tumor-informed mode
  - [ ] Transformer architecture (if needed)
  
- [ ] **Phase 5**: Clinical Validation (Months 9-12)
  - [ ] Prospective validation
  - [ ] Clinical QC optimization
  - [ ] Documentation for clinical use

---

## Technology Stack

### Core Libraries
```python
# Scientific computing
numpy, scipy, pandas

# Statistical modeling
statsmodels      # GLMs, GAMs
pymc            # Bayesian inference
scikit-learn    # Baseline ML

# Advanced ML
xgboost, lightgbm    # Gradient boosting
skranger             # Quantile regression forests
torch                # Transformers (if needed)

# Genomics
pysam           # BAM/CRAM I/O
pyBigWig        # Coverage tracks
scikit-bio      # Bioinformatics utilities

# Signal processing
scipy.fft       # Nucleosome spacing (FFT)
scipy.signal    # Filtering, smoothing

# Uncertainty quantification
mapie           # Conformal prediction

# Interpretability
shap            # Feature importance
lime            # Local explanations

# Performance
numba           # JIT compilation
dask            # Parallel processing
```

---

## Installation

```bash
# Clone repository
git clone https://github.com/yourusername/cfTF.git
cd cfTF

# Create conda environment
conda env create -f environment.yml
conda activate cfTF

# Install package
pip install -e .

# Run tests
pytest tests/
```

---

## Quick Start

```python
from cfTF import TumorFractionEstimator

# Initialize model
estimator = TumorFractionEstimator(
    coverage_threshold=0.5,
    model_type='gradient_boosting',
    uncertainty_method='conformal'
)

# Load sample
sample = estimator.load_sample('sample.bam', reference='hg38')

# Estimate tumor fraction
result = estimator.predict(sample)

print(f"Tumor Fraction: {result.TF:.3f}")
print(f"95% CI: [{result.CI_lower:.3f}, {result.CI_upper:.3f}]")
print(f"Quality Score: {result.quality_score:.2f}")
print(f"Reliable: {result.reliable}")

# Get detailed QC report
qc_report = result.quality_control
print(qc_report)

# Visualize feature contributions
result.plot_feature_importance()
result.plot_uncertainty_breakdown()
```

---

## Contributing

We welcome contributions! Areas where help is particularly valuable:

- **Validation datasets**: Samples with orthogonal TF measurements (ddPCR, targeted sequencing)
- **Pediatric tumor profiles**: Age-stratified normal references, tumor-specific CNA signatures
- **Algorithm improvements**: Better feature engineering, model architectures
- **Clinical testing**: Real-world validation in clinical cohorts
- **Documentation**: Tutorials, use cases, best practices

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

---

## Citation

*Manuscript in preparation*

For now, please cite this repository:

```bibtex
@software{cfTF2026,
  title = {cfTF: Cell-Free Tumor Fraction Estimation},
  author = {Your Name},
  year = {2026},
  url = {https://github.com/yourusername/cfTF}
}
```

---

## License

MIT License - see [LICENSE](LICENSE) file for details.

---

## Acknowledgments

This project builds on pioneering work from:

- **ichorCNA**: Adalsteinsson et al., Nature Communications 2017
- **DELFI**: Cristiano et al., Nature 2019; van 't Erve et al., Nature Communications 2024
- The broader liquid biopsy and cfDNA fragmentomics research community

Special thanks to the pediatric oncology community for their insights into clinical needs and challenges.

---

## Contact

- **Issues**: [GitHub Issues](https://github.com/yourusername/cfTF/issues)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/cfTF/discussions)
- **Email**: your.email@institution.edu

---

## Frequently Asked Questions

**Q: Why Python instead of R?**  
A: Better ecosystem for modern ML, easier GPU integration if needed, broader clinical software compatibility, and more accessible for software engineering best practices.

**Q: Will this work for my cancer type?**  
A: cfTF is designed to be pan-cancer, but is optimized for pediatric cancers. Performance may vary by tumor type depending on CNA burden and cfDNA shedding patterns.

**Q: What if I only have 0.3× coverage like ichorCNA?**  
A: The tool will work but with reduced performance on fragmentomic features. For best results, we recommend 1-5× coverage.

**Q: Can I use this for early detection/screening?**  
A: cfTF is designed for tumor fraction estimation in known cancer patients. Early detection requires different models and thresholds. However, the fragmentomic features could be adapted for screening applications.

**Q: How does this handle clonal hematopoiesis (CH)?**  
A: CH can confound mutation-based TF estimates but has minimal impact on fragmentomic features. Our multi-modal approach is more robust to CH than mutation-based methods.

---

**Last Updated**: February 2026
