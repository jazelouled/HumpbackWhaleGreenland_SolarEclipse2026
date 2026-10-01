# Humpback whale solar-eclipse analysis

<p align="center">
  <strong>Reproducible analysis pipeline for a satellite-tagged humpback whale observed during the 12 August 2026 solar eclipse in East Greenland.</strong>
</p>

<p align="center">
  <a href="#quick-start">Quick start</a> ·
  <a href="#analysis-overview">Analysis overview</a> ·
  <a href="#outputs">Outputs</a> ·
  <a href="#reproducibility-notes">Reproducibility</a>
</p>

> [!NOTE]
> This repository accompanies the manuscript **“Shorter-than-expected dives by a humpback whale during a solar eclipse in East Greenland.”** One self-contained R script regenerates the processed datasets, fitted models, manuscript Figure 4, supplementary figures, Tables S1--S8, audit files, and a session-information record.

<p align="center">
  <img src="https://img.shields.io/badge/language-R-276DC3?logo=r&logoColor=white" alt="R" />
  <img src="https://img.shields.io/badge/analysis-reproducible-16803c" alt="Reproducible analysis" />
  <img src="https://img.shields.io/badge/timezone-UTC-4B5563" alt="Timezone UTC" />
</p>

---

## What this pipeline does

The analysis uses transmitted dive records and Argos locations from a Wildlife Computers SPLASH10-373 tag deployed on humpback whale ID 260680 in East Greenland. It tests whether dives recorded during the eclipse differed from comparable non-eclipse dives while accounting for diel solar conditions, dive depth, and temporal dependence.

| Component | Question addressed | Main output |
| --- | --- | --- |
| **Data processing** | How are dive duration, maximum depth, daily position and solar elevation derived? | Table S1; processed dive record |
| **Principal GAMMs** | Did eclipse dives differ in duration or depth? | Tables S2--S3; Figure 4C |
| **Sensitivity analyses** | Does the duration result depend on reasonable analytical choices? | Table S4 |
| **Counterfactual predictions** | How long would the 13 eclipse dives be expected to last without the eclipse effect? | Table S5 |
| **125-min episode analysis** | How unusual was the complete eclipse period? | Table S6; Figure 4D |
| **Depth-sequence analysis** | Was there an unusual shallow interruption within the eclipse sequence? | Table S7; Figure S2 |
| **Obscuration analysis** | Did behaviour vary continuously with mean solar obscuration? | Table S8 |

---

## Repository layout

```text
.
├── eclipse_final_analysis.R          # Complete analysis and figure-production script
├── 00input/
│   └── 260680/                       # Raw Wildlife Computers exports (not distributed)
│       ├── 260680-Behavior.csv
│       └── 260680-Locations.csv
├── outputs_final/                    # Created when the pipeline is run
│   ├── data/                         # Processed analysis datasets
│   ├── tables/                       # Tables S1--S8 as CSV files
│   ├── figures/                      # Figure 4 and Figures S1--S2
│   ├── models/                       # Fitted model objects (.rds)
│   ├── audit/                        # Supporting checks and diagnostics
│   ├── FINAL_ANALYSIS_REPORT.txt     # Compact results summary
│   └── sessionInfo.txt               # Exact R-session information
└── supplementary_latex/              # LaTeX source for the supplementary material
    ├── supplementary.tex
    ├── figures/
    └── README.md
```

> [!IMPORTANT]
> Raw tag data are not included in the repository. To reproduce the analysis, place the two original Wildlife Computers CSV exports in `00input/260680/` using the filenames shown above.

---

## Quick start

### 1. Install R packages

The script checks for the required packages and installs any that are missing:

```r
required_packages <- c(
  "mgcv", "nlme", "suncalc", "dplyr", "ggplot2", "grid",
  "gridExtra", "svglite", "scales"
)
```

### 2. Set the project path

At the start of [`eclipse_final_analysis.R`](eclipse_final_analysis.R), set `project_dir` to the local path of this repository:

```r
project_dir <- "/path/to/your/local/repository"
```

### 3. Run the full pipeline

```r
source("eclipse_final_analysis.R")
```

The submitted analysis uses a moving-block bootstrap with 10,000 replicates. For a quick development run only, set:

```r
run_block_bootstrap <- FALSE
```

> [!TIP]
> A full run writes to a new `outputs_final/` directory. Existing files with the same names are replaced, so rename the output directory first if you want to retain an earlier run.

---

## Analysis overview

### 1. Data processing and solar elevation

Dive duration is calculated from the transmitted start and end times. Maximum dive depth is the midpoint of the transmitted depth interval. Each dive is associated with the median daily latitude and longitude of Argos quality classes 1--3. Solar elevation is calculated at each dive midpoint using `suncalc::getSunlightPosition` in UTC.

### 2. Principal analysis

The principal comparison includes dives occurring within the exact solar-elevation range observed during the eclipse (25.57--35.06°), excluding the first and last incomplete recording days. It contains 151 dives: 13 eclipse dives and 138 non-eclipse dives.

Separate generalized additive mixed models (GAMMs) are fitted for duration and maximum depth:

\[
D_i = \beta_0 + \beta_E E_i + f_{solar}(S_i) + f_{depth}(Z_i) + \varepsilon_i
\]

\[
Z_i = \gamma_0 + \gamma_E E_i + g_{solar}(S_i) + \eta_i
\]

where \(E_i\) denotes eclipse status, \(S_i\) solar elevation and \(Z_i\) maximum dive depth. Smooth terms use thin-plate regression splines with \(k = 8\).

Candidate temporal residual structures are compared by maximum-likelihood AIC: no temporal correlation, continuous-time AR(1), discrete AR(1), discrete AR(2), and discrete ARMA(1,1). The selected structure is refitted by REML for inference:

- **Dive duration:** discrete ARMA(1,1)
- **Maximum dive depth:** discrete AR(2)

### 3. Complementary analyses

<details>
<summary><strong>Open details</strong></summary>

#### Sensitivity analyses

The duration model is rerun while varying eclipse definition, depth representation, solar-elevation control range, partial-day inclusion, residual-correlation structure, spline flexibility, and short-term dependence terms.

#### Counterfactual predictions

For each eclipse dive, observed solar elevation and depth are retained while eclipse status is set to control. These predictions quantify the model-based, adjusted reduction in duration.

#### Eclipse-episode windows

A reference duration model without eclipse status generates observed-minus-expected residuals. The mean residual over the 125-min eclipse is compared descriptively with the same statistic calculated in non-eclipse 125-min windows.

#### Depth sequences

Each 125-min sequence is assessed for a shallow interior excursion relative to diving at both its beginning and end. Overlapping control sequences are descriptive comparisons, not independent replicates.

#### Obscuration

Mean eclipse obscuration is calculated at 30-s intervals over each dive and replaces the binary eclipse term in exploratory duration and depth models.

</details>

---

## Outputs

| Directory | Contents |
| --- | --- |
| `outputs_final/data/` | Complete processed dive record, principal-analysis dataset, window-analysis dataset and obscuration-analysis dataset |
| `outputs_final/tables/` | `Table_S1` through `Table_S8` in CSV format |
| `outputs_final/figures/` | Manuscript Figure 4 (`.png`, `.pdf`, `.svg`) and supplementary Figures S1--S2 (`.png`, `.pdf`) |
| `outputs_final/models/` | All principal, baseline, window and obscuration model objects in `supplementary_analysis_objects.rds` |
| `outputs_final/audit/` | Input-resolution summaries, eclipse-definition checks, common-support diagnostics, nearest-control checks and moving-block bootstrap results |

### Principal reported result

> [!SUCCESS]
> After adjustment for solar elevation, dive depth and temporal dependence, eclipse dives were estimated to be **2.89 min shorter** than comparable non-eclipse dives (95% CI: -5.05 to -0.74 min; \(P = 0.009\)).

---

## Reproducibility notes

- **Times:** all timestamps are parsed and analysed in UTC.
- **Position:** one daily representative position is calculated from Argos quality classes 1--3.
- **Model fitting:** candidate residual structures are compared by maximum likelihood; selected models are refitted by REML.
- **Randomness:** the moving-block bootstrap fixes a deterministic seed for each tested block size (`20260822 + block_n`).
- **Software record:** the exact package and R-session versions are saved in `outputs_final/sessionInfo.txt`.
- **Traceability:** every reported supplementary table is exported as a standalone CSV, and the fitted objects used to produce them are retained as an `.rds` file.

---

## Supplementary material

The submission-ready supplementary material is in [`supplementary_latex/`](supplementary_latex/). It includes the supplementary methods, Tables S1--S8 and Figures S1--S2.

To compile it locally:

```bash
cd supplementary_latex
latexmk -pdf supplementary.tex
```

---

## Citation

If you use or adapt this workflow, please cite the associated manuscript:

> Ouled-Cheikh, J., Ugarte, F., Huffeldt, N. P., Heide-Jørgensen, M. P., and colleagues. *Shorter-than-expected dives by a humpback whale during a solar eclipse in East Greenland.*

---

<p align="center">
  <sub>Analysis code for the East Greenland humpback whale solar-eclipse study.</sub>
</p>
