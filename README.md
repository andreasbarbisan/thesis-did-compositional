# Difference-in-Differences in Educational Policies: Evidence on Full-Day Education under Compositional Changes

**Andreas Azambuja Barbisan**  
Insper Institute of Education and Research, São Paulo, Brazil  
andreas.barbisan@gmail.com

---

## Overview

This repository contains the replication code and paper source for my undergraduate thesis (Insper, 2026), which evaluates proficiency gains associated with the Full-Time Education Program (*Programa Ensino Integral*, PEI) in São Paulo state schools.

The central empirical challenge is that PEI is adopted in a staggered fashion and is documented to change the composition of the student body after conversion. Standard DiD estimators for repeated cross-section data rely on a stationarity assumption that may be strained when treatment itself alters who is observed in the sample. The paper uses this as a laboratory to compare four estimators of increasing rigor:

| Estimator | What it addresses |
|---|---|
| TWFE | Baseline benchmark |
| Callaway & Sant'Anna (2021) | Staggered adoption and treatment effect heterogeneity |
| Sant'Anna & Zhao (2020) | Observable covariates via doubly robust DiD |
| Aggregate S&X approximation | Sensitivity to observable compositional change at school level |

Main estimates indicate average post-treatment gains of **0.50–0.68 standard deviations** in standardized SARESP scores (Portuguese Language and Mathematics), for 9th grade and 3rd grade of high school. The compositional correction does not eliminate the results, though data limitations prevent a full individual-level decomposition.

---

## Papers

- **Working paper** (`paper/`): English-language working paper adapted from the thesis — *Difference-in-Differences in Educational Policies: Evidence on Full-Day Education under Compositional Changes* (May 2026).
- **Thesis** (`thesis/`): Full undergraduate thesis in Portuguese, submitted to Insper (2026).

Both compiled PDFs (`main.pdf`) are included in their respective folders.

---

## Data

Raw data are **not included** in this repository. The analysis uses three publicly available sources:

| Source | Description | Access |
|---|---|---|
| SARESP | Aggregate school-level proficiency by year, grade, and subject | [SP Education Open Data](https://dados.educacao.sp.gov.br/story/saresp) |
| Enrollment microdata | Student-level administrative records from São Paulo state schools | [SP Education Open Data](https://dados.educacao.sp.gov.br/) |
| IBGE municipal data | Municipal GDP per capita, population, and age structure | [IBGE](https://www.ibge.gov.br/estatisticas/economicas/contas-nacionais/9088-produto-interno-bruto-dos-municipios.html) |

The PEI first-treatment year by school (`escolas_pei_first_treat.dta`) was constructed from SARESP school records and is available upon request.

The main estimation window is **2011–2018**. The unit of analysis is the school × year × grade × subject panel.

---

## Repository Structure

```
thesis-did-compositional/
├── README.md
├── code/
│   ├── 00_prep_controles.do          # Aggregate enrollment microdata → school × year × grade covariates
│   ├── 01_prep.do                    # Stack SARESP aggregates, merge first_treat, create DiD variables
│   ├── 01b_merge_controles.do        # Merge composition covariates; impute missing values; save panel
│   ├── 02_estimacao.do               # Main estimation: TWFE, CS21, DR-DiD, S&X approximation
│   ├── 03_tabelas_resultados.do      # Export main results tables to LaTeX
│   ├── 04_densidade_trat_controle.do # Density diagnostic: proficiency distribution post-treatment
│   ├── 04b_densidade_pre_english.do  # English-language version of density figures (for working paper)
│   ├── 05_diagnosticos_revisao.do    # Diagnostics: mean reversion, spillovers, imputation sensitivity
│   ├── 06_apendices_estimacoes.do    # Appendix tables: ATT by relative period for all specifications
│   ├── 06b_export_appendix_figures_english.do  # English-language appendix figures (for working paper)
│   ├── 07_participacao_saresp.do     # Diagnostic: SARESP participation rate relative to enrollment
│   ├── 08_hausman_mde.do             # Hausman-type test (DR-DiD vs S&X) and minimum detectable effect
│   └── exploratory/                  # Auxiliary and diagnostic scripts (not part of main pipeline)
│       ├── 00b_prep_endereco_controles_municipais.do
│       ├── 00c_inspeciona_controles_municipais.do
│       ├── 00d_diagnostico_covariadas_estimacao.do
│       ├── 00e_testa_csdid_controles_municipais.do
│       ├── 02b_sx_nivel_b.do
│       └── 02c_hausman_patch.do
├── paper/                            # Working paper (English) — LaTeX source + compiled PDF
│   ├── main.tex
│   ├── main.pdf
│   ├── chapters/
│   ├── appendix/
│   ├── files/                        # Tables (.tex) and figures (.pdf/.png)
│   └── referencias.bib
└── thesis/                           # Undergraduate thesis (Portuguese) — LaTeX source + compiled PDF
    ├── main.tex
    ├── main.pdf
    ├── 2-textuais/                   # Chapter source files
    ├── tabelas/                      # Tables and figures
    └── referencias.bib
```

---

## Pipeline

Scripts are numbered and should be run in order. Each script reads from and writes to a `checkpoints/` directory (not included — created automatically on first run).

```
00_prep_controles      →  checkpoints/controles_escola_ano_serie.dta
01_prep                →  checkpoints/painel_escola_ano.dta
01b_merge_controles    →  checkpoints/painel_completo.dta
02_estimacao           →  checkpoints/estimacao_resumo_tese.csv
                          checkpoints/estimacao_eventstudy_long.csv
                          checkpoints/estimacao_sx_detalhe.csv
03_tabelas_resultados  →  checkpoints/*.tex  (LaTeX tables)
04_densidade_*         →  figuras/*.pdf
05_diagnosticos_*      →  checkpoints/diagnostico_*.dta
06_apendices_*         →  figuras/appendix/*.pdf
07_participacao_saresp →  checkpoints/participacao_saresp.csv
08_hausman_mde         →  checkpoints/tabela_hausman_composicao.tex
```

---

## Requirements

**Stata 16 or later** is required. The following user-written packages must be installed before running the pipeline:

```stata
ssc install csdid        // Callaway & Sant'Anna (2021) estimator
ssc install drdid        // Doubly robust DiD (Sant'Anna & Zhao 2020)
ssc install reghdfe      // High-dimensional fixed effects regression
ssc install coefplot     // Coefficient plots and event study figures
ssc install ftools       // Fast Stata tools (dependency of reghdfe)
```

---

## Key Results

Average post-treatment effect on standardized SARESP proficiency (mean of k=0 to k=4 relative to conversion):

| Grade–Subject | CS21 | DR-DiD | S&X approx. | TWFE |
|---|---|---|---|---|
| 9th grade – Portuguese | 0.491 | 0.510 | 0.502 | 0.475 |
| 9th grade – Mathematics | 0.585 | 0.601 | 0.603 | 0.583 |
| HS 3rd grade – Portuguese | 0.526 | 0.571 | 0.601 | 0.507 |
| HS 3rd grade – Mathematics | 0.605 | 0.627 | 0.675 | 0.604 |

A Hausman-type diagnostic comparing DR-DiD and the S&X approximation yields high p-values across all cells, but the minimum detectable effect at 80% power ranges from 0.20 to 0.25 standard deviations. The test is best read as a sensitivity diagnostic, not as proof of compositional robustness.

---

## Notes on the S&X Approximation

The aggregate approximation inspired by Sant'Anna & Xu (2026) does **not** fully reproduce the individual-level estimator. A literal implementation would require linking individual enrollment records to individual SARESP scores — a merge that is not feasible with the publicly available data due to inconsistent student identifiers across databases.

The adaptation operates at the school × year × grade level: for each cohort and relative time, a four-category multinomial logit is estimated using baseline composition covariates, Hájek weights are constructed from the resulting propensity scores, and three outcome regressions adjust each counterfactual cell. Cells with insufficient common support are excluded rather than replaced by a simpler estimator.

---

## Citation

> Barbisan, A. A. (2026). *Difference-in-Differences in Educational Policies: Evidence on Full-Day Education under Compositional Changes*. Undergraduate Thesis, Insper Institute of Education and Research.
