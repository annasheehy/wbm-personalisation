# wbm-personalisation

MATLAB code for testing how physiological data, blood metabolomics, fasting and physical activity change the predictions of the Harvey and Harvetta whole-body metabolic models (WBMs).

- varies physiological parameters (body weight, heart rate, haematocrit, glomerular filtration rate, height, urine creatinine) across the ranges tested and measures their effect on predicted fluxes and basal metabolic rate (BMR);
- personalises the models with synthetic blood profiles of 26 metabolites, tests feasibility, finds why models fail and tests a fix;
- builds fasted and active versions of the models and compares their predictions with fluxes measured in humans.

## Requirements

- MATLAB (R2020b was used)
- [COBRA Toolbox](https://github.com/opencobra/cobratoolbox), which provides Harvey and Harvetta 1.03d (`loadPSCMfile`) and the PSCM toolbox functions
- Gurobi (13.0.0 was used) for LP and QP

The COBRA Toolbox functions `persWBM`, `persWBMmetabolomics`, `physiologicalConstraintsHMDBbased` and `calcOrganFract` are called unchanged. Two unit errors found in this work are fixed in the COBRA Toolbox (branch `ch5-bugfixes`, to be merged into `develop`): `individualizedLabReport` cardiac output options 3 to 5 returned L/day where mL/min is used, and `persWBMmetabolomics` read urine concentrations in µmol/L as µmol per mmol creatinine. `heightCardiacOutput.m` uses the corrected option 4; the chapter uses blood metabolites only.

## Setup

```matlab
initCobraToolbox(false)
setupPersWBM('path/to/project')   % adds the code folders to the path and sets the project folder
```

The project folder holds the inputs and receives the outputs:

| Folder | Contents |
|---|---|
| `references/` | `healthyReferenceValues.csv` (blood metabolite values from the published studies) and `healthyPhysiologyRanges.csv` (ranges of the physiological inputs), provided in the Supplementary Material |
| `results/` | checkpoints, CSV files and figures written by the scripts |
| `ThesisChapter/LateX` | the chapter (`main.tex`), its appendix (`appendix_chapter5.tex`, written by `supplementaryMaterial.m`) and `Figures` |
| `results/latex` | LaTeX tables written by the scripts, to be pasted into `main.tex` |

## Code

| Folder | Contents |
|---|---|
| `common/` | shared functions: BMR conditions (`setBMRConditions`), BMR from a flux distribution (`bmrFromFlux`, `calculateModelBMR`), flux ranges under BMR conditions (`fastingTargetFVA`), fed models rebuilt with other parameters (`buildFedWBM`, `onLiftedModel`), dependent physiological parameters (`recalcDependentParams`), body surface area (`bsaDuBois`), net organ exchanges (`netExchange`), the 11 read-outs (`oatReadouts`), paths (`persWBMpaths`) |
| `fasting/` | control (fed) and fasted models, the tests that set the fasted state, comparison with fluxes measured in fasting humans, urine creatinine |
| `metabolomics/` | personalisation with synthetic blood metabolite profiles: bounds changed, feasibility, crossed bounds |
| `inputEffects/` | effect of each physiological parameter and metabolite on its own; height through cardiac output and blood volume |
| `active/` | one hour of running |
| `supplement/` | supplementary tables (LaTeX appendix) and Excel workbooks of the inputs and results |

## Run order

All scripts start from the models built by `fastingWBM.m`.

| Step | Script | Output in the chapter |
|---|---|---|
| 1 | `fasting/fastingWBM.m` | control and fasted models, fed and fasted BMR |
| 2 | `fasting/fastingSetupTests.m` | brain and heart ATP supply of the fasted models for each amino-acid release limit and with phosphate uptake capped |
| 3 | `fasting/fastingTables.m` | tables of fasted BMR and gas exchange, and of the comparison with fasting humans (computes the FVA of the benchmarks on its first run) |
| 4 | `fasting/fastingUrineCreatinine.m` | fasted models with urine creatinine at the low and high end of its range: feasibility, BMR and urinary urea excretion |
| 5 | `metabolomics/metabolomicsPersonalisation.m` | metabolite panel table; figures of bounds per metabolite, uptake limits and feasibility |
| 6 | `inputEffects/metabolomicsOAT.m` | table of ranges tested; figure of effects on organ exchanges. Run again with `oatMode = 'l2'` for the figure of effects on the basal flux distribution and BMR |
| 7 | `inputEffects/oatBMRdecomposition.m` | BMR change with heart rate and haematocrit, by reaction and organ |
| 8 | `inputEffects/heightCardiacOutput.m` | models with cardiac output from height and personalised blood volume |
| 9 | `active/activeWBM.m` | table of extra muscle work and the run |
| 10 | `supplement/supplementaryMaterial.m` | appendix Tables S1-S9 (`appendix_chapter5.tex`) and the Excel workbooks of inputs and results, from the saved results of steps 1-9 |

`metabolomicsOAT.m` and `heightCardiacOutput.m` save one file per case and skip finished cases, so they can be run in several MATLAB sessions at once (see the header of each script).

## Author

Anna Sheehy, University of Galway
