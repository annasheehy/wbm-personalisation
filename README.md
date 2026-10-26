# wbm-personalisation

MATLAB code for testing how physiological data, blood metabolomics, fasting and physical activity change the predictions of the Harvey and Harvetta whole-body metabolic models (WBMs).

This is the code for Chapter 5 of the PhD thesis of Anna Sheehy (University of Galway), *Testing the personalisation of whole-body metabolic models: physiology, blood metabolomics, fasting and physical activity* (Anna Sheehy and Ines Thiele). The chapter:

- varies physiological parameters (body weight, heart rate, haematocrit, glomerular filtration rate, height) across the ranges tested and measures their effect on predicted fluxes and basal metabolic rate (BMR);
- personalises the models with synthetic blood profiles of 26 metabolites, tests feasibility, finds why models fail and tests a fix;
- builds fasted and active versions of the models and compares their predictions with fluxes measured in humans.

## Requirements

- MATLAB (R2020b was used)
- [COBRA Toolbox](https://github.com/opencobra/cobratoolbox) (commit 91c8cbf27 was used), which provides Harvey and Harvetta 1.03d (`loadPSCMfile`) and the PSCM toolbox functions
- Gurobi (13.0.0 was used) for LP and QP

The COBRA Toolbox functions `persWBM`, `persWBMmetabolomics`, `physiologicalConstraintsHMDBbased` and `calcOrganFract` are called unchanged.

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
| `ThesisChapter/LateX/Tables`, `Figures` | LaTeX tables and figures of the chapter |

## Code

| Folder | Contents |
|---|---|
| `common/` | BMR from the basal flux distribution (`calculateModelBMR`), flux ranges under BMR conditions (`fastingTargetFVA`), dependent physiological parameters (`recalcDependentParams`), paths (`persWBMpaths`) |
| `fasting/` | control (fed) and fasted models, comparison with fluxes measured in fasting humans, tests of adipose glycerol release |
| `metabolomics/` | personalisation with synthetic blood metabolite profiles: bounds changed, feasibility, crossed bounds |
| `inputEffects/` | effect of each physiological parameter and metabolite on its own; height through cardiac output and blood volume |
| `active/` | one hour of running |

## Run order

All scripts start from the models built by `fastingWBM.m`.

| Step | Script | Output in the chapter |
|---|---|---|
| 1 | `fasting/fastingWBM.m` | control and fasted models, fasted BMR |
| 2 | `fasting/fastingFigures.m` | flux ranges of the fasting comparison (used by step 3) |
| 3 | `fasting/fastingTables.m` | tables of fasted BMR and gas exchange, and of the comparison with fasting humans |
| 4 | `metabolomics/metabolomicsPersonalisation.m` | metabolite panel table; figures of bounds per metabolite, uptake limits and feasibility |
| 5 | `inputEffects/metabolomicsOAT.m` | table of ranges tested; figure of effects on organ exchanges. Run again with `oatMode = 'l2'` for the figure of effects on the basal flux distribution and BMR |
| 6 | `inputEffects/oatBMRdecomposition.m` | BMR change with heart rate and haematocrit, by reaction and organ |
| 7 | `inputEffects/heightCardiacOutput.m` | models with cardiac output from height and personalised blood volume |
| 8 | `active/activeWBM.m` | table of extra muscle work and the run |
| 9 | `fasting/fastingAdiposeFix.m` | tests of adipose glycerol release in the fasted models |

`fasting/fastingWBM_investigations.m` and `fasting/femaleBrainATP_investigation.m` contain the tests of amino acid release and phosphate uptake that set the fasted state.

`metabolomicsOAT.m` and `heightCardiacOutput.m` save one file per case and skip finished cases, so they can be run in several MATLAB sessions at once (see the header of each script).

## Author

Anna Sheehy, University of Galway
