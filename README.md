# Replication and Extension of Hurst and Sintov (2022, Study 2)

This repository contains the materials, the power simulation and the analysis code of a preregistered very close replication and extension of Hurst and Sintov (2022, Study 2). The study examines how negative feedback about the carbon footprint of a chosen company affects guilt, anticipated emotions and the decision to donate part of a possible bonus to an environmental organization.

**Status:** in preparation. The preregistration and the data will be added here.

## Contents

| Folder | What it contains |
|---|---|
| `materials/` | The full questionnaire, the informed consent form and the debriefing |
| `pilot/` | The report of the pilot study |
| `analysis/power/` | The power simulation (`.md`) and the simulated data (`.rds`)|
| `analysis/preregistered/` | The preregistered analysis, tested on one simulated data set (`.md`) |
| `data/` | The anonymized data (added after data collection) |

## Viewing the Results

The files (`.md`) contain the text, the code, the tables and the figures.

## Running the Code

The analyses run in R.

1. Install the R packages `effectsize`, `MASS`, `ggplot2`, `patchwork` and `knitr`.
2. Download PROCESS for R (Hayes) from https://www.processmacro.org and save it as `analysis/PROCESS/process.R`. PROCESS is not included here, because its license does not allow redistribution.
3. Open a `.md` file in RStudio and click **Render**.

The power simulation reuses its saved results (`power_results_simulated.rds`) and renders in a few seconds. Without the saved results, it runs the full simulation, which takes about 75 minutes.

## Original Study

Hurst, K. F., & Sintov, N. D. (2022). Guilt consistently motivates pro-environmental outcomes while pride depends on context. *Journal of Environmental Psychology, 80*, 101776. https://doi.org/10.1016/j.jenvp.2022.101776

## License

The materials in this repository are licensed under the Creative Commons Attribution-ShareAlike 4.0 International License (CC BY-SA 4.0).

## Contact

Ewa Duda, Department of Sociology, Faculty of Social Sciences, Charles University, Prague (54855986@fsv.cuni.cz)
