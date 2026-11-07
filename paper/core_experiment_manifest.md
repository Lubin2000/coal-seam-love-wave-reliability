# Paper 1 core experiment manifest

## Current manuscript

`manuscript_en_v0.3.md`

## Execution order

Run from any directory:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_core_experiments.ps1
```

Use `-SkipWavefields` to regenerate the analytical, tracking, inversion,
calibration, bootstrap and QEDispInv results without rerunning the expensive
finite-difference wavefields. Use `-FiguresOnly` to redraw all publication
figures from the existing summaries.

## Core experiments

| Script | Purpose | Main outputs |
|---|---|---|
| `analysis/e13_acquisition_information.py` | Fisher-information budget with fixed physical mode order | E13 CSV and JSON |
| `analysis/e14_mode_availability.py` | Random-model information under partial mode availability | E14 CSV and JSON |
| `analysis/e15_automatic_mode_picking.py` | Independent peak baseline | E15 CSV and JSON |
| `analysis/e16_robust_mode_tracking.py` | Greedy and fixed-slot global tracking | E16 CSV and JSON |
| `analysis/e17_robust_likelihood.py` | Gaussian, Huber and Student-t inversion | E17 CSV and JSON |
| `analysis/e18_robust_posterior_calibration.py` | Local and split-conformal interval evaluation | E18 CSV and JSON |
| `analysis/e22_love_forward_validation.py` | Transfer-matrix and 2D SH forward checks | E22 CSV, JSON, NPZ and Figure 6 |
| `analysis/e24_wavefield_tracking_inversion.py` | End-to-end tracking of wavefield images | E24 CSV and JSON |
| `analysis/e25_qedispinv_validation.py` | Independent QEDispInv root comparison | E25 CSV, JSON and Figure 7 |
| `analysis/e26_additional_picker_baseline.py` | Connected-component segmentation baseline | E26 CSV and JSON |
| `analysis/e27_bootstrap_uncertainty.py` | Model-level bootstrap intervals | E27 JSON |
| `analysis/e28_sls_attenuation_validation.py` | Five-mechanism SLS attenuation check | E28 CSV, JSON, NPZ and Figure 8 |

## Reproducibility boundaries

- All results are synthetic or analytical; no field validation is claimed.
- Analytical mode identity is fixed by order \(n=0,1,2\), not by sorting
  velocities at each frequency.
- The SLS validation is quantified at a 280-Hz source frequency. Its material
  fit spans 80-550 Hz, but this is not a multi-source broadband inversion.
- The global tracker is a transparent baseline, not a line-by-line
  reproduction of a published segmentation or deep-learning method.
- `QEDispInv` is an external program. Its local build notes and the upstream
  license discrepancy must be reviewed before redistribution.
