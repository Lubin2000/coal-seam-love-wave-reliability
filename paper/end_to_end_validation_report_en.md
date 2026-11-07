# End-to-End Model Tracking on Independent Wavefield Images

## Purpose

E22 demonstrated that analytical Love branches produce coherent peaks in a finite-difference wavefield. This end-to-end test asks a stricter question:

> What happens when the image produced by the independent wavefield is passed directly to the same automatic tracker used by the paper?

## Data

- Three seams: 1.8, 2.8 and 3.8 m.
- Frequencies: 100-560 Hz at 20-Hz spacing.
- Images: E22 normalized phase-shift coherence.
- Comparison: analytical images generated on the same frequency-velocity grid with Gaussian modal ridges, a false ridge and random noise.
- Trackers: independent peak selection, greedy tracking and global order-constrained tracking.
- Inversion: robust Student-t fit for thickness and coal shear velocity using global-tracked picks.

The inversion used a broad global search followed by local refinement. It did not use the true thickness to initialize the search.

## Tracking results

| Image | Tracker | Mode 0 | Mode 1 | Mode 2 | False picks |
|---|---|---:|---:|---:|---:|
| Analytical | Independent peaks | 98.5% | 100.0% | 59.3% | 46.5% |
| Analytical | Greedy | 90.2% | 87.8% | 51.9% | 12.6% |
| Analytical | Global | 81.9% | 100.0% | 0.0% | 46.4% |
| Finite difference | Independent peaks | 47.5% | 9.0% | 3.7% | 79.7% |
| Finite difference | Greedy | 8.7% | 0.0% | 0.0% | 90.7% |
| Finite difference | Global | 21.1% | 9.0% | 3.7% | 84.3% |

## Global-track inversion

| Image | Thickness MAE | Shear-velocity MAE |
|---|---:|---:|
| Analytical images | 2.53 m | 107.3 m/s |
| Finite-difference images | 2.25 m | 76.7 m/s |

The finite-difference velocity MAE is not evidence of better performance: too few modes were retained for a stable three-parameter interpretation, and mode-label errors can compensate for one another in a fit.

## Interpretation

The mode-availability and information-budget conclusions remain valid because they are calculated from known modal branches and controlled image perturbations.

The automatic-tracking and inversion performance observed in synthetic images is optimistic. Direct transfer to the finite-difference image causes:

- broad and weak low-frequency ridges;
- several local peaks on a single ridge;
- body-wave and leaky-energy ridges competing with true modes;
- mode labels changing across branch gaps;
- a large false-pick fraction.

A physical order constraint helps relative to greedy tracking, but it does not eliminate non-modal coherent energy. The next tracking version must include amplitude-depth polarization, local wavenumber coherence, mode-confidence probabilities and a rejection class for non-modal arrivals.

## Artifacts

- `e24_wavefield_tracking_inversion.py`
- `e24_wavefield_tracking_inversion.csv`
- `e24_wavefield_tracking_inversion_summary.json`
