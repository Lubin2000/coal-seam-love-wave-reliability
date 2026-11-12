# Reviewer-Facing Evidence Addendum

## 1. Additional image-segmentation baseline

A connected-component segmentation picker was implemented on the same 50
synthetic dispersion images used for independent-peak, greedy and global
tracking.

| Method | Mode 0 | Mode 1 | Mode 2 | False picks | Median thickness interval |
|---|---:|---:|---:|---:|---:|
| Independent peaks | 60.6% | 79.6% | 93.8% | 47.6% | 0.0343 m |
| Greedy | 83.4% | 71.8% | 42.1% | 29.5% | 0.0384 m |
| Global fixed-slot | 97.2% | 92.3% | 32.2% | 40.8% | 0.0336 m |
| Connected-component segmentation | 91.5% | 88.2% | 78.5% | 40.0% | 0.0339 m |

The segmentation baseline recovers more second-higher picks than the global
tracker, but its false-pick rate and median parameter interval are similar.
Raw mode count is therefore not a sufficient quality criterion.

## 2. Bootstrap uncertainty

Two thousand model-level bootstrap resamples were used, with the model as the
resampling unit.

Selected 95% intervals:

- Global second-higher detection: 25.2%-39.5%.
- Global false-pick rate: 38.0%-44.3%.
- Global median thickness interval: 0.0296-0.0369 m.
- Segmentation second-higher detection: 73.8%-83.7%.
- Segmentation false-pick rate: 37.6%-43.0%.
- Student-t thickness RMSE at 30% outliers, phase only: 0.442-0.647 m.
- Student-t velocity RMSE at 30% outliers, phase only: 37.2-61.4 m/s.
- Student-t thickness RMSE at 30% outliers with attenuation: 0.252-0.490 m.
- Student-t relative Q error at 30% outliers with attenuation: 3.5%-6.5%.
- Conformal thickness coverage: 66.4%-90.7%.
- Conformal velocity coverage: 55.7%-82.9%.
- Conformal Q coverage: 80.0%-95.0%.

## 3. Cross-domain conformal calibration

Random same-family calibration reaches approximately 79% thickness, 70%
velocity and 89% Q coverage, with much wider intervals than the raw local
posterior. This does not transfer reliably when the
calibration and test regimes differ.

| Transfer scenario | Raw h / Vs coverage | Conformal h / Vs coverage | Conformal mean width h / Vs |
|---|---|---|---|
| Thin to thick | 20.7% / 42.4% | 76.1% / 77.2% | 0.859 m / 53.2 m/s |
| Thick to thin | 12.5% / 18.8% | 100.0% / 56.3% | 1.549 m / 102.1 m/s |
| Low to high velocity | 14.3% / 27.4% | 41.7% / 44.0% | 0.560 m / 63.2 m/s |
| High to low velocity | 23.2% / 44.6% | 100.0% / 100.0% | 1.338 m / 51.8 m/s |
| Low to high outlier score | 4.7% / 15.6% | 35.9% / 35.9% | 0.373 m / 56.3 m/s |

Conformal calibration is therefore useful for exchangeable data but cannot be
presented as a universal field guarantee. Some transfers under-cover; others
achieve coverage only through uninformatively wide intervals.

## 4. Broadband viscoelastic attenuation check

The five-mechanism SLS fit has a mean relative Q error of 4.0% over 80-550 Hz
for both coal Q=40 and Q=300. At a 280-Hz source frequency, the modal-energy
model predicts an effective Q of 40.92, while the paired SLS wavefield gives
37.35, an 8.7% difference. This is a forward consistency check, not a
broadband Q inversion or a field calibration.

## 5. Independent phase-shift image correction

A manual audit of Figure 6 exposed two errors in the earlier phase-shift
coherence implementation:

1. The image multiplied receiver amplitudes by the trial steering phasor but
   omitted the observed complex trace spectrum. The image therefore primarily
   represented an amplitude taper rather than dispersion coherence.
2. After restoring the complex spectrum, the steering sign was inconsistent
   with NumPy's positive-frequency FFT convention.

The corrected transform keeps the complex spectrum and applies
`exp(-i 2 pi f x / c)`. The independent-wavefield results were regenerated:

- Fundamental: 7/70 strict detections (10.0%).
- First higher: 3/34 strict detections (8.8%).
- Second higher: 0/11 strict detections (0.0%).
- Median matched velocity error: 17.9 and 23.7 m/s for the first two modes.

The corrected finite-difference images were also passed through the
end-to-end tracker. Global tracking detected 18.6%, 8.8% and 0.0% of modes
0-2 at an 87.3% mean false-pick rate, with thickness and velocity MAEs of
2.22 m and 57.5 m/s.

These results supersede the Figure 6 and end-to-end values in archived
software version v0.3.0.

## Artifacts

- `e26_additional_picker_baseline.py`
- `e26_additional_picker_runs.csv`
- `e26_additional_picker_summary.json`
- `e27_bootstrap_uncertainty.py`
- `e27_bootstrap_uncertainty_summary.json`
- `e18_robust_posterior_calibration_summary.json`
- `e28_sls_attenuation_validation.py`
- `e28_sls_attenuation_validation_summary.json`
- `e28_sls_wavefields.npz`
