# Revision note: fixed modal identity and regenerated results

## Reason for revision

The original higher-mode tracker represented a state as the list of currently active velocities. When the fundamental or first-higher branch disappeared, the remaining velocities shifted to lower list positions and could be assigned the wrong mode number.

This was detected during the independent wavefield end-to-end test, where a single visible high-velocity branch was repeatedly labelled as mode 1 because a low-velocity body-wave ridge occupied the second slot.

## Code change

`e16_robust_mode_tracking.py` now uses three fixed mode slots:

- mode 0: fundamental;
- mode 1: first higher;
- mode 2: second higher.

Each slot can independently contain a velocity or be absent. Continuity is evaluated only when the same mode is present at adjacent frequencies. An appearance penalty applies when a mode starts at a frequency during high-to-low tracking.

Candidate peaks are separated before dynamic programming to prevent several minor maxima on the same broad ridge from being labelled as different modes.

## Regenerated outputs

- `e16_robust_mode_tracking.py`
- `e16_robust_tracking_runs.csv`
- `e16_robust_tracking_summary.json`
- `e17_robust_likelihood.py`
- `e17_robust_likelihood_runs.csv`
- `e17_robust_likelihood_summary.json`
- `e18_robust_posterior_calibration.py`
- `e18_robust_posterior_calibration_runs.csv`
- `e18_robust_posterior_calibration_summary.json`
- `paper/figures/fig3_mode_picking_en.svg`
- `paper/figures/fig4_robust_likelihood_en.svg`
- `paper/figures/fig5_calibration_en.svg`

## Updated principal results

- Global tracking fundamental detection: 93.8%.
- Global tracking first-higher detection: 99.6%.
- Global tracking second-higher detection: 69.9%.
- Global tracking false-pick rate: 35.4%.
- Student-t thickness RMSE with 30% outliers and phase only: 0.102 m.
- Student-t thickness RMSE with 30% outliers and attenuation: 0.101 m.
- Raw local 90% coverage across 140 held-out synthetic cases: 37.9% thickness, 33.6% velocity, 85.0% Q.
- Global split-conformal coverage: 83.6% thickness, 80.0% velocity, 92.9% Q.

## Interpretation

The corrected mode-slot formulation improves higher-mode retention but does not solve direct transfer to finite-difference images. Broad ridges, non-modal energy and label changes across gaps remain substantial. The paper therefore keeps the synthetic-image results as a controlled benchmark and reports the end-to-end wavefield degradation separately.

## v0.2 additions

- Added a connected-component image-segmentation baseline on the same 50 images.
- Added 2,000 model-level bootstrap confidence intervals for picking, inversion and calibration metrics.
- Added cross-domain conformal calibration tests across thin/thick seams, low/high velocities and low/high outlier regimes.
- Added an independent QEDispInv root comparison at 34 coal-mode points.
- Added an independent 2D SH finite-difference forward and end-to-end tracking check.
- Updated Figure 3 to include the segmentation baseline.
- Added Figure 5 (calibration), Figure 6 (independent wavefield) and Figure 7 (QEDispInv).

The v0.2 English draft is retained as a historical audit version.

## v0.3 correction: physical mode order and viscoelastic validation

### Additional branch-order defect

The fixed-slot tracker was not sufficient by itself. The analytical root
function returned all available roots sorted by decreasing phase velocity.
Before a higher branch appeared, the fundamental branch occupied root index 0.
After a higher branch appeared near its cutoff, the descending order could
place that branch at index 0 and shift the fundamental branch to index 1.

This produced an artificial mode relabeling in the forward responses, group
velocity derivatives, attenuation calculations and all downstream scripts.
It was exposed when the SLS attenuation test isolated the high-velocity
near-cutoff branch and overestimated attenuation by about 61%.

The analytical and independent transfer-matrix root lists are now sorted by
physical mode order. Mode 0 is the \(n=0\) fundamental branch, mode 1 is the
\(n=1\) first higher branch and mode 2 is the \(n=2\) second higher branch.
Velocity-based reordering is no longer used.

### SLS implementation correction

The first stress-memory update used a first-order Euler approximation. A
harmonic constitutive test showed that it reproduced Q=24.8 when Q=40 was
requested. The solver now uses auxiliary SLS memory variables and a
Crank-Nicolson update. The same test reproduces Q=38.50 at 280 Hz for the
Q=40 fit and Q=288.32 for the Q=300 fit.

Layer moduli are normalized at 280 Hz so that changing Q does not intentionally
change the reference phase velocity. The five-mechanism fit has approximately
4% mean relative Q error over 80-550 Hz.

### Regenerated headline results

- First-higher theoretical thickness benefit: 33.9%, not 44%.
- Joint higher-mode/group/attenuation benefit: 62.2%, not 72%.
- Global fixed-slot detection: 97.2% / 92.3% / 32.2% for modes 0/1/2.
- Global false-pick rate: 40.8%.
- Segmentation detection: 91.5% / 88.2% / 78.5% for modes 0/1/2.
- Student-t phase-only thickness RMSE at 30% outliers: 0.548 m.
- Student-t thickness RMSE with attenuation at 30% outliers: 0.378 m.
- Relative Q error with attenuation at 30% outliers: 4.9%.
- Raw local coverage: 17.9% thickness, 34.3% velocity and 56.4% Q.
- Global conformal coverage: 79.3% thickness, 70.0% velocity and 88.6% Q.
- SLS effective-Q check at 280 Hz: 40.92 predicted and 37.35 observed, an 8.7% difference.

### New artifacts

- `e28_sls_attenuation_validation.py`
- `e28_sls_attenuation_validation.csv`
- `e28_sls_attenuation_validation_summary.json`
- `e28_sls_wavefields.npz`
- `paper/figures/fig8_sls_attenuation_en.svg`
- `paper/manuscript_en_v0.3.md`

The previous v0.2 draft is retained for audit, but its numerical results are
superseded where they depend on analytical mode indices or SLS attenuation.

## Author and funding metadata

- Author and corresponding author: Lu Bin.
- Affiliation: Wenzhou University of Technology, Wenzhou, Zhejiang, China.
- Corresponding email: lubin2000@163.com.
- Funding: National Natural Science Foundation of China, Grant No. 41974209.
- Author contribution: sole-author CRediT statement added.
