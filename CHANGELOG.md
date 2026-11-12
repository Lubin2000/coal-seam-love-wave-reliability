# Changelog

## 0.4.0 - 2026-11-10

- Corrected the independent phase-shift coherence calculation to retain the
  observed complex receiver spectrum.
- Corrected the steering-sign convention for NumPy's positive-frequency FFT.
- Corrected Figure 6 velocity-cell coordinate mapping.
- Regenerated the independent-wavefield and end-to-end tracking results.
- Added manuscript v0.4 and the JSE submission package.

## 0.3.0 - 2026-11-06

- Adopted MIT for original software and CC BY 4.0 for original data, figures,
  and documentation.
- Public release archived by Zenodo: version DOI
  `10.5281/zenodo.24054134`; concept DOI `10.5281/zenodo.24054133`.
- Fixed analytical mode ordering so that modes 0, 1 and 2 retain their physical
  branch identity across cutoff frequencies.
- Updated the tracker, information analysis, robust inversion and calibration
  experiments after the branch-order correction.
- Replaced first-order stress-memory integration with an auxiliary-memory,
  Crank-Nicolson SLS formulation.
- Added the 80-550 Hz SLS material fit and 280-Hz paired-wavefield attenuation
  validation.
- Added the current English manuscript, eight publication figures, evidence
  reports and a one-command experiment runner.

## 0.2.0

- Added connected-component segmentation, bootstrap uncertainty,
  cross-domain conformal calibration, QEDispInv and independent SH-wavefield
  checks.

## 0.1.0

- Initial controlled synthetic benchmark for multimodal Love-wave picking and
  robust inversion.
