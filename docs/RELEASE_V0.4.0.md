# GitHub release v0.4.0

## Suggested tag

`v0.4.0`

## Suggested title

`v0.4.0 - Corrected Phase-Shift Validation`

## Suggested release notes

```markdown
Correction release accompanying manuscript v0.4.

## Corrections

- Retained the observed complex receiver spectrum in the phase-shift image.
- Applied the steering sign consistent with NumPy's positive-frequency FFT.
- Replaced 2-m receiver spacing with 1-m spacing to avoid high-frequency
  spatial aliasing while retaining the 200-m aperture.
- Corrected the Figure 6 velocity-coordinate mapping.
- Regenerated the independent-wavefield and end-to-end inversion results.

## Presentation

- Main Figure 6 now shows only analytical branches and detected peaks.
- The raw coherence image is retained as Supplementary Figure S1 because it
  contains strong nonmodal energy.

## Submission package

The Journal of Seismic Exploration Word, PDF and graphical materials are
included under `paper/jse_submission/`.

Version DOI: to be assigned by Zenodo for v0.4.0.
Concept DOI: https://doi.org/10.5281/zenodo.24054133
```

After publishing:

1. Wait for Zenodo to archive the release and issue the version DOI.
2. Add the DOI to `README.md`, `CITATION.cff`, and manuscript v0.4.
3. Update the submission PDFs.
4. Commit and push the DOI metadata update.
