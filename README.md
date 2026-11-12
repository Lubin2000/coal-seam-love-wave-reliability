# Coal-seam Love-wave mode reliability analysis

[![DOI](https://zenodo.org/badge/1449400437.svg)](https://doi.org/10.5281/zenodo.24054133)

This repository contains the code, synthetic results, figures and manuscript
materials for a controlled study of higher-mode Love-wave information,
automatic dispersion tracking, robust inversion and interval reliability in a
three-layer coal-seam model.

**Author:** Lu Bin  
**Affiliation:** Wenzhou University of Technology, Wenzhou, Zhejiang, China  
**Corresponding email:** lubin2000@163.com  
**Funding:** National Natural Science Foundation of China, Grant No. 41974209

**Public repository:** <https://github.com/Lubin2000/coal-seam-love-wave-reliability>

**Previous version DOI (v0.3.0):** <https://doi.org/10.5281/zenodo.24054134>  
**Concept DOI (all versions):** <https://doi.org/10.5281/zenodo.24054133>  
**Current version:** v0.4.0, which corrects the independent phase-shift image
calculation. Its version DOI will be added after archival.

## Scientific scope

- The study is entirely analytical or synthetic.
- No field-data accuracy or mine-scale validation is claimed.
- Analytical modal identity is fixed by physical order \(n=0,1,2\). Roots are
  not reordered by phase velocity at each frequency.
- The main inversion uses a first-order energy-weighted modal-Q model.
- A five-mechanism SLS forward check is quantified at a 280-Hz source
  frequency over material fits spanning 80-550 Hz.

## Repository layout

```text
analysis/       Core Python scripts and their numerical outputs
paper/          Current manuscript, figures, references and evidence reports
paper/jse_submission/  Current Journal of Seismic Exploration upload package
docs/           Licensing, third-party and release documentation
run_core_experiments.ps1
requirements.txt
CITATION.cff
```

## Environment

Tested on Windows with Python 3.13 and NumPy 2.x.

```powershell
py -m pip install -r requirements.txt
```

The complete finite-difference validation takes several minutes. To redraw
figures from existing results:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_core_experiments.ps1 -FiguresOnly
```

To rerun all analytical and inversion experiments while skipping the expensive
wavefields:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_core_experiments.ps1 -SkipWavefields
```

To rerun everything:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_core_experiments.ps1
```

## Current manuscript

`paper/manuscript_en_v0.4.md`

The current submission package is in `paper/jse_submission/`.

## License

- Original software code is licensed under the MIT License. See
  `LICENSE-CODE`.
- Original data, figures and documentation are licensed under CC BY 4.0. See
  `LICENSE-DATA`.

Third-party material is excluded. See `LICENSE.md` and `docs/THIRD_PARTY.md`.

## Publication status

The repository is public at
<https://github.com/Lubin2000/coal-seam-love-wave-reliability>. Version v0.3.0
remains archived at <https://doi.org/10.5281/zenodo.24054134>. Version v0.4.0
corrects the phase-shift coherence implementation and is the release candidate
for the submitted manuscript. The concept DOI
<https://doi.org/10.5281/zenodo.24054133> identifies the version family.

## Third-party software

QEDispInv is not redistributed in this repository. Its input configuration,
reference output and build notes are included for audit. See
`docs/THIRD_PARTY.md`.
