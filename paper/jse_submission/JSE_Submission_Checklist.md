# Journal of Seismic Exploration submission package

Generated from `manuscript_en_v0.4.md`, including the corrected phase-shift
coherence calculation and regenerated Figure 6.

## Prepared files

- `JSE_Title_Page_and_Back_Matter_v1.docx`
- `JSE_Main_Manuscript_v1.docx`
- `JSE_Cover_Letter_v1.docx`
- `JSE_Graphical_Abstract_v1.png`
- `JSE_Supplementary_Figure_S1.docx`
- `assets/figure_1.png` through `assets/figure_8.png`
- `build_jse_docx.py` (reproducible Word conversion)
- `rendered/JSE_Main_Manuscript_v1.pdf`
- `rendered/JSE_Title_Page_and_Back_Matter_v1.pdf`
- `rendered/JSE_Cover_Letter_v1.pdf`
- `rendered/JSE_Supplementary_Figure_S1.pdf`
- `JSE_Submission_Package_v1.zip`
- `JSE_Submission_Package_v1.sha256`

## Current compliance checks

- Article type: original research article.
- Title: 10 words; journal recommendation is no more than 16.
- Abstract: exactly 250 words; the author guide states no more than 250.
- Keywords: seven.
- Main sections: Introduction; Physical model and methods; Numerical
  experiments and results; Discussion; Conclusion.
- Equations: 11 numbered display equations.
- Figures: eight, inserted after their first textual citation.
- Tables: seven editable Word tables with captions above.
- References: 51, all cited and numbered in first-appearance order.
- Graphical abstract: 2656 x 1062 pixels, above the minimum size stated in the
  author guide.
- Author, affiliation, corresponding email and funding are present.
- Conflict, CRediT author contribution, data availability, ethics, consent and
  generative-AI disclosures are present.

## Journal facts requiring author confirmation

- The journal is Gold Open Access.
- The listed article processing charge is USD 1,460.
- An ORCID is suggested but was not supplied.
- The final submission portal may request suggested reviewers.
- The author should open the Word files once in Microsoft Word or a compatible
  editor to verify fonts, equation placement, image cropping and page breaks
  before upload.

## LibreOffice export check

- The three DOCX files were opened and exported by LibreOffice in headless mode.
- Main manuscript: 32 PDF pages with eight embedded images.
- Title page and back matter: two PDF pages.
- Cover letter: one PDF page.
- The generated PDFs retained all eight image objects and no conversion command
  reported an error.
- The corrected Figure 6 uses the observed complex receiver spectrum and the
  negative steering sign required by NumPy's positive-frequency FFT convention.
- The potentially misleading raw coherence heatmap is moved to Supplementary
  Figure S1; main Figure 6 now shows only theory-consistent detected peaks.

## Scientific checks already completed

- All analytical, tracking, inversion, calibration and bootstrap experiments
  were rerun from the public repository.
- The eight source SVGs passed XML validation and were rendered to PNG at a
  resolution sufficient for the intended Word layout.
- No field-data accuracy is claimed.
