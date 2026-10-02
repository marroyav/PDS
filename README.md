# Photon-detector background model

Analytical background estimates for ProtoDUNE-VD, DUNE FD-VD, and the DUNE
Near Detector 2x2 demonstrator.

## Build

Install a TeX distribution such as TeX Live or MacTeX with `latexmk`,
`pdflatex`, and the packages `geometry`, `amsmath`, `amssymb`, `bm`,
`booktabs`, `siunitx`, and `hyperref`. Then run:

```sh
make
```

The PDF is written to `output/pdf/background_model.pdf`. The main document
includes `nd_2x2_background.tex`; that file is a section, not a standalone
LaTeX document. `latexmk` runs the passes needed to resolve cross-references.
The build ignores user-level latexmk configuration for reproducibility.

To remove auxiliary build files while keeping the PDF:

```sh
make clean
```

Generated output is excluded from Git.

## Reproduce the numerical reasoning

With Python 3, NumPy, and SciPy available:

```sh
python3 tools/check_background_model.py
```

The script reproduces the centered-XA direct-light integrals, finite-aperture
corrections, atmospheric Ar-39 normalization, and the ND numerical illustration.
It also checks the solid angle against independent surface integration,
limiting cases, and quadrature convergence. The review was run with NumPy
2.1.3 and SciPy 1.15.2; no experimental datasets are used.

The geometry-only FD/ProtoDUNE ratio is 2.912 for the stated reference boxes
and centered small apertures (2.923 with finite 0.6 m square apertures).
This is not a measured total-background ratio. Absolute PE rates require
scintillation yield, actual optical geometry, transport, and calibrated PDE.
The note distinguishes primary PE, charge including correlated SiPM noise,
tagged-event light, and trigger rates, and states the assumptions behind
the ND coverage illustration.
