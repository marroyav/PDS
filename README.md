# Photon-detector background model

Analytical background estimates for ProtoDUNE-VD, one DUNE FD-VD module,
full ND-LAr, and the ND 2x2 demonstrator.

## Per-channel comparison

The PDF opens with numerical rates for every detector and readout technology.
The reference **Ar-39 primary scintillation** estimates are:

| Detector | PDS electronic channels | Primary PE/s/channel |
| --- | ---: | ---: |
| ProtoDUNE-VD | 32 | 39,000 wall; 77,000 cathode |
| FD-VD | 1,344 | 110,000 wall; 225,000 cathode |
| Full ND-LAr | 8,400 | 44 ArCLight; 133 LCM; **89 mean** |
| ND 2x2 | 384 | 15 ArCLight; 48 LCM; **31 mean** |

Full ND-LAr has 35 modules with 240 light channels each; 2x2 has four
smaller modules with 96 each. Counts follow the cited designs, not connected
channels in a particular run. The ND source masses are 146.5 t and 2.4 t.

These are conditional reference estimates, not measured trigger rates.
VD values use representative positions, not installed-array averages.
ND values use uniform first-hit wall flux. The common reference yield is
4,400 photons/Ar-39 decay; whole-device PDEs are 4.5% XA, 0.2% ArCLight,
and 0.6% LCM. The PDF cites the inputs, explains channel sharing, and gives
the rescaling factors. It also separates total Bq divided by channels from
detected light and gives site-dependent, track-only cosmic estimates.

## Build

Install a TeX distribution such as TeX Live or MacTeX with `latexmk`,
`pdflatex`, and the packages `geometry`, `amsmath`, `amssymb`, `bm`,
`booktabs`, `siunitx`, and `hyperref`. Then run:

```sh
make
```

The PDF is written to `output/pdf/background_model.pdf`. The main document
includes `per_channel_estimates.tex`, `cosmic_estimates.tex`, and
`nd_2x2_background.tex` (the latter covers both ND sizes). These are sections,
not standalone LaTeX documents. `latexmk` resolves the cross-references.
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
corrections, atmospheric Ar-39 normalization, and all per-channel tables.
It also checks the solid angle against independent surface integration,
limiting cases, and quadrature convergence. The review was run with NumPy
2.1.3 and SciPy 1.15.2; no experimental datasets are used.

The geometry-only FD/ProtoDUNE ratio is 2.912 for the stated reference boxes
and centered small apertures (2.923 with finite 0.6 m square apertures).
This is not a measured total-background ratio. The absolute estimates use
explicit reference yield, geometry, transport, and PDE assumptions.
The note distinguishes primary PE, charge including correlated SiPM noise,
tagged-event light, and trigger rates, and states the assumptions behind
the ND coverage illustration.
