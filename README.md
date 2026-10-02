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

Generated output is excluded from Git. Compilation verifies document syntax
and references; it does not validate the physical assumptions or estimates.
