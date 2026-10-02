.PHONY: all pdf clean

all: pdf

pdf:
	latexmk -norc -pdf -interaction=nonstopmode -halt-on-error -outdir=output/pdf background_model.tex

clean:
	latexmk -norc -c -outdir=output/pdf background_model.tex
