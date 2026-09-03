.PHONY: all clean

FIGURE_ASSETS := $(wildcard assets/*.png)

all: solvai_main.pdf

appendix.pdf: appendix.tex $(FIGURE_ASSETS)
	latexmk -pdf -interaction=nonstopmode -halt-on-error appendix.tex

solvai_main.pdf: main.tex checklist.tex appendix.pdf neurips_2026.sty ../references.bib $(FIGURE_ASSETS)
	latexmk -pdf -jobname=solvai_main -interaction=nonstopmode -halt-on-error main.tex

clean:
	latexmk -C -jobname=solvai_main main.tex
	latexmk -C appendix.tex
