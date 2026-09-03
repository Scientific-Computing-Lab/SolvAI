.PHONY: all clean

all: solvai_main.pdf

appendix.pdf: appendix.tex
	latexmk -pdf -interaction=nonstopmode -halt-on-error appendix.tex

solvai_main.pdf: main.tex checklist.tex appendix.pdf neurips_2026.sty ../references.bib
	latexmk -pdf -jobname=solvai_main -interaction=nonstopmode -halt-on-error main.tex

clean:
	latexmk -C -jobname=solvai_main main.tex
	latexmk -C appendix.tex
