.PHONY: all clean

all: main.pdf

appendix.pdf: appendix.tex
	latexmk -pdf -interaction=nonstopmode -halt-on-error appendix.tex

main.pdf: main.tex checklist.tex appendix.pdf neurips_2026.sty ../references.bib
	latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex

clean:
	latexmk -C main.tex
	latexmk -C appendix.tex
