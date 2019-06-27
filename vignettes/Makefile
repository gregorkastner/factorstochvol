all: pdf

pdf: paper.tex
	pdflatex paper
	bibtex paper
	pdflatex paper
	pdflatex paper

clean:
	rm -rfv ./*.{aux,bbl,blg,fdb_latexmk,fls,log,out,tex}

distclean: clean
	rm -rfv Figures
	rm -fv paper.Rnw
	rm -fv ref.bib
