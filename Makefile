# run `make vignettes && make clean` before `R CMD build ...`

vignettes: vignettes/paper.Rtex

vignettes/paper-original.Rnw: paper.Rnw
	cp paper.Rnw vignettes/paper-original.Rnw

vignettes/paper-Sweave.Rtex: vignettes/paper-original.Rnw
	R -e 'knitr::knit("vignettes/paper-original.Rnw")'
	mv paper-original.tex vignettes/paper-original.tex
	mv vignettes/paper-original.tex vignettes/paper-Sweave.Rtex
	mv Figures vignettes/Figures

vignettes/paper.Rtex: vignettes/paper-knitr.Rtex
	cp vignettes/paper-knitr.Rtex vignettes/paper.Rtex

vignettes/paper-knitr.Rtex: vignettes/paper-Sweave.Rtex
	R -e 'knitr::Sweave2knitr("vignettes/paper-Sweave.Rtex")'
	mv vignettes/paper-Sweave-knitr.Rtex vignettes/paper-knitr.Rtex

clean:
	rm -rf \
		vignettes/paper-original.Rnw \
		vignettes/paper-knitr.Rtex \
		vignettes/paper-Sweave.Rtex \
		vignettes/Sweave.sty

distclean:
	rm -rf \
		vignettes/paper-original.Rnw \
		vignettes/paper-knitr.Rtex \
		vignettes/paper-Sweave.Rtex \
		vignettes/Sweave.sty \
		vignettes/paper.Rtex \
		vignettes/Figures

