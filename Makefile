all: inc-3.pdf

inc-3.pdf: inc-3.tex inc-3.bbl reliabilityanalysis/prismmodel.py
	pdflatex inc-3 && \
	sage inc-3.sagetex.sage && \
	pdflatex inc-3 && \
	pdflatex inc-3

reliabilityanalysis/prismmodel.py: reliabilityanalysis/prismmodel.sage
	(cd reliabilityanalysis && make prismmodel.py)

inc-3.bbl: inc-3.bib
	pdflatex inc-3 && \
	bibtex inc-3 && \
	pdflatex inc-3 && \
	pdflatex inc-3

.PHONY: clean

clean:
	rm -vf inc-3.log inc-3.pdf \
	inc-3.bbl inc-3.blg inc-3.aux \
	inc-3.sagetex.py inc-3.sagetex.sage \
	inc-3.sagetex.scmd inc-3.sagetex.sout \
	sage-plots-for-inc-3.tex/*.png \
	sage-plots-for-inc-3.tex/*.pdf
