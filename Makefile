all:
	pdflatex report.tex
	bibtex report
	pdflatex report.tex
	xpdf report.pdf
clean:
	rm -rf *~ *.bbl *.aux *.blg *.dvi *.log 
