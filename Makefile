CC=gcc
PROG=chamhash
FLAGS=-O3 -g
LIB=-lgmp -lbsd -lm

all:
	pdflatex report.tex
	bibtex report
	pdflatex report.tex
	xpdf report.pdf
claw_free_square:
	${CC} ${FLAGS} src/mod_math.c src/claw_free_square.c -o ${PROG} ${LIB}
claw_free_rsa:
	${CC} ${FLAGS} src/mod_math.c src/claw_free_rsa.c -o ${PROG} ${LIB}
clean:
	rm -rf *~ *.bbl *.aux *.blg *.dvi *.log 
