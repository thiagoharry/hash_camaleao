CC=gcc
PROG=chamhash
FLAGS=-O2 -g -Wall
LIB=-lgmp -lbsd -lm

help:
	@echo "make doc: Build documentation"
	@echo "make claw_free_square: Build c. hash based in permutation of quadratic residues"
	@echo "make claw_free_rsa: Build chameleon hash based in RSA permutation"
	@echo "make krawczyk_log: Build Krawxzyk's chameleon hash based in discrete log"
	@echo "make nyberg_rueppel: Build  chameleon hash based in Nyberg-Rueppel"
	@echo "make fiat_shamir: Build chameleon hash based in Fiat-Shamir sigma protocol"
doc:
	pdflatex report.tex
	bibtex report
	pdflatex report.tex
	xpdf report.pdf
claw_free_square:
	${CC} ${FLAGS} src/mod_math.c src/claw_free_square.c -o ${PROG} ${LIB}
claw_free_rsa:
	${CC} ${FLAGS} src/mod_math.c src/claw_free_rsa.c -o ${PROG} ${LIB}
krawczyk_log:
	${CC} ${FLAGS} src/mod_math.c src/krawczyk_log.c -o ${PROG} ${LIB}
nyberg_rueppel:
	${CC} ${FLAGS} src/mod_math.c src/nyberg_rueppel.c -o ${PROG} ${LIB}
fiat_shamir:
	${CC} ${FLAGS} src/mod_math.c src/fiat_shamir.c -o ${PROG} ${LIB}
clean:
	rm -rf *~ *.bbl *.aux *.blg *.dvi *.log 
