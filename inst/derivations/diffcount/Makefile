# A Makefile for generating the PDF from the Latex.
# 
# Simply type 'make' in this directory to compile the Latex into a
# single PDF.
#
# Type 'make clean' to discard the generated PDF, and all accompanying
# output.
#

# RULES
# -----
all: diffcount.pdf

# Create the main PDF file.
diffcount.pdf : diffcount.tex diffcount.bib
	pdflatex diffcount
	bibtex diffcount
	pdflatex diffcount
	pdflatex diffcount

clean:
	rm -f diffcount.aux diffcount.pdf diffcount.log
	rm -f diffcount.blg diffcount.out diffcount.bbl
	rm -f diffcount.thm

