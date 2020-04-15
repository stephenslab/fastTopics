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
all: tclr.pdf

# Create the main PDF file.
tclr.pdf : tclr.tex tclr.bib
	pdflatex tclr
	bibtex tclr
	pdflatex tclr
	pdflatex tclr

clean:
	rm -f tclr.aux tclr.pdf tclr.log
	rm -f tclr.blg tclr.out tclr.bbl
	rm -f tclr.thm

