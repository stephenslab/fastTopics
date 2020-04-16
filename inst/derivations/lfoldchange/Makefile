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
all: lfoldchange.pdf

# Create the main PDF file.
lfoldchange.pdf : lfoldchange.tex lfoldchange.bib
	pdflatex lfoldchange
	bibtex lfoldchange
	pdflatex lfoldchange
	pdflatex lfoldchange

clean:
	rm -f lfoldchange.aux lfoldchange.pdf lfoldchange.log
	rm -f lfoldchange.blg lfoldchange.out lfoldchange.bbl
	rm -f lfoldchange.thm

