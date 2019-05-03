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
all: fast-topics.pdf

# Create the main PDF file.
fast-topics.pdf : fast-topics.tex fast-topics.bib
	pdflatex fast-topics
	bibtex fast-topics
	pdflatex fast-topics
	pdflatex fast-topics

clean:
	rm -f fast-topics.aux fast-topics.pdf fast-topics.log
	rm -f fast-topics.blg fast-topics.out fast-topics.bbl
	rm -f fast-topics.thm

