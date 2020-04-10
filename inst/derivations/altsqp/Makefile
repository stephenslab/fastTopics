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
all: altsqp.pdf

# Create the main PDF file.
altsqp.pdf : altsqp.tex altsqp.bib
	pdflatex altsqp
	bibtex altsqp
	pdflatex altsqp
	pdflatex altsqp

clean:
	rm -f altsqp.aux altsqp.pdf altsqp.log
	rm -f altsqp.blg altsqp.out altsqp.bbl
	rm -f altsqp.thm

