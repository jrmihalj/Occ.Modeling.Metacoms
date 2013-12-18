# Makefile for our metacommunity occupancy model paper

all: output/article_text.docx output/article_text.pdf 

output/article_text.docx: manuscript/references.bib manuscript/manuscript.md
	pandoc -H format.sty -V fontsize=12pt --bibliography manuscript/references.bib --csl=ecology.csl manuscript/manuscript.md -o output/article_text.docx

output/article_text.pdf: manuscript/references.bib manuscript/manuscript.md
	pandoc -H format.sty -V fontsize=12pt --bibliography manuscript/references.bib --csl=ecology.csl manuscript/manuscript.md -o output/article_text.pdf
