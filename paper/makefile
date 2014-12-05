FILE=final_paper

all: $(FILE).pdf

.PHONY: clean

clean:
		\rm *.aux *.log *.bbl *.blg *.pdf

$(FILE).pdf: $(FILE).tex *.tex
		pdflatex $(FILE)
		bibtex $(FILE)
		pdflatex $(FILE)
		pdflatex $(FILE)
