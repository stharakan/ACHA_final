FILE=presentation

all: $(FILE).pdf

.PHONY: clean

clean:
		\rm *.aux *.nav *.out *.snm *.toc *.log

$(FILE).pdf: $(FILE).tex
		pdflatex $(FILE).tex
