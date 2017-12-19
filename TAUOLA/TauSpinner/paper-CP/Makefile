PAPER=CP-spiner

all:
	latex $(PAPER).tex
	bibtex $(PAPER)
	bibtex $(PAPER)
	latex $(PAPER).tex
	latex $(PAPER).tex
	dvipdf $(PAPER).dvi

clean:
	rm -f $(PAPER).dvi $(PAPER).pdf $(PAPER).log $(PAPER).aux $(PAPER).bbl $(PAPER).blg
