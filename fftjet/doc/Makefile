all: fftjet_manual.pdf

fftjet_manual.pdf: fftjet_manual.ps
	ps2pdf $<

%.ps : %.dvi
	dvips -t letter -o $@ $<

%.dvi : %.tex
	latex $<
	latex $<

fftjet_manual.dvi: fftjet_manual.tex
	latex fftjet_manual.tex
	makeindex fftjet_manual.idx
	latex fftjet_manual.tex
	latex fftjet_manual.tex

clean:
	rm -f fftjet_manual.ps
	rm -f *~ *.log *.aux *.dvi *.bak *.ind *.idx *.toc *.ilg
