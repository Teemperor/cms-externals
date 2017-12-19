all:
	latex tauspinner-application-studies.tex
	bibtex tauspinner-application-studies
	latex tauspinner-application-studies
	latex tauspinner-application-studies
	dvipdf tauspinner-application-studies.dvi
	@echo "##################################################################"
	@echo "WARNING: You have processed your analysis into Apenndix A.        "
	@echo "         In this way you have modified the content of your copy   "
	@echo "         of paper IFJPAN-IV-2013-19. Your results collected in    "
	@echo "         Appendix A may differ from the one of paper if your      "
	@echo "         sample was of distinct physics nature.                   "
	@echo "         main text of  tauspinner-application-studies.pdf will    "
	@echo "         include some useful explanations also in this case.      "
	@echo "##################################################################"

user-analysis:
	latex user-analysis.tex
	latex user-analysis.tex
	dvipdf user-analysis.dvi
	@echo "##################################################################"
	@echo "WARNING: You have processed your analysis into user-analysis.pdf  "
	@echo "         It is to systematize your results. Alternatively         "
	@echo "         tauspinner-application-studies.pdf  can be constructed   "
	@echo "         with 'make all'. Then your results will be placed in     "
	@echo "         appendix A of paper IFJPAN-IV-2013-19                    "
	@echo "##################################################################"

clean:
	rm -f *.toc *.aux *.dvi *.pdf *.out *.log *.blg *.bbl
