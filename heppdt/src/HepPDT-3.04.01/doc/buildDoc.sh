#! /bin/sh

# this script builds the documents
# the bootstrap process executes this script

set -x
# run doxygen if present
if [ `whereis doxygen | wc -w` -gt 1 ]; then
   doxygen doxygen.conf
   if [ `whereis latex | wc -w` -gt 1 ]; then
      cd latex
      pdflatex refman.tex
      pdflatex refman.tex
      pdflatex refman.tex
      pdflatex refman.tex
      makeindex refman.idx
      pdflatex refman.tex
      mv refman.pdf ../HepPDT/HepPDT_reference_manual.pdf
      cd ..
      rm -rf latex
   fi
fi
# now run latex
if [ `whereis latex | wc -w` -gt 1 ]; then
   cd HepPDT
   latex HepPDT.tex 
   latex HepPDT.tex
   dvips HepPDT
   dvipdf HepPDT
   cd ../HepPID
   latex HepPID.tex 
   latex HepPID.tex
   dvips HepPID
   dvipdf HepPID
   cd ..
fi
# cleanup
rm -f H*/*.aux H*/*.dvi H*/*.log H*/*.toc
