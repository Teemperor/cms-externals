#! /bin/sh
# tests/HepPDT/listPythiaNames.sh.  Generated from listPythiaNames.sh.in by configure.

rm -f listPythiaNames.out

./listPythiaNames 
diff -q -b listPythiaNames.out ./listPythiaNames.output > /dev/null

