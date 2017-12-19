#! /bin/sh
# tests/HepPDT/listPDGNames.sh.  Generated from listPDGNames.sh.in by configure.

rm -f listPDGNames.out

./listPDGNames 
diff -q -b listPDGNames.out ./listPDGNames.output > /dev/null

