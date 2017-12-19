#! /bin/sh
# tests/HepPID/listPythiaTranslation.sh.  Generated from listPythiaTranslation.sh.in by configure.

rm -f listPythiaTranslation.out

./listPythiaTranslation 

cmd=`diff -q -b listPythiaTranslation.out ./listPythiaTranslation.output`

if [ -n "$cmd" ]; then
  echo "listPythiaTranslation.out and ./listPythiaTranslation.output differ"
  exit 1;
fi
