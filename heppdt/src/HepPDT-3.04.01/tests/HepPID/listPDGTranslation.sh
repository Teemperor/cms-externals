#! /bin/sh
# tests/HepPID/listPDGTranslation.sh.  Generated from listPDGTranslation.sh.in by configure.

rm -f listPDGTranslation.out

./listPDGTranslation 

cmd=`diff -q -b listPDGTranslation.out ./listPDGTranslation.output`

if [ -n "$cmd" ]; then
  echo "listPDGTranslation.out and ./listPDGTranslation.output differ"
  exit 1;
fi

