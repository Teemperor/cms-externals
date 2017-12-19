#! /bin/sh
# tests/HepPID/listHerwigTranslation.sh.  Generated from listHerwigTranslation.sh.in by configure.

rm -f listHerwigTranslation.out

./listHerwigTranslation 

cmd=`diff -q -b listHerwigTranslation.out ./listHerwigTranslation.output`

if [ -n "$cmd" ]; then
  echo "listHerwigTranslation.out and ./listHerwigTranslation.output differ"
  exit 1;
fi

