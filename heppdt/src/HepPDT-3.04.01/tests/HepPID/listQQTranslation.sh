#! /bin/sh
# tests/HepPID/listQQTranslation.sh.  Generated from listQQTranslation.sh.in by configure.

rm -f listQQTranslation.out

./listQQTranslation 
cmd=`diff -q -b listQQTranslation.out ./listQQTranslation.output`

if [ -n "$cmd" ]; then
  echo "listQQTranslation.out and ./listQQTranslation.output differ"
  exit 1;
fi
