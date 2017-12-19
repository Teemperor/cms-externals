#! /bin/sh
# tests/HepPID/listEvtGenTranslation.sh.  Generated from listEvtGenTranslation.sh.in by configure.

rm -f listEvtGenTranslation.out

./listEvtGenTranslation 

cmd=`diff -q -b listEvtGenTranslation.out ./listEvtGenTranslation.output`

if [ -n "$cmd" ]; then
  echo "listEvtGenTranslation.out and ./listEvtGenTranslation.output differ"
  exit 1;
fi

