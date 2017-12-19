#! /bin/sh
# tests/HepPDT/listEvtGenNames.sh.  Generated from listEvtGenNames.sh.in by configure.

rm -f listEvtGenNames.out

./listEvtGenNames 

if ( ! `diff -q -b listEvtGenNames.out ./listEvtGenNames.output > /dev/null` )
then
  echo "listEvtGenNames.out and ./listEvtGenNames.output differ"
  exit 1;
fi

exit 0;
