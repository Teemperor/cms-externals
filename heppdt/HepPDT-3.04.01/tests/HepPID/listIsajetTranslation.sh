#! /bin/sh
# tests/HepPID/listIsajetTranslation.sh.  Generated from listIsajetTranslation.sh.in by configure.

rm -f listIsajetTranslation.out

./listIsajetTranslation 

cmd=`diff -q -b listIsajetTranslation.out ./listIsajetTranslation.output`

if [ -n "$cmd" ]; then
  echo "listIsajetTranslation.out and ./listIsajetTranslation.output differ"
  exit 1;
fi

