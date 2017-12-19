#! /bin/sh
# tests/HepPDT/testReadEvtGen.sh.  Generated from testReadEvtGen.sh.in by configure.

rm -f testReadEvtGen.out

./testReadEvtGen

if ( ! `sed 's/e-0\([0-9][0-9]\)/e-\1/g' testReadEvtGen.out | \
  sed 's/e+0\([0-9][0-9]\)/e+\1/g'  | \
  diff -q -b - ./testReadEvtGen.output > /dev/null` )
then
  echo "testReadEvtGen.out and ./testReadEvtGen.output differ"
  exit 1;
fi

exit 0;

