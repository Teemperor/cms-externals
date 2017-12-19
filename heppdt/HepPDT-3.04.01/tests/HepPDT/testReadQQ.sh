#! /bin/sh
# tests/HepPDT/testReadQQ.sh.  Generated from testReadQQ.sh.in by configure.

rm -f testReadQQ.out

./testReadQQ

if ( ! `sed 's/e-0\([0-9][0-9]\)/e-\1/g' testReadQQ.out | sed 's/e+0\([0-9][0-9]\)/e+\1/g' | diff -q -b - ./testReadQQ.output > /dev/null` )
then
  echo "testReadQQ.out and ./testReadQQ.output differ"
  exit 1;
fi

exit 0;
