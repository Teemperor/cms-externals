#! /bin/sh
# tests/HepPDT/testReadIsajet.sh.  Generated from testReadIsajet.sh.in by configure.

rm -f testReadIsajet.out

./testReadIsajet
if ( ! `sed 's/e-0\([0-9][0-9]\)/e-\1/g' testReadIsajet.out | \
  sed 's/e+0\([0-9][0-9]\)/e+\1/g'  | \
  diff -q -b - ./testReadIsajet.output > /dev/null` )
then
  echo "testReadIsajet.out and ./testReadIsajet.output differ"
  exit 1;
fi

exit 0;
