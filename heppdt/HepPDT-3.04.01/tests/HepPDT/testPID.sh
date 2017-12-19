#! /bin/sh
# tests/HepPDT/testPID.sh.  Generated from testPID.sh.in by configure.

rm -f testPID.out

./testPID  >& testPID.out
if ( !  `diff -q -b testPID.out ./testPID.output > /dev/null` )
then
  echo "testPID.out and ./testPID.output differ"
  exit 1;
fi

exit 0;
