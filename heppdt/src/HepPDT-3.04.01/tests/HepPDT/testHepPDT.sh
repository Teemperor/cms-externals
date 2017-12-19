#! /bin/sh
# tests/HepPDT/testHepPDT.sh.  Generated from testHepPDT.sh.in by configure.

rm -f testHepPDT.out testHepPDTtable.out testHepPDTfragment.out testHepPDTstatus.out

./testHepPDT < testHepPDT.input 

if ( ! `sed 's/e-0\([0-9][0-9]\)/e-\1/g' testHepPDT.out | \
  sed 's/e+0\([0-9][0-9]\)/e+\1/g'  | \
  diff -q -b - ./testHepPDT.output > /dev/null` )
then
  echo "testHepPDT.out and ./testHepPDT.output differ"
  exit 1;
fi

if ( ! `sed 's/e-0\([0-9][0-9]\)/e-\1/g' testHepPDTtable.out | \
  sed 's/e+0\([0-9][0-9]\)/e+\1/g'  | \
  diff -q -b - ./testHepPDTtable.output > /dev/null` )
then
  echo "testHepPDTtable.out and ./testHepPDTtable.output differ"
  exit 1;
fi

if ( ! `sed 's/e-0\([0-9][0-9]\)/e-\1/g' testHepPDTfragment.out | \
  sed 's/e+0\([0-9][0-9]\)/e+\1/g'  | \
  diff -q -b - ./testHepPDTfragment.output > /dev/null` )
then
  echo "testHepPDTfragment.out and ./testHepPDTfragment.output differ"
  exit 1;
fi

if ( ! `sed 's/e-0\([0-9][0-9]\)/e-\1/g' testHepPDTstatus.out | \
  sed 's/e+0\([0-9][0-9]\)/e+\1/g'  | \
  diff -q -b - ./testHepPDTstatus.output > /dev/null` )
then
  echo "testHepPDTstatus.out and ./testHepPDTstatus.output differ"
  exit 1;
fi

exit 0;
