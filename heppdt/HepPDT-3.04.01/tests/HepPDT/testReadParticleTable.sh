#! /bin/sh
# tests/HepPDT/testReadParticleTable.sh.  Generated from testReadParticleTable.sh.in by configure.

rm -f testReadParticleTable.out

./testReadParticleTable

if ( ! `sed 's/e-0\([0-9][0-9]\)/e-\1/g' testReadParticleTable.out | \
  sed 's/e+0\([0-9][0-9]\)/e+\1/g'  | \
  diff -q -b - ./testReadParticleTable.output > /dev/null` )
then
  echo "testReadParticleTable.out and ./testReadParticleTable.output differ"
  exit 1;
fi

if ( ! `sed 's/e-0\([0-9][0-9]\)/e-\1/g' testReadParticleTableStatus.out | \
  sed 's/e+0\([0-9][0-9]\)/e+\1/g'  | \
  diff -q -b - ./testReadParticleTableStatus.output > /dev/null` )
then
  echo "testReadParticleTableStatus.out and ./testReadParticleTableStatus.output differ"
  exit 1;
fi

exit 0;
