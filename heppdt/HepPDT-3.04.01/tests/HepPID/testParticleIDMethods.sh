#! /bin/sh
# tests/HepPID/testParticleIDMethods.sh.  Generated from testParticleIDMethods.sh.in by configure.

./testParticleIDMethods 

cmd=`diff -q -b testParticleIDMethods.out ./testParticleIDMethods.output`

if [ -n "$cmd" ]; then
  echo "testParticleIDMethods.out and ./testParticleIDMethods.output differ"
  exit 1;
fi
