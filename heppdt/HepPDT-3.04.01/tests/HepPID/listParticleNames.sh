#! /bin/sh
# tests/HepPID/listParticleNames.sh.  Generated from listParticleNames.sh.in by configure.

rm -f listParticleNames.out

./listParticleNames 

cmd=`diff -q -b listParticleNames.out ./listParticleNames.output`

if [ -n "$cmd" ]; then
  echo "listParticleNames.out and ./listParticleNames.output differ"
  exit 1;
fi

