#!/bin/sh
nohup ./test-ipol.exe /tmp/ipol.unpolarized.hepmc      ipol.unpolarized.root      0 1000000 &
nohup ./test-ipol.exe /tmp/ipol.polarized.hepmc        ipol.polarized.root        1 1000000 &
nohup ./test-ipol.exe /tmp/ipol.modified.plzap0.hepmc  ipol.modified.plzap0.root  2 1000000 &
nohup ./test-ipol.exe /tmp/ipol.modified.costhe0.hepmc ipol.modified.costhe0.root 3 1000000 &
