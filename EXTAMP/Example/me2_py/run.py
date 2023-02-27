#!/usr/bin/env python

from subprocess import check_call
from filecmp import cmp as compare

check_call(["./test.py", "ME_SIGNAL_GENERATOR=Comix"])
check_call(["./test.py", "ME_SIGNAL_GENERATOR=External OpenLoops"])

from json import load

with open('./Comix') as infile:
    comix = load(infile)
    
with open('./External') as infile:
    external = load(infile)

err = abs((comix-external)/comix)
print("Relative difference: {0}".format(err))

exit(0 if abs < 1.0e-9 else 1)
    
