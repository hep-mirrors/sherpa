#!/usr/bin/env python

"""
# before running this script, create corrections.txt by running a simple setup like
# the following with
# $ Sherpa 2> corrections.txt

BEAMS: 2212
BEAM_ENERGIES: 6500
MI_HANDLER: None

HADRON_DECAYS:
  AlwaysIntegrate: true

PROCESSES:
- 93 93 -> 11 -12 93{0}:
    Order: {QCD: 0, EW: 2}
"""

import math
import ruamel.yaml
yaml = ruamel.yaml.YAML()
yaml.indent(mapping=2, sequence=4, offset=2)

decaydata = yaml.load(open("Decaydata.yaml"))

## # Example format for corrections.txt written by HADRONS++:
## # kfc channel xs dxs max
## 421 321,-321,-311,111 1.22179e-20 1.22173e-23 3.41777e-18
## 413 421,211 8.40687e-07 0 3.37993e-06
## 413 411,111 7.81938e-07 0 3.14373e-06
## 4122 2212,211,211,-211,-211 1.34865e-11 3.10543e-15 3.69694e-10
for line in open("corrections.txt"):
    words = line.split(" ")
    if not len(words)==5:
        continue
    oldmax = decaydata['HADRON_DECAYS']['Channels'][int(words[0])][words[1]]['IntResults'][2]
    newmax = float(words[4])

    if not math.isclose(oldmax, newmax, rel_tol=0.05):
        print(f"Updating {words[0]} -> {words[1]} \t\t(diff = {int(200*abs(newmax-oldmax)/(newmax+oldmax))} %)")
        decaydata['HADRON_DECAYS']['Channels'][int(words[0])][words[1]]['IntResults'][0]=float(words[2])
        decaydata['HADRON_DECAYS']['Channels'][int(words[0])][words[1]]['IntResults'][1]=float(words[3])
        decaydata['HADRON_DECAYS']['Channels'][int(words[0])][words[1]]['IntResults'][2]=float(words[4])

yaml.dump(decaydata,open("Decaydata.yaml","w"))
