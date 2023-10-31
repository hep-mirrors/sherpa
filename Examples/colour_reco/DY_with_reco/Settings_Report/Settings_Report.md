---
title: Sherpa run-time settings
date: Thu Sep 28 17:19:15 2023
...

Unused settings
-------------------
Parameters that have never been read by Sherpa during its run are listed here. If you did expect the setting to be used, check its spelling, and note that Sherpa setting names are case-sensitive.

### DY.yml
- INTRINSIC_KPERP:2212:FORM
- INTRINSIC_KPERP:2212:RECOIL
- INTRINSIC_KPERP:2212:REFERENCE_ENERGY
- INTRINSIC_KPERP:2212:ENERGY_SCALING_EXPO
- INTRINSIC_KPERP:2212:SHOWER_INITIATOR_MEAN
- INTRINSIC_KPERP:2212:SHOWER_INITIATOR_SIGMA
- INTRINSIC_KPERP:2212:SHOWER_INITIATOR_Q2
- INTRINSIC_KPERP:2212:SHOWER_INITIATOR_KTMAX
- INTRINSIC_KPERP:2212:SHOWER_INITIATOR_KTEXPO
- INTRINSIC_KPERP:2212:BEAM_SPECTATOR_MEAN
- INTRINSIC_KPERP:2212:BEAM_SPECTATOR_SIGMA
- INTRINSIC_KPERP:2212:BEAM_SPECTATOR_Q2
- INTRINSIC_KPERP:2212:BEAM_SPECTATOR_KTMAX
- INTRINSIC_KPERP:2212:BEAM_SPECTATOR_KTEXPO
- CSS_FORCED_DECAYS
- CSS_FORCED_GLUON_SCALING


Customised settings
-------------------
The parameters listed here have been customised in some way and they have been read by Sherpa during its run. The last column lists the actual value used after taking into account all setting sources (default values, overrides, input files and the command line).

In some cases, an alternative default value is being used. These alternatives will be separated by "`-- AND --`" from the standard default, which will always be listed on top.

Note that parameters that can take on different values because they are set within a list, for example `param: [{x: 1}, {x: 2}, ...]`, will not appear in the config-file or command-line columns. They will be listed in the final-value column, with each different value separated by an "`-- AND --`" line.

| parameter | default value | override by SHERPA | /Applications/sherpa/sherpa/build/share/SHERPA-MC//Decaydata.yaml | command line | final value |
|-|-|-|-|-|-|
ALPHAS:USE\_PDF| 1 |  |  | 0 |  | 0 |  |
ALPHAS\(MZ\)| 0\.118 |  |  | 0\.1188 |  | 0\.1188 |  |
ANALYSIS|  |  |  | Rivet |  | Rivet |  |
BEAMS| 0, 0 |  |  | 2212 |  | 2212 |  |
BEAM\_ENERGIES| 0, 0 |  |  | 3500 |  | 3500 |  |
COLOUR\_RECONNECTIONS:MODE| 0 |  |  | true |  | 1 |  |
FRAGMENTATION| Ahadic |  |  | None |  | None |  |
ME\_GENERATORS| Comix, Amegic, Internal |  |  | Comix, Amegic |  | Comix, Amegic |  |
ORDER\_ALPHAS| 2 |  |  | 1 |  | 1 |  |
PROCESSES:5 \-5 \-\> 13 \-13:Order:EW| \-1 |  |  |  |  | 2 |  |
PROCESSES:5 \-5 \-\> 13 \-13:Order:QCD| \-1 |  |  |  |  | 0 |  |
RIVET:\-\-analyses|  |  |  | ATLAS\_2011\_S9131140 |  | ATLAS\_2011\_S9131140 |  |
RUNDATA|  |  |  |  | DY\.yml | DY\.yml |  |
SCALES| METS\{MU\_F2\}\{MU\_R2\}\{MU\_Q2\} |  |  | METS\{H\_T2\}\{H\_T2\} |  | METS\{H\_T2\}\{H\_T2\} |  |
SELECTORS|  |  |  | Mass, 13, \-13, 66, E\_CMS |  | Mass, 13, \-13, 66, 7000 |  |
Settings kept at their default value
-------------------
The parameter listed here have not been customised, but they have been read by Sherpa during its run.

| parameter | default value |
|-|-|
1/ALPHAQED\(0\)| 137\.03599976 |  |
ABS\_ERROR| 0 |  |
ALPHAS:FREEZE\_VALUE| 1 |  |
ALPHAS:PDF\_SET\_VERSION| 0 |  |
AMEGIC:DEFAULT\_GAUGE| 1 |  |
AMEGIC:PARTIAL\_COMMIT| 0 |  |
AMISIC:E\(ref\)| 7000 |  |
AMISIC:Eta| 0\.08 |  |
AMISIC:MATTER\_FORM| Single\_Gaussian |  |
AMISIC:MATTER\_FRACTION1| 0\.5 |  |
AMISIC:MATTER\_RADIUS1| 1 |  |
AMISIC:MATTER\_RADIUS2| 2 |  |
AMISIC:MU\_F\_FACTOR| 1 |  |
AMISIC:MU\_R\_FACTOR| 0\.5 |  |
AMISIC:MU\_R\_SCHEME| PT |  |
AMISIC:PT\_0| 2\.2 |  |
AMISIC:PT\_0\(IR\)| 0\.5 |  |
AMISIC:PT\_0\(ref\)| 2\.2 |  |
AMISIC:PT\_Min| 2\.5 |  |
AMISIC:PT\_Min\(ref\)| 2\.5 |  |
AMISIC:PomeronIntercept| 0\.0808 |  |
AMISIC:PomeronSlope| 0\.25 |  |
AMISIC:ReggeonIntercept| \-0\.4525 |  |
AMISIC:SIGMA\_ND\_NORM| 1 |  |
AMISIC:TriplePomeronCoupling| 0\.318 |  |
AMISIC:nMC\_points| 1000 |  |
AMISIC:nPT\_bins| 200 |  |
AMISIC:nS\_bins| 100 |  |
ANALYSIS\_OUTPUT| Analysis/ |  |
ANALYSIS\_WRITEOUT\_INTERVAL| 1\.84467440737e\+19 |  |
ASSOCIATED\_CONTRIBUTIONS\_VARIATIONS|  |  |
AS\_FORM| Smooth |  |
BATCH\_MODE| 1 |  |
BBR\_PDF\_LIBRARY|  |  |
BEAM\_1| 0 |  |
BEAM\_2| 0 |  |
BEAM\_ENERGY\_1| 0 |  |
BEAM\_ENERGY\_2| 0 |  |
BEAM\_MODE| Collider |  |
BEAM\_POLARIZATIONS| 0, 0 |  |
BEAM\_REMNANTS| 1 |  |
BEAM\_RESCATTERING| None |  |
BEAM\_SPECTRA| Monochromatic, Monochromatic |  |
BEAM\_SPECTRUM\_1| Monochromatic |  |
BEAM\_SPECTRUM\_2| Monochromatic |  |
BRH\_VMODE| 0 |  |
BUNCHES|  |  |
CHECK\_LIBLOCK| 0 |  |
CHECK\_SETTINGS| 1 |  |
CHECK\_WEIGHT| 0 |  |
CITATION\_DEPTH| 1 |  |
CI\_OMODE| 1 |  |
CKM:Order| 0 |  |
CKM:Output| 0 |  |
COLOUR\_RECONNECTIONS:PMODE| 0 |  |
COLOUR\_RECONNECTIONS:Q\_0| 1 |  |
COLOUR\_RECONNECTIONS:R\_0| 100 |  |
COLOUR\_RECONNECTIONS:Reshuffle| 0\.111111111111 |  |
COLOUR\_RECONNECTIONS:etaQ| 0\.1 |  |
COLOUR\_RECONNECTIONS:etaR| 0\.16 |  |
COLOUR\_RECONNECTIONS:kappa| 2 |  |
COLOUR\_SCHEME| 0 |  |
COMIX:AEXP| 0\.9 |  |
COMIX:BMODE| 1 |  |
COMIX:CHECK\_POLES| 0 |  |
COMIX:ECMODE| 2 |  |
COMIX:ITMAX| 1000000 |  |
COMIX:ITMIN| 1000 |  |
COMIX:MFAC| 1 |  |
COMIX:N\_GPL| 3 |  |
COMIX:OMODE| 3 |  |
COMIX:PARTIAL\_COMMIT| 0 |  |
COMIX:PG\_MODE| 0 |  |
COMIX:PMODE| D |  |
COMIX:PS\_ADD\_ANISO| 0 |  |
COMIX:PS\_CHTH| 0\.01 |  |
COMIX:SEXP| 1 |  |
COMIX:SPEAK| 1 |  |
COMIX:SRBASE| 1\.05 |  |
COMIX:STEXP| 0\.001 |  |
COMIX:TEXP| 0\.9 |  |
COMIX:THEXP| 1\.5 |  |
COMIX:TMODE| 1 |  |
COMIX:VINTS| 8 |  |
COMIX:VL\_MODE| 0 |  |
COMIX:VMODE| 1 |  |
COMIX:VSOPT| 1 |  |
COMIX:WF\_MODE| 0 |  |
COMIX:ZMODE| 0 |  |
COMIX\_DEFAULT\_GAUGE| 1 |  |
COMPRESS\_PARTONIC\_DECAYS| 1 |  |
COUPLINGS| Alpha\_QCD 1 |  |
CSS\_CKFMODE| 1 |  |
CSS\_ENHANCE|  |  |
CSS\_EVOLUTION\_SCHEME| 3030 |  |
CSS\_EW\_MODE| 0 |  |
CSS\_FS\_AS\_FAC| 1 |  |
CSS\_FS\_PT2MIN| 1 |  |
CSS\_IS\_AS\_FAC| 0\.5 |  |
CSS\_IS\_PT2MIN| 2 |  |
CSS\_KFACTOR\_SCHEME| 1 |  |
CSS\_KIN\_SCHEME| 0 |  |
CSS\_KMODE| 2 |  |
CSS\_MASS\_THRESHOLD| 0 |  |
CSS\_MAXEM| 18446744073709551615 |  |
CSS\_MAXPART| 2147483647 |  |
CSS\_MAX\_REWEIGHT\_FACTOR| 1000 |  |
CSS\_PDFCHECK| 1 |  |
CSS\_PDF\_FAC| 1 |  |
CSS\_PDF\_MIN| 0\.0001 |  |
CSS\_PDF\_MIN\_X| 0\.01 |  |
CSS\_QCD\_MODE| 1 |  |
CSS\_RECO\_CHECK| 0 |  |
CSS\_RECO\_DECAYS| 0 |  |
CSS\_RESPECT\_Q2| 0 |  |
CSS\_REWEIGHT| 1 |  |
CSS\_REWEIGHT\_SCALE\_CUTOFF| 5 |  |
CSS\_SCALE\_FACTOR| 1 |  |
CSS\_SCALE\_SCHEME| 0 |  |
CSS\_SCALE\_VARIATION\_SCHEME| 1 |  |
DEBUG\_INTERVAL| 0 |  |
DEBUG\_STEP| \-1 |  |
DECAYER| 0 |  |
DECOMPOSE\_4G\_VERTEX| 1 |  |
DIPOLES:ALPHA| 1 |  |
DIPOLES:ALPHA\_FF| 1 |  |
DIPOLES:ALPHA\_FI| 1 |  |
DIPOLES:ALPHA\_IF| 1 |  |
DIPOLES:ALPHA\_II| 1 |  |
DIPOLES:AMIN| 1e\-08 |  |
DIPOLES:KAPPA| 0\.666666666667 |  |
DIPOLES:KT2MAX| 49000000 |  |
DIPOLES:NF\_GSPLIT| 5 |  |
DIPOLES:SCHEME| CSS |  |
DM\_ENERGY\_DISTRIBUTION| 1 |  |
DM\_RELATIVISTIC| 1 |  |
DM\_TEMPERATURE| 1 |  |
DM\_beam\_weighted| 1 |  |
ENHANCE\_XS| 0 |  |
ERROR| 0\.01 |  |
EVENTS| 100 |  |
EVENT\_DISPLAY\_INTERVAL| 100 |  |
EVENT\_GENERATION\_MODE| PartiallyUnweighted |  |
EVENT\_INPUT|  |  |
EVENT\_OUTPUT|  |  |
EVENT\_SEED\_FILE| ran\.stat\.1234 |  |
EVENT\_SEED\_MODE| 0 |  |
EVENT\_TYPE| StandardPerturbative |  |
EVT\_FILE\_PATH| \. |  |
EVT\_OUTPUT| 2 |  |
EVT\_OUTPUT\_START| 0 |  |
EW\_REN\_SCHEME| Gmu |  |
EW\_SCHEME| Gmu |  |
EXTERNAL\_RNG| None |  |
EXTRAXS\_CSS\_APPROX\_ME| 0 |  |
FACTORIZATION\_SCALE\_FACTOR| 1 |  |
FINISH\_OPTIMIZATION| 1 |  |
FLAG\_PARTONIC\_DECAYS| 1 |  |
FREEZE\_PDF\_FOR\_LOW\_Q| 0 |  |
GENERATE\_RESULT\_DIRECTORY| 1 |  |
GF| 1\.16639e\-05 |  |
GLOBAL\_KFAC| 0 |  |
GMU\_CMS\_AQED\_CONVENTION| 0 |  |
HADRON\_DECAYS:Spin\_Correlations| 0 |  |
HARD\_DECAYS:Enabled| 0 |  |
HARD\_DECAYS:Spin\_Correlations| 0 |  |
HELICITY\_SCHEME| 1 |  |
HEPMC\_EXTENDED\_WEIGHTS| 0 |  |
HEPMC\_TREE\_LIKE| 0 |  |
HEPMC\_USE\_NAMED\_WEIGHTS| 1 |  |
HISTOGRAM\_OUTPUT\_PRECISION| 6 |  |
IB\_THRESHOLD\_KILL| \-1e\+12 |  |
IB\_WHBINS| 100 |  |
INIT\_ONLY| 0 |  |
INTEGRATION\_ERROR| 0\.01 |  |
INTEGRATOR| Default |  |
INTRINSIC\_KPERP:CUT\_EXPO| 5 |  |
INTRINSIC\_KPERP:FORM| gauss\_limited |  |
INTRINSIC\_KPERP:MAX| 3 |  |
INTRINSIC\_KPERP:MEAN| 0 |  |
INTRINSIC\_KPERP:Q2| 0\.77 |  |
INTRINSIC\_KPERP:REFE| 7000 |  |
INTRINSIC\_KPERP:SCALE\_EXPO| 0\.08 |  |
INTRINSIC\_KPERP:SIGMA| 1\.5 |  |
INT\_MINSIJ\_FACTOR| 0 |  |
JET\_MASS\_THRESHOLD| 10 |  |
KFACTOR| None |  |
KFACTOR\_ALLOW\_MAPPING| 1 |  |
LHAPDF:DISALLOW\_FLAVOUR|  |  |
LHAPDF:USE\_Q2LIMIT| 1 |  |
LHEF\_PDF\_NUMBER| \-1 |  |
LOG\_FILE|  |  |
MASSIVE\_PS|  |  |
MASSLESS\_PS|  |  |
MCNLO\_DADS| 1 |  |
MEH\_EWADDMODE| 0 |  |
MEH\_NLOADD| 1 |  |
MEH\_QCDADDMODE| 0 |  |
MEMLEAK\_WARNING\_THRESHOLD| 16777216 |  |
MENLOPS\_MAX\_KFAC| 10 |  |
MEPS:ALLOW\_SCALE\_UNORDERING| 0 |  |
MEPS:CLUSTER\_MODE| 360 |  |
MEPS:CORE\_SCALE| Default |  |
MEPS:MEPS\_COLORSET\_MODE| 0 |  |
MEPS:NLO\_COUPLING\_MODE| 2 |  |
MEPS:NLO\_NMAX\_ALLCONFIGS| \-1 |  |
MEPS:NMAX\_ALLCONFIGS| \-1 |  |
MEPS:UNORDERED\_SCALE| None |  |
MEPSNLO\_PDFCT| 1 |  |
METS:CLUSTER\_MODE| 0 |  |
METS\_BBAR\_MODE| EnabledExclCluster |  |
ME\_QED:CLUSTERING\_ENABLED| 1 |  |
ME\_QED:CLUSTERING\_THRESHOLD| 10 |  |
ME\_QED:ENABLED| 1 |  |
ME\_QED:INCLUDE\_RESONANCES| 0 |  |
MI\_CSS\_FS\_AS\_FAC| 0\.66 |  |
MI\_CSS\_FS\_PT2MIN| 1 |  |
MI\_CSS\_IS\_AS\_FAC| 0\.66 |  |
MI\_CSS\_IS\_PT2MIN| 4 |  |
MI\_CSS\_KFACTOR\_SCHEME| 0 |  |
MI\_CSS\_KIN\_SCHEME| 1 |  |
MI\_HANDLER| Amisic |  |
MODEL| SM |  |
MPI\_OUTPUT| 0 |  |
MPI\_PDF\_LIBRARY|  |  |
MPI\_PDF\_SET|  |  |
MPI\_PDF\_SET\_VERSIONS|  |  |
MPI\_PT\_MAX| 1e\+12 |  |
MPI\_PT\_Max\_Fac| 1 |  |
NLO\_IMODE| IKP |  |
NLO\_MUR\_COEFFICIENT\_FROM\_VIRTUAL| 1 |  |
NLO\_NF\_CONVERSION\_TERMS| None |  |
NLO\_SMEAR\_POWER| 0\.5 |  |
NLO\_SMEAR\_THRESHOLD| 0 |  |
NLO\_SUBTRACTION\_MODE| QCD |  |
NO\_ZERO\_PDF| 0 |  |
NUM\_ACCURACY| 1e\-10 |  |
OUTPUT| 2 |  |
OUTPUT\_ME\_ONLY\_VARIATIONS| 1 |  |
OVERWEIGHT\_THRESHOLD| 1e\+12 |  |
PB\_USE\_FMM| 0 |  |
PDF\_LIBRARY|  |  |
PDF\_SET|  |  |
PDF\_SET\_VERSIONS|  |  |
PRETTY\_PRINT| On |  |
PRINT\_PS\_POINTS| 0 |  |
PRINT\_VERSION\_INFO| 0 |  |
PROCESSES:5 \-5 \-\> 13 \-13:CKKW|  |  |
PROCESSES:5 \-5 \-\> 13 \-13:Cut\_Core| 0 |  |
PROCESSES:5 \-5 \-\> 13 \-13:Decay|  |  |
PROCESSES:5 \-5 \-\> 13 \-13:DecayOS|  |  |
PROCESSES:5 \-5 \-\> 13 \-13:No\_Decay|  |  |
PROCESSES:5 \-5 \-\> 13 \-13:Sort\_Flavors| 3 |  |
PS\_PT\_FILE|  |  |
Q2\_AS| 1 |  |
RANDOM\_SEED| \-1, \-1, \-1, \-1 |  |
RANDOM\_SEED1| \-1 |  |
RANDOM\_SEED2| \-1 |  |
RANDOM\_SEED3| \-1 |  |
RANDOM\_SEED4| \-1 |  |
RELIC\_DENSITY\_EMAX| 1000000 |  |
REMNANTS:DELTA\_MASS| 1\.5 |  |
REMNANTS:SOFT\_ETA\_RANGE| 7\.5 |  |
REMNANTS:SOFT\_MASS| 5 |  |
REMNANTS:SOFT\_X\_EXPONENT| \-2 |  |
RENORMALIZATION\_SCALE\_FACTOR| 1 |  |
RESPECT\_MASSIVE\_FLAG| 0 |  |
RESULT\_DIRECTORY| Results |  |
RESUMMATION\_SCALE\_FACTOR| 1 |  |
REWEIGHT\_SPLITTING\_ALPHAS\_SCALES| 1 |  |
REWEIGHT\_SPLITTING\_PDF\_SCALES| 1 |  |
RIVET:\-l| 20 |  |
RIVET:HISTO\_INTERVAL| 0 |  |
RIVET:IGNORE\_BEAMS| 0 |  |
RIVET:JETCONTS| 0 |  |
RIVET:MATCH\_WEIGHTS|  |  |
RIVET:NLO\_SMEARING| 0 |  |
RIVET:NOMINAL\_WEIGHT|  |  |
RIVET:SKIP\_WEIGHTS| 0 |  |
RIVET:SPLITCOREPROCS| 0 |  |
RIVET:SPLITPM| 0 |  |
RIVET:SPLITSH| 0 |  |
RIVET:UNMATCH\_WEIGHTS| ^EXTRA\_\_\.\*,^IRREG\_\_\.\* |  |
RIVET:USE\_HEPMC\_EXTENDED\_WEIGHTS| 0 |  |
RIVET:USE\_HEPMC\_NAMED\_WEIGHTS| 1 |  |
RIVET:USE\_HEPMC\_SHORT| 0 |  |
RIVET:USE\_HEPMC\_TREE\_LIKE| 0 |  |
RIVET:WEIGHT\_CAP| 0 |  |
RLIMIT\_AS| 2042626048 |  |
RLIMIT\_BY\_CPU| 0 |  |
RUN\_MASS\_BELOW\_POLE| 0 |  |
SAVE\_STATUS|  |  |
SCALE\_FACTOR| 1 |  |
SELECTION\_WEIGHT\_MODE| 0 |  |
SHERPA\_CPP\_PATH|  |  |
SHERPA\_LDADD|  |  |
SHERPA\_LIB\_PATH|  |  |
SHERPA\_VERSION|  |  |
SHOWER\_GENERATOR| CSS |  |
SHOW\_ANALYSIS\_SYNTAX| 0 |  |
SHOW\_FILTER\_SYNTAX| 0 |  |
SHOW\_KFACTOR\_SYNTAX| 0 |  |
SHOW\_ME\_GENERATORS| 0 |  |
SHOW\_MODEL\_SYNTAX| 0 |  |
SHOW\_NLOMC\_GENERATORS| 0 |  |
SHOW\_NTRIALS| 0 |  |
SHOW\_PDF\_SETS| 0 |  |
SHOW\_PS\_GENERATORS| 0 |  |
SHOW\_SCALE\_SYNTAX| 0 |  |
SHOW\_SELECTOR\_SYNTAX| 0 |  |
SHOW\_SHOWER\_GENERATORS| 0 |  |
SHOW\_VARIABLE\_SYNTAX| 0 |  |
SOFT\_COLLISIONS| None |  |
SP:ADD\_DOC| 0 |  |
SP:SET\_COLORS| 0 |  |
STATUS\_PATH|  |  |
THRESHOLD\_ALPHAS| 1 |  |
TIMEOUT| \-1 |  |
USERHOOKS| None |  |
USR\_WGT\_MODE| 1 |  |
VARIATIONS\_INCLUDE\_CV| 0 |  |
VEGAS\_MODE| 2 |  |
VIRTUAL\_EVALUATION\_FRACTION| 1 |  |
WIDTH\_SCHEME| CMS |  |
WRITE\_REFERENCES\_FILE| 1 |  |
YFS:1/ALPHAQED| 0 |  |
YFS:CHECK\_FIRST| 0 |  |
YFS:DRCUT| 1\.79769313486e\+308 |  |
YFS:FF\_RECOIL\_SCHEME| 2 |  |
YFS:FI\_RECOIL\_SCHEME| 2 |  |
YFS:INCREASE\_MAXIMUM\_WEIGHT| 1 |  |
YFS:IR\_CUTOFF| 0\.001 |  |
YFS:IR\_CUTOFF\_FRAME| Multipole\_CMS |  |
YFS:MAXEM| 2147483647 |  |
YFS:MINEM| 0 |  |
YFS:MODE| Full |  |
YFS:PHOTON\_SPLITTER\_ENHANCE\_FACTOR| 1 |  |
YFS:PHOTON\_SPLITTER\_MAX\_HADMASS| 0\.5 |  |
YFS:PHOTON\_SPLITTER\_MODE| 15 |  |
YFS:PHOTON\_SPLITTER\_ORDERING\_SCHEME| 2 |  |
YFS:PHOTON\_SPLITTER\_SPECTATOR\_SCHEME| 0 |  |
YFS:PHOTON\_SPLITTER\_STARTING\_SCALE\_SCHEME| 1 |  |
YFS:REDUCE\_MAXIMUM\_ENERGY| 1 |  |
YFS:STRICTNESS| 0 |  |
YFS:USE\_ME| 1 |  |
YFS:USE\_RUNNING\_PARAMETERS| 0 |  |
YFS:UV\_CUTOFF| 1\.79769313486e\+308 |  |
YUKAWA\_MASSES| Running |  |
