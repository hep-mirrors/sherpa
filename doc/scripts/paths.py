"""The folders the scripts in doc/scripts read from and write to.

Reports (each with its .tex, figures/ and the run cards it used in runcards/):
    lund_tune/                 (anti)deuterons with Sherpa + Lund strings: tunedLund.tex
    lund_tuning_report/        how Sherpa with Lund strings was tuned: lund_tuning_report.tex
    deuteron_report/           the first deuteron report and the plots of all settings
    deuteron_tune_report/      the ScanTune settings s6_008 and s6_010: deuteron_tune_report.tex
    deuteron_review/           review of the deuteron data (ALICE, BaBar, H1)
Tunes and runs:
    minimum_bias_tuning/                         the minimum-bias tune: cards/, data/, runs/, fit results
    minimum_bias_tuning/non_diffractive_events/  the tune with non-diffractive events only: fits, figures/
    examples/RemnantTests/SimpleLEP/             the deuteron runs (DeuteronFormationFromXSecs/), the
                                                 notebooks (XSecs/), the tune scans (TuneScan/), LEP runs
                                                 (LEPRivet/) and the Pythia reference (PythiaReference/)
"""
import os

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
DOC = os.path.dirname(SCRIPTS)
SHERPA_LUND = os.path.join(os.path.dirname(DOC), "build_lund", "bin", "Sherpa")   # Sherpa built with Pythia 8

# reports
LUND_TUNE = os.path.join(DOC, "lund_tune")
LUND_TUNING_REPORT = os.path.join(DOC, "lund_tuning_report")
DEUTERON_REPORT = os.path.join(DOC, "deuteron_report")
DEUTERON_TUNE_REPORT = os.path.join(DOC, "deuteron_tune_report")

# minimum-bias tune
TUNING = os.path.join(DOC, "minimum_bias_tuning")
NON_DIFFRACTIVE = os.path.join(TUNING, "non_diffractive_events")

# deuteron runs, notebooks and scans
SIMPLE_LEP = os.path.join(DOC, "examples", "RemnantTests", "SimpleLEP")
DEUTERON_DATA = os.path.join(SIMPLE_LEP, "DeuteronFormationFromXSecs")
XSECS = os.path.join(SIMPLE_LEP, "XSecs")
TUNE_SCAN = os.path.join(SIMPLE_LEP, "TuneScan")
APPRENTICE = os.path.join(TUNE_SCAN, "apprentice")
LEP_RIVET = os.path.join(SIMPLE_LEP, "LEPRivet")
LUND_RUNS = os.path.join(SIMPLE_LEP, "LundRuns")
PYTHIA_REFERENCE = os.path.join(SIMPLE_LEP, "PythiaReference")

# the report folder of every set of run cards: <report>/runcards/<set>
CARD_SET_REPORT = {}
CARD_SET_REPORT["Lund_tune1"] = LUND_TUNE
CARD_SET_REPORT["Lund_tune4"] = LUND_TUNE
CARD_SET_REPORT["Lund_mbtune"] = LUND_TUNE
CARD_SET_REPORT["Lund_tune2"] = LUND_TUNING_REPORT
CARD_SET_REPORT["Lund_tune3"] = LUND_TUNING_REPORT
CARD_SET_REPORT["Lund_as0.35"] = LUND_TUNING_REPORT
CARD_SET_REPORT["Lund_as0.5"] = LUND_TUNING_REPORT
CARD_SET_REPORT["Lund_s6_008"] = LUND_TUNING_REPORT
CARD_SET_REPORT["ScanTune_s6_008"] = DEUTERON_TUNE_REPORT
CARD_SET_REPORT["ScanTune_s6_010"] = DEUTERON_TUNE_REPORT
CARD_SET_REPORT["ScanTune_s3_053"] = DEUTERON_TUNE_REPORT
CARD_SET_REPORT["ScanTune_s4_005"] = DEUTERON_TUNE_REPORT
CARD_SET_REPORT["scan_dopt"] = DEUTERON_TUNE_REPORT
CARD_SET_REPORT["master_CRon"] = DEUTERON_REPORT
CARD_SET_REPORT["master_ourTune"] = DEUTERON_REPORT
CARD_SET_REPORT["retuneQ7"] = DEUTERON_REPORT
CARD_SET_REPORT["old_branch"] = DEUTERON_REPORT


def card_set(name):
    """The folder of the run-card set <name> (e.g. Lund_tune1).  A new set that is not in
    CARD_SET_REPORT yet goes to lund_tune/runcards/."""
    if name in CARD_SET_REPORT:
        report = CARD_SET_REPORT[name]
    else:
        report = LUND_TUNE
    return os.path.join(report, "runcards", name)
