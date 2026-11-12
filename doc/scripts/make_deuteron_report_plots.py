"""Runs the code of every notebook in SimpleLEP/XSecs (as written by make_notebooks.py) twice, once
for each proposed Sherpa setting, and saves all figures and printed tables into
doc/deuteron_report/figures/.

    python3 make_deuteron_report_plots.py            # the two AHADIC settings and the setting-independent checks
    python3 make_deuteron_report_plots.py --lund     # Sherpa with Lund strings: s6_008's activity -> figures/lund,
                                                     # re-tuned activity (Lund_tune1) -> figures/lund_tuned

Each pass shows the data, Pythia 8.318 and ONE Sherpa setting:
    figures/s6_008/<notebook>_<nn>.{png,pdf}   and   figures/s6_008/<notebook>.txt
    figures/s6_010/...
The .txt files hold everything the notebook prints (fitted 1/sigma_0, chi2, MC/data ratios).
Model_Validation and Coalescence_Models read fixed earlier runs (not the settings), so they
are run once, into figures/common/.  Colours as in deuteron_tune_report/deuteron_tune_report.pdf."""
import io, os, re, sys, contextlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paths import XSECS, DEUTERON_REPORT
os.chdir(XSECS)      # the notebooks' folder, as when they run
import deuteron_common as dc
import make_notebooks

PYTHIA = "Pythia 8.318"
SETTINGS = ["Sherpa s6_008", "Sherpa s6_010"]
COMMON = ["Model_Validation.ipynb", "Coalescence_Models.ipynb"]
dc.STYLE["Sherpa s6_008"] = dict(color="#D55E00", ls="-", lw=2.0)
dc.STYLE["Sherpa s6_010"] = dict(color="#009E73", ls="-", lw=2.0)
dc.STYLE[PYTHIA] = dict(color="#0072B2", ls="--", lw=1.6)
LUND = "--lund" in sys.argv
wanted = [a for a in sys.argv[1:] if a != "--lund"] or list(make_notebooks.NOTEBOOKS)

if LUND:
    passes = [("lund", ["Sherpa Lund", PYTHIA], [n for n in wanted if n not in COMMON]),
              ("lund_tuned", ["Sherpa Lund tuned", PYTHIA], [n for n in wanted if n not in COMMON]),
              ("lund_tune2", ["Sherpa Lund tune2", "Sherpa Lund tuned", PYTHIA], [n for n in wanted if n not in COMMON]),
              ("lund_as", ["Sherpa Lund FS_AS_FAC 0.5", "Sherpa Lund FS_AS_FAC 0.35", "Sherpa Lund tuned", PYTHIA],
               [n for n in wanted if n in ("LEP_ALEPH_OPAL.ipynb", "BaBar_ARGUS_CLEO.ipynb")])]
else:
    passes = [(setting.split()[-1], [setting, PYTHIA], [n for n in wanted if n not in COMMON]) for setting in SETTINGS]
    passes.append(("common", SETTINGS[:1] + [PYTHIA], [n for n in wanted if n in COMMON]))
for folder, shown, notebooks in passes:
    dc.SHOWN = shown
    dc.MODELS = {m: dc.ALL_MODELS[m] for m in shown}
    outdir = os.path.join(DEUTERON_REPORT, "figures", folder)
    os.makedirs(outdir, exist_ok=True)
    for name in notebooks:
        stem = name.replace(".ipynb", "")
        counter = [0]
        def show(*args, **kwargs):
            "replaces plt.show: write the current figure as <stem>_<nn>"
            counter[0] += 1
            fig = plt.gcf()
            # values read off the paper's figures (Herwig++, Pythia 8.186) are not shown
            for ax in fig.axes:
                # no "cf. Fig. N" (of the paper) in the titles
                title = re.sub(r"\s*\(cf\. [^)]*\)|cf\. Fig\. \d+:?\s*|,?\s*cf\. Fig\. \d+", "", ax.get_title())
                ax.set_title(title.strip(" ,"))
                paper = [l for l in ax.get_lines() if "(paper)" in l.get_label()]
                if not paper: continue
                for l in paper: l.remove()
                ax.relim(); ax.autoscale_view(); ax.set_ylim(0, None); ax.legend(fontsize=8)
            for ext in ("png", "pdf"):
                fig.savefig(os.path.join(outdir, f"{stem}_{counter[0]:02d}.{ext}"), dpi=130 if ext == "png" else None)
            plt.close(fig)
        plt.show = show
        namespace = {"__name__": "__notebook__"}
        printed = io.StringIO()
        with contextlib.redirect_stdout(printed):
            exec("from deuteron_common import *", namespace)
            for cell in make_notebooks.NOTEBOOKS[name]["cells"]:
                if cell["cell_type"] != "code": continue
                lines = []
                for line in cell["source"]:
                    if line.startswith("%"): continue
                    if line.startswith("COMPARE = ["): line = f"COMPARE = {shown!r}\n"
                    lines.append(line)
                try:
                    exec("".join(lines), namespace)
                except Exception as error:
                    print(f"\n!! cell failed: {type(error).__name__}: {error}\n")
                    plt.close("all")
        open(os.path.join(outdir, stem + ".txt"), "w").write(printed.getvalue())
        print(f"{folder}: {stem}: {counter[0]} figures", file=sys.stderr)
