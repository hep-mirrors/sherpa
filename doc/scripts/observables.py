"""The minimum-bias measurements used in the tune, one entry per histogram.

Every entry is a measurement in pp collisions that Rivet ships as reference data
(rivetbootstrap/local/share/Rivet/<ANALYSIS>.yoda.gz).  mb_analysis.py fills the same
quantity from Sherpa or Pythia events, with the same cuts as the Rivet analysis
(the .cc file named in 'analysis').  No deuteron measurement is used anywhere.

Keys of an entry:
  name       short name, used in file names and plots
  analysis   Rivet analysis (= the paper)
  histogram  histogram in the reference file
  energy     sqrt(s) [GeV]
  kind       what is measured:
               eta     1/N_ev dN_ch/deta                      (per selected event)
               pt      1/N_ev 1/(2 pi pT) d2N_ch/deta dpT      (per selected event, per unit eta)
               nch     1/N_ev dN_ev/dn_ch                     (fraction of selected events)
               ptnch   <pT> of the charged particles vs n_ch   (average over particles)
               ue_nch  transverse-region charged-particle density vs leading-track pT
               ue_ptsum transverse-region scalar pT-sum density vs leading-track pT
               ident   1/N_ev d2N/dy dpT of pi, K or p at |y| < 0.5 (per inelastic event)
  ptmin      lowest charged-particle pT [GeV]
  etamax     |eta| acceptance of the counted particles
  sel_etamax |eta| acceptance of the particles that select the event (= etamax, except ALICE 13 TeV)
  nchmin     lowest number of charged particles in the acceptance (event selection)
  no_strange_baryons  True: Sigma+-, Xi-, Omega- are not counted (ATLAS 13 TeV)
  tune       True: used in the fit.  False: shown in the plots only.
  reason     why it is (not) used
"""

# Diffraction: Sherpa's diffractive events (single and double dissociation) are mixed in
# with Sherpa's own cross sections (see combine.py).  Measurements that need many charged
# particles (n_ch >= 6 with pT > 500 MeV, n_ch >= 20, the underlying event with a leading
# track above 1 GeV) contain almost no diffractive events; the ones with n_ch >= 1 or 2
# do, and there the result depends on how well diffraction is modelled.  The tune uses
# the first kind only; the second kind is shown as a check.

ATLAS_MB_2010 = "ATLAS_2010_I882098"     # arXiv:1012.5104, charged particles at 0.9 and 7 TeV
ATLAS_UE_2010 = "ATLAS_2010_I879407"     # arXiv:1012.0791, underlying event at 0.9 and 7 TeV
ATLAS_MB_2016 = "ATLAS_2016_I1419652"    # arXiv:1602.01633, charged particles at 13 TeV
ATLAS_UE_2017 = "ATLAS_2017_I1509919"    # arXiv:1701.05390, underlying event at 13 TeV
ALICE_MB_2010 = "ALICE_2010_I852264"     # arXiv:1004.3514, dN/deta (INEL>0) at 0.9, 2.36, 7 TeV
ALICE_MB_2015 = "ALICE_2015_I1395253"    # arXiv:1509.08734, dN/deta and pT (INEL>0) at 13 TeV
ALICE_PIKP_2015 = "ALICE_2015_I1357424"  # arXiv:1504.00024, pi, K, p at 7 TeV (used for Lund_tune1)

DIFFRACTION_SUPPRESSED = "few diffractive events pass this selection"
DIFFRACTION_SENSITIVE = "depends on diffraction (n_ch >= 1 or 2): check only"
UNDERLYING_EVENT = "underlying event: activity next to the hardest track, set by the multiple scatterings"

OBSERVABLES = []


def add(name, analysis, histogram, energy, kind, ptmin, etamax, nchmin, tune, reason,
        no_strange_baryons=False, species=None, sel_etamax=None):
    entry = {}
    entry["name"] = name
    entry["analysis"] = analysis
    entry["histogram"] = histogram
    entry["energy"] = energy
    entry["kind"] = kind
    entry["ptmin"] = ptmin
    entry["etamax"] = etamax
    entry["nchmin"] = nchmin
    entry["tune"] = tune
    entry["reason"] = reason
    entry["no_strange_baryons"] = no_strange_baryons
    entry["species"] = species
    if sel_etamax is None:
        entry["sel_etamax"] = etamax
    else:
        entry["sel_etamax"] = sel_etamax
    OBSERVABLES.append(entry)


# ---------------------------------------------------------------------------------------
# 0.9 TeV
# ---------------------------------------------------------------------------------------
add("atlas900_nch6_eta", ATLAS_MB_2010, "d06-x01-y01", 900, "eta", 0.5, 2.5, 6, True, DIFFRACTION_SUPPRESSED)
add("atlas900_nch6_pt", ATLAS_MB_2010, "d13-x01-y01", 900, "pt", 0.5, 2.5, 6, True, DIFFRACTION_SUPPRESSED)
add("atlas900_nch6_nch", ATLAS_MB_2010, "d20-x01-y01", 900, "nch", 0.5, 2.5, 6, True, DIFFRACTION_SUPPRESSED)
add("atlas900_nch1_eta", ATLAS_MB_2010, "d01-x01-y01", 900, "eta", 0.5, 2.5, 1, False, DIFFRACTION_SENSITIVE)
add("atlas900_nch1_pt", ATLAS_MB_2010, "d08-x01-y01", 900, "pt", 0.5, 2.5, 1, False, DIFFRACTION_SENSITIVE)
add("atlas900_nch1_nch", ATLAS_MB_2010, "d15-x01-y01", 900, "nch", 0.5, 2.5, 1, False, DIFFRACTION_SENSITIVE)
add("atlas900_nch1_ptnch", ATLAS_MB_2010, "d22-x01-y01", 900, "ptnch", 0.5, 2.5, 1, True,
    "<pT> vs n_ch: each n_ch bin is a fixed multiplicity, so the diffractive fraction does not matter much")
add("atlas900_nch2_eta", ATLAS_MB_2010, "d04-x01-y01", 900, "eta", 0.1, 2.5, 2, False, DIFFRACTION_SENSITIVE)
add("atlas900_nch2_nch", ATLAS_MB_2010, "d18-x01-y01", 900, "nch", 0.1, 2.5, 2, False, DIFFRACTION_SENSITIVE)
add("atlas900_ue_nch", ATLAS_UE_2010, "d01-x01-y01", 900, "ue_nch", 0.5, 2.5, 1, True, UNDERLYING_EVENT)
add("atlas900_ue_ptsum", ATLAS_UE_2010, "d03-x01-y01", 900, "ue_ptsum", 0.5, 2.5, 1, True, UNDERLYING_EVENT)
add("alice900_inel0_eta", ALICE_MB_2010, "d04-x01-y01", 900, "eta", 0.0, 1.0, 1, False, DIFFRACTION_SENSITIVE)

# ---------------------------------------------------------------------------------------
# 7 TeV
# ---------------------------------------------------------------------------------------
add("atlas7000_nch6_eta", ATLAS_MB_2010, "d07-x01-y01", 7000, "eta", 0.5, 2.5, 6, True, DIFFRACTION_SUPPRESSED)
add("atlas7000_nch6_pt", ATLAS_MB_2010, "d14-x01-y01", 7000, "pt", 0.5, 2.5, 6, True, DIFFRACTION_SUPPRESSED)
add("atlas7000_nch6_nch", ATLAS_MB_2010, "d21-x01-y01", 7000, "nch", 0.5, 2.5, 6, True, DIFFRACTION_SUPPRESSED)
add("atlas7000_nch20_eta", ATLAS_MB_2010, "d27-x01-y01", 7000, "eta", 0.1, 2.5, 20, True, DIFFRACTION_SUPPRESSED)
add("atlas7000_nch20_nch", ATLAS_MB_2010, "d35-x01-y01", 7000, "nch", 0.1, 2.5, 20, True, DIFFRACTION_SUPPRESSED)
add("atlas7000_nch1_eta", ATLAS_MB_2010, "d03-x01-y01", 7000, "eta", 0.5, 2.5, 1, False, DIFFRACTION_SENSITIVE)
add("atlas7000_nch1_pt", ATLAS_MB_2010, "d10-x01-y01", 7000, "pt", 0.5, 2.5, 1, False, DIFFRACTION_SENSITIVE)
add("atlas7000_nch1_nch", ATLAS_MB_2010, "d17-x01-y01", 7000, "nch", 0.5, 2.5, 1, False, DIFFRACTION_SENSITIVE)
add("atlas7000_nch1_ptnch", ATLAS_MB_2010, "d23-x01-y01", 7000, "ptnch", 0.5, 2.5, 1, True,
    "<pT> vs n_ch: each n_ch bin is a fixed multiplicity, so the diffractive fraction does not matter much")
add("atlas7000_nch2_eta", ATLAS_MB_2010, "d05-x01-y01", 7000, "eta", 0.1, 2.5, 2, False, DIFFRACTION_SENSITIVE)
add("atlas7000_nch2_nch", ATLAS_MB_2010, "d19-x01-y01", 7000, "nch", 0.1, 2.5, 2, False, DIFFRACTION_SENSITIVE)
add("atlas7000_ue_nch", ATLAS_UE_2010, "d02-x01-y01", 7000, "ue_nch", 0.5, 2.5, 1, True, UNDERLYING_EVENT)
add("atlas7000_ue_ptsum", ATLAS_UE_2010, "d04-x01-y01", 7000, "ue_ptsum", 0.5, 2.5, 1, True, UNDERLYING_EVENT)
add("alice7000_inel0_eta", ALICE_MB_2010, "d06-x01-y01", 7000, "eta", 0.0, 1.0, 1, False, DIFFRACTION_SENSITIVE)
add("alice7000_inel0_nch", ALICE_MB_2010, "d03-x01-y01", 7000, "nch", 0.0, 1.0, 1, False, DIFFRACTION_SENSITIVE)
add("alice7000_pions", ALICE_PIKP_2015, "d01-x01-y01", 7000, "ident", 0.0, 0.0, 0, False,
    "per inelastic event: used for Lund_tune1, shown as a check", species=211)
add("alice7000_kaons", ALICE_PIKP_2015, "d01-x01-y02", 7000, "ident", 0.0, 0.0, 0, False,
    "per inelastic event: used for Lund_tune1, shown as a check", species=321)
add("alice7000_protons", ALICE_PIKP_2015, "d01-x01-y03", 7000, "ident", 0.0, 0.0, 0, False,
    "per inelastic event: used for Lund_tune1, shown as a check", species=2212)

# ---------------------------------------------------------------------------------------
# 13 TeV
# ---------------------------------------------------------------------------------------
add("atlas13000_ue_nch", ATLAS_UE_2017, "d21-x01-y01", 13000, "ue_nch", 0.5, 2.5, 1, True, UNDERLYING_EVENT,
    no_strange_baryons=True)
add("atlas13000_ue_ptsum", ATLAS_UE_2017, "d24-x01-y01", 13000, "ue_ptsum", 0.5, 2.5, 1, True, UNDERLYING_EVENT,
    no_strange_baryons=True)
add("atlas13000_nch1_eta", ATLAS_MB_2016, "d03-x01-y01", 13000, "eta", 0.5, 2.5, 1, False, DIFFRACTION_SENSITIVE,
    no_strange_baryons=True)
add("atlas13000_nch1_pt", ATLAS_MB_2016, "d04-x01-y01", 13000, "pt", 0.5, 2.5, 1, False, DIFFRACTION_SENSITIVE,
    no_strange_baryons=True)
add("atlas13000_nch1_nch", ATLAS_MB_2016, "d05-x01-y01", 13000, "nch", 0.5, 2.5, 1, False, DIFFRACTION_SENSITIVE,
    no_strange_baryons=True)
add("atlas13000_nch1_ptnch", ATLAS_MB_2016, "d06-x01-y01", 13000, "ptnch", 0.5, 2.5, 1, True,
    "<pT> vs n_ch: each n_ch bin is a fixed multiplicity, so the diffractive fraction does not matter much",
    no_strange_baryons=True)
add("alice13000_inel0_eta", ALICE_MB_2015, "d01-x01-y02", 13000, "eta", 0.0, 1.8, 1, False, DIFFRACTION_SENSITIVE,
    sel_etamax=1.0)


def observables_at(energy):
    "the entries measured at sqrt(s) = energy"
    result = []
    for entry in OBSERVABLES:
        if entry["energy"] == energy:
            result.append(entry)
    return result


ENERGIES = [900, 7000, 13000]
