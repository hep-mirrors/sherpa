#include "ATOOLS/Phys/Variations.H"

#include <iterator>
#include <numeric>
#include <algorithm>

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Strings.H"
#include "ATOOLS/Phys/Blob.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/PDF_Base.H"
#ifdef USING__LHAPDF
  #include "LHAPDF/LHAPDF.h"
#endif

namespace ATOOLS {

  Variations* s_variations;

  // utitilies for parsing variations
  struct ScaleFactor_Pair : public std::pair<double, double> {
    ScaleFactor_Pair(const double& a, const double& b)
        : std::pair<double, double>(a, b)
    {
    }
  };
  struct ExpandableVariation {
    ExpandableVariation(std::string raw_var)
    {
      if (raw_var.back() == '*') {
        raw_var = raw_var.substr(0, raw_var.size() - 1);
        expand = true;
      }
      var = raw_var;
    }
    std::string var;
    bool expand {false};
  };
}

using namespace ATOOLS;

bool Variations::NeedsLHAPDF6Interface()
{
  // go through all parameters and return true if any PDF variation is
  // requested
  Settings& s = Settings::GetMainSettings();
  for (auto single_variation_settings : s["VARIATIONS"].GetItems()) {
    if (single_variation_settings["PDF"].IsSetExplicitly()) {
      return true;
    }
  }
  for (auto single_variation_settings : s["PDF_VARIATIONS"].GetItems()) {
    if (single_variation_settings.Get<std::string>() != "None") {
      return true;
    }
  }
  return false;
}

void Variations::CheckConsistencyWithBeamSpectra(BEAM::Beam_Spectra_Handler *beamspectra)
{
  // make sure we deal with (anti-)protonic beam, otherwise warn
  bool shouldwarn(false);
  for (int beam(0); beam <= 1; ++beam) {
    if (beamspectra->GetBeam(beam)->Bunch().Kfcode() != kf_p_plus) {
      shouldwarn = true;
      break;
    }
  }
  if (shouldwarn) {
    msg_Error()<<"WARNING in "<<METHOD<<": "<<std::endl
    <<"   The internal reweighting is only tested"<<std::endl
    <<"   for hadronic beams of (anti-)protons."<<std::endl
    <<"   Will continue and hope for the best."<<std::endl;
  }
}

Variations::Variations(Variations_Mode mode)
{
  m_enabled = true;
  if (mode == Variations_Mode::nominal_only)
    return;

  ReadDefaults();
  int lhapdfverbosity(0);
  const bool needslhapdf(NeedsLHAPDF6Interface());
#ifndef USING__LHAPDF
  if (needslhapdf) {
    THROW(fatal_error, "LHAPDF has not been enabled during configuration and "
                       "is required for PDF variations");
  }
#else
  if (needslhapdf) {
    if (!s_loader->LibraryIsLoaded("LHAPDFSherpa")) {
      THROW(fatal_error, "LHAPDF interface is not initialised.");
    }
    lhapdfverbosity = LHAPDF::verbosity();
    LHAPDF::setVerbosity(0);
  }
  InitialiseParametersVector();

  if (!m_parameters_vector.empty() || !m_qcut_parameters_vector.empty()) {
    rpa->gen.AddCitation(1, "The Sherpa-internal reweighting is published in "
                            "\\cite{Bothmann:2016nao}.");
  }

  if (needslhapdf) {
    LHAPDF::setVerbosity(lhapdfverbosity);
  }
#endif
}


Variations::~Variations()
{
  for (auto params : m_parameters_vector)
    delete params;
}

std::string
Variations::GetVariationNameAt(Variations::Parameters_Vector::size_type i,
                               Variations_Type t, Variations_Source s,
                               Variations_Name_Type name_type) const {
  switch (t) {
  case Variations_Type::qcd:
    return m_parameters_vector.at(i)->Name(s, name_type);
  case Variations_Type::qcut:
    return m_qcut_parameters_vector.at(i).Name(s);
  case Variations_Type::custom:
    THROW(fatal_error, "Variations does not manage custom variations.");
  }
  return std::string();
}

size_t Variations::Size(Variations_Type t) const
{
  if (!m_enabled)
    return 0;
  switch (t) {
  case Variations_Type::qcd:
    return m_parameters_vector.size();
  case Variations_Type::qcut:
    return m_qcut_parameters_vector.size();
  case Variations_Type::custom:
    THROW(fatal_error, "Variations does not manage custom variations.");
  }
  return 0;
}

void Variations::ReadDefaults()
{
  Settings& s = Settings::GetMainSettings();
  m_includecentralvaluevariation =
    s["VARIATIONS_INCLUDE_CV"].SetDefault(false).Get<bool>();
  const bool cssreweight = s["SHOWER"]["REWEIGHT"].Get<bool>();
  m_reweightsplittingalphasscales =
    s["REWEIGHT_SPLITTING_ALPHAS_SCALES"].SetDefault(cssreweight).Get<bool>();
  m_reweightsplittingpdfsscales =
    s["REWEIGHT_SPLITTING_PDF_SCALES"].SetDefault(cssreweight).Get<bool>();
  if (!cssreweight && m_reweightsplittingalphasscales) {
    msg_Out() << Strings::OverrideAlphaSSplittingScalesWarning;
    s["REWEIGHT_SPLITTING_ALPHAS_SCALES"].OverrideScalar(false);
    m_reweightsplittingalphasscales = false;
  }
  if (!cssreweight && m_reweightsplittingpdfsscales) {
    msg_Out() << Strings::OverridePDFSplittingScalesWarning;
    s["REWEIGHT_SPLITTING_PDF_SCALES"].OverrideScalar(false);
    m_reweightsplittingpdfsscales = false;
  }
}


void Variations::PrintStatistics(std::ostream &str)
{
  // collect all warnings into a (warning -> number of reports) map and a
  // (warning -> (reporter that reported most, number of reports)) map
  std::map<std::string, unsigned long> warningnums;
  std::map<std::string, std::pair<std::string, unsigned long> > warningreps;
  for (Parameters_Vector::const_iterator paramsit = m_parameters_vector.begin();
       paramsit != m_parameters_vector.end();
       ++paramsit) {
    for (std::map<std::string, unsigned long>::const_iterator
         warningsit((*paramsit)->m_warnings.begin());
         warningsit != (*paramsit)->m_warnings.end(); ++warningsit) {
      if (warningsit->second > 0) {
        warningnums[warningsit->first] += warningsit->second;
        if (warningreps[warningsit->first].second < warningsit->second) {
          warningreps[warningsit->first].first  = (*paramsit)->Name();
          warningreps[warningsit->first].second = warningsit->second;
        }
      }
    }
  }

  if (m_warnings.empty() && warningnums.empty())
    return;

  str << "Reweighting statistics:\n";
  if (!m_warnings.empty()) {
    // output own warnings (warning -> number of reports)
    for (std::map<std::string, unsigned long>::const_iterator
         it(m_warnings.begin()); it != m_warnings.end(); ++it) {
      str << "  " << it->first << ": " << it->second << std::endl;
    }
  }
  if (!warningnums.empty()) {
    // output parameter-specific warnings (warning -> number of reports)
    for (std::map<std::string, unsigned long>::const_iterator it(warningnums.begin());
        it != warningnums.end(); ++it) {
      str << "  " << it->first << ": " << it->second << std::endl;
    }
    if (msg_LevelIsDebugging()) {
      // output parameter-specific warnings
      // (warning -> (reporter that reported largest number of reports))
      str << "  Parameter variation that reported most warnings:" << std::endl;
      for (std::map<std::string, std::pair<std::string, unsigned long> >::const_iterator
          it(warningreps.begin());
          it != warningreps.end(); ++it) {
        str << "    " << it->first << ": " << it->second.first;
        str << " (" << it->second.second << ")" << std::endl;
      }
    }
  }
}


void Variations::InitialiseParametersVector()
{
  // All variations can be specified either using VARIATIONS (which allows for
  // correlated variations) ...

  Settings& s = Settings::GetMainSettings();
  for (auto single_variation_settings : s["VARIATIONS"].GetItems()) {
    AddParameters(single_variation_settings);
  }

  // ... or by using the specialised (and easier to remember) versions below.

  for (auto single_variation_settings : s["PDF_VARIATIONS"].GetItems()) {
    PDFs_And_AlphaS_List pdfsandalphaslist =
        PDFsAndAlphaSList(single_variation_settings.Get<std::string>());
    AddParameterExpandingScaleFactors({"1.0", "1.0", "1.0"},
                                      ScaleFactorExpansions::None,
                                      pdfsandalphaslist);
  }

  for (auto single_variation_settings : s["SCALE_VARIATIONS"].GetItems()) {
    std::vector<std::string> scalestringparams;
    ScaleFactorExpansions::code scalefactorexpansions(
        ScaleFactorExpansions::None);
    if (single_variation_settings.IsScalar()) {
      const auto val = single_variation_settings.Get<std::string>();
      if (val == "None")
        continue;
      ExpandableVariation var {val};
      if (var.expand)
        scalefactorexpansions |= ScaleFactorExpansions::SevenPoint;
      scalestringparams = {var.var, var.var, "1.0"};
    } else {
      ExpandableVariation mufvar {
          single_variation_settings[0].Get<std::string>()};
      if (mufvar.expand)
        scalefactorexpansions |= ScaleFactorExpansions::MuF;
      ExpandableVariation murvar {
          single_variation_settings[1].Get<std::string>()};
      if (murvar.expand)
        scalefactorexpansions |= ScaleFactorExpansions::MuR;
      scalestringparams = {mufvar.var, murvar.var, "1.0"};
    }
    AddParameterExpandingScaleFactors(scalestringparams, scalefactorexpansions,
                                      {});
  }

  for (auto single_variation_settings : s["QCUT_VARIATIONS"].GetItems()) {
    std::vector<std::string> scalestringparams;
    ScaleFactorExpansions::code scalefactorexpansions(
        ScaleFactorExpansions::None);
    ExpandableVariation var {single_variation_settings.Get<std::string>()};
    if (var.expand)
      scalefactorexpansions |= ScaleFactorExpansions::QCUT;
    scalestringparams = {"1.0", "1.0", var.var};
    AddParameterExpandingScaleFactors(scalestringparams, scalefactorexpansions,
                                      {});
  }
}

void Variations::AddParameters(Scoped_Settings& s)
{
  if (s.IsScalar()) {
    return;
  }

  // parse scale factors and QCUT factors, and check whether they are requested
  // to be expanded to x -> [x, 1/x] via an appended asterisk
  std::vector<std::string> scalestringparams;
  ScaleFactorExpansions::code scalefactorexpansions(ScaleFactorExpansions::None);
  if (s["ScaleFactors"].IsMap()) {
    auto mufac = s["ScaleFactors"]["Mu2"].SetDefault("None").Get<std::string>();
    if (mufac != "None") {
      ExpandableVariation var {mufac};
      if (var.expand)
        scalefactorexpansions |= ScaleFactorExpansions::SevenPoint;
      scalestringparams = {var.var, var.var};
    } else {
      // "Mu2" not given, or left explicitly at its default value; try the
      // individual scale factors MuF2 and MuR2 instead
      auto muffac = s["ScaleFactors"]["MuF2"].SetDefault("1.0").Get<std::string>();
      ExpandableVariation mufvar {muffac};
      if (mufvar.expand)
        scalefactorexpansions |= ScaleFactorExpansions::MuF;
      auto murfac = s["ScaleFactors"]["MuR2"].SetDefault("1.0").Get<std::string>();
      ExpandableVariation murvar {murfac};
      if (murvar.expand)
        scalefactorexpansions |= ScaleFactorExpansions::MuR;
      scalestringparams = {mufvar.var, murvar.var};
    }
    auto qcutfac = s["ScaleFactors"]["QCUT"].SetDefault("1.0").Get<std::string>();
    ExpandableVariation var {qcutfac};
    if (var.expand)
      scalefactorexpansions |= ScaleFactorExpansions::QCUT;
    scalestringparams.push_back(var.var);
  } else {
    // just a scalar is given for s["ScaleFactors"], interprete it though it
    // was given as s["ScaleFactors"]["Mu2"]
    auto mufac = s["ScaleFactors"].SetDefault("1.0").Get<std::string>();
    ExpandableVariation var {mufac};
    if (var.expand)
      scalefactorexpansions |= ScaleFactorExpansions::SevenPoint;
    scalestringparams = {var.var, var.var, "1.0"};
  }

  // parse PDF set, and check whether it is requested to be expanded to all its
  // members via an appended asterisk
  PDFs_And_AlphaS_List pdfsandalphaslist =
      PDFsAndAlphaSList(s["PDF"].SetDefault("None").Get<std::string>());

  // parse AlphaS(MZ) variations and append them
  auto alphasmz = s["AlphaS(MZ)"].SetDefault(-1.0).Get<double>();
  if (alphasmz != -1.0) {
    pdfsandalphaslist.items.push_back(PDFs_And_AlphaS(alphasmz));
  }

  AddParameterExpandingScaleFactors(scalestringparams,
                                    scalefactorexpansions,
                                    pdfsandalphaslist);
}

void Variations::AddParameterExpandingScaleFactors(
    std::vector<std::string> scalestringparams, ScaleFactorExpansions::code expansions,
    PDFs_And_AlphaS_List pdfsandalphaslist)
{
  bool scalevariationrequested(!scalestringparams.empty());

  const double muF2fac {ToType<double>(scalestringparams[0])};
  const double muR2fac {ToType<double>(scalestringparams[1])};
  const double qcutfac {ToType<double>(scalestringparams[2])};

  if (qcutfac != 1.0) {
    if (muF2fac != 1.0 || muR2fac != 1.0 || !pdfsandalphaslist.items.empty()) {
      THROW(not_implemented,
            "Simultaneous variations of QCUT and QCD"
            " parameters (muF2, muR2, PDFs/AlphaS(mZ)) are not supported.")
    }
    m_qcut_parameters_vector.push_back({qcutfac});
    if (expansions & ScaleFactorExpansions::QCUT) {
      if (m_includecentralvaluevariation)
        m_qcut_parameters_vector.push_back({1.0});
      m_qcut_parameters_vector.push_back({1.0 / qcutfac});
    } else {
      m_qcut_parameters_vector.push_back({qcutfac});
    }
    return;
  }

  // translate muF2fac and muR2fac into scale-factor pairs, expanding if
  // necessary

  std::vector<ScaleFactor_Pair> scalefactorpairs;
  bool did_expand {true};
  if (expansions & ScaleFactorExpansions::SevenPoint) {
    scalefactorpairs.push_back(ScaleFactor_Pair(muR2fac, 1.0));
    scalefactorpairs.push_back(ScaleFactor_Pair(1.0, muF2fac));
    scalefactorpairs.push_back(ScaleFactor_Pair(muR2fac, muF2fac));
    scalefactorpairs.push_back(ScaleFactor_Pair(1.0, 1.0));
    scalefactorpairs.push_back(ScaleFactor_Pair(1.0 / muR2fac, 1.0));
    scalefactorpairs.push_back(ScaleFactor_Pair(1.0, 1.0 / muF2fac));
    scalefactorpairs.push_back(ScaleFactor_Pair(1.0 / muR2fac, 1.0 / muF2fac));
  } else if (expansions & ScaleFactorExpansions::MuR &&
             !(expansions & ScaleFactorExpansions::MuF)) {
    scalefactorpairs.push_back(ScaleFactor_Pair(muR2fac, muF2fac));
    scalefactorpairs.push_back(ScaleFactor_Pair(1.0, muF2fac));
    scalefactorpairs.push_back(ScaleFactor_Pair(1.0 / muR2fac, muF2fac));
  } else if (expansions & ScaleFactorExpansions::MuF &&
             !(expansions & ScaleFactorExpansions::MuR)) {
    scalefactorpairs.push_back(ScaleFactor_Pair(muR2fac, muF2fac));
    scalefactorpairs.push_back(ScaleFactor_Pair(muR2fac, 1.0));
    scalefactorpairs.push_back(ScaleFactor_Pair(muR2fac, 1.0 / muF2fac));
  } else if (expansions & ScaleFactorExpansions::MuF &&
             expansions & ScaleFactorExpansions::MuR) {
    scalefactorpairs.push_back(ScaleFactor_Pair(muR2fac, 1.0));
    scalefactorpairs.push_back(ScaleFactor_Pair(1.0, muF2fac));
    scalefactorpairs.push_back(ScaleFactor_Pair(muR2fac, muF2fac));
    scalefactorpairs.push_back(ScaleFactor_Pair(1.0 / muR2fac, 1.0));
    scalefactorpairs.push_back(ScaleFactor_Pair(1.0, 1.0));
    scalefactorpairs.push_back(ScaleFactor_Pair(1.0, 1.0 / muF2fac));
    scalefactorpairs.push_back(ScaleFactor_Pair(1.0 / muR2fac, 1.0 / muF2fac));
    scalefactorpairs.push_back(ScaleFactor_Pair(1.0 / muR2fac, muF2fac));
    scalefactorpairs.push_back(ScaleFactor_Pair(muR2fac, 1.0 / muF2fac));
  } else {
    scalefactorpairs.push_back(ScaleFactor_Pair(muR2fac, muF2fac));
    did_expand = false;
  }

  // if there is no explicit PDF requested, we use the nominal one
  if (pdfsandalphaslist.items.empty()) {
    pdfsandalphaslist.items.push_back(PDFs_And_AlphaS());
    pdfsandalphaslist.did_expand = true;
  }

  // now add all possible combinations to our list of parameters
  for (std::vector<PDFs_And_AlphaS>::const_iterator pdfasit(pdfsandalphaslist.items.begin());
        pdfasit != pdfsandalphaslist.items.end(); pdfasit++) {
    bool assignedownershipofpdfsandalphas {false};
    for (const auto sfpair : scalefactorpairs) {
      if (!m_includecentralvaluevariation && did_expand &&
          pdfsandalphaslist.did_expand &&
          IsQCDVariationTrivial(sfpair.first, sfpair.second, pdfasit->m_pdfs[0],
                                pdfasit->m_pdfs[1], pdfasit->p_alphas)) {
        continue;
      }
      int shoulddeletepdfmask {0};
      bool shoulddeletealphas {false};
      if (!assignedownershipofpdfsandalphas) {
        shoulddeletepdfmask = pdfasit->m_shoulddeletepdfmask;
        shoulddeletealphas = pdfasit->m_shoulddeletealphas;
        assignedownershipofpdfsandalphas = true;
      }
      AddParameters(
          sfpair.first,
          sfpair.second,
          pdfasit,
          shoulddeletepdfmask,
          shoulddeletealphas);
    }
  }
}


void Variations::AddParameters(double muR2fac, double muF2fac,
                               std::vector<PDFs_And_AlphaS>::const_iterator pdfsandalphas,
                               int deletepdfmask,
                               bool deletealphas)
{
  const bool showermuR2enabled {m_reweightsplittingalphasscales};
  const bool showermuF2enabled {m_reweightsplittingpdfsscales};
  QCD_Variation_Params *params =
    new QCD_Variation_Params(
	muR2fac, muF2fac, showermuR2enabled, showermuF2enabled,
        pdfsandalphas->m_pdfs[0], pdfsandalphas->m_pdfs[1],
        pdfsandalphas->p_alphas,
        deletepdfmask, deletealphas);
  m_parameters_vector.push_back(params);
}

Variations::PDFs_And_AlphaS_List
Variations::PDFsAndAlphaSList(const std::string &pdf) const
{
  // parse PDF set, and check whether it is requested to be expanded to all its
  // members via an appended asterisk
  PDFs_And_AlphaS_List pdfsandalphaslist;
  if (pdf != "None") {
    Settings& s = Settings::GetMainSettings();
    if (s["ALPHAS"]["USE_PDF"].Get<int>() == 0) {
      THROW(fatal_error,
            "`ALPHAS: {USE_PDF: 0}` is incompatible with doing PDF/AlphaS "
            "variations.");
    }
    int alphasbeam {
      s["PDF_VARIATION_ALPHAS_BEAM"].SetDefault(0).Get<int>() - 1};
    int beammask {0};
    for (int i : s["PDF_VARIATION_BEAMS"].SetDefault({1, 2}).GetVector<int>()) {
      if (rpa->gen.Bunch(i - 1).IsHadron() || rpa->gen.Bunch(i - 1).IsPhoton()) {
        beammask |= 1 << (i - 1);
        if (alphasbeam < 0)
          alphasbeam = i - 1;
      }
    }

    // translate PDF string parameter into actual AlphaS and PDF objects
    ExpandableVariation var{pdf};
    pdfsandalphaslist =
        PDFsAndAlphaSList(var.var, var.expand, beammask, alphasbeam);
  }
  return pdfsandalphaslist;
}

Variations::PDFs_And_AlphaS_List Variations::PDFsAndAlphaSList(
    std::string pdfstringparam,
    bool expandpdf,
    int beammask,
    int alphasbeam) const
{
  // parse PDF member(s)
  PDFs_And_AlphaS_List pdfsandalphaslist;
  size_t firstmember(0);
  size_t lastmember(0);
  if (expandpdf) {
    pdfsandalphaslist.did_expand = true;
    // determine the number of set members to load
    bool lhapdflookupsuccessful(false);
    #ifdef USING__LHAPDF
      if (s_loader->LibraryIsLoaded("LHAPDFSherpa")) {
        const std::vector<std::string>& availablepdfsets(LHAPDF::availablePDFSets());
        if (std::find(availablepdfsets.begin(), availablepdfsets.end(), pdfstringparam)
            != availablepdfsets.end()) {
          LHAPDF::PDFSet set(pdfstringparam);
          lastmember = set.size() - 1; // this assumes members are labelled 0...(set.size()-1)
          lhapdflookupsuccessful = true;
        }
      }
      if (!lhapdflookupsuccessful) {
        // the LHAPDF interface is not available or does not know about this set
        // provide a fallback at least for the default PDF set
        if (pdfstringparam == "PDF4LHC21_40_pdfas") {
          lastmember = 42;
        } else {
          THROW(not_implemented,
                "Full PDF set reweightings only work with LHAPDF6 sets or the internal default set."
                + std::string(" Otherwise specify explicitly."));
        }
    }
    #endif
  } else {
    // single PDF member: "Set/i" or just "Set"
    int member(0);
    if (pdfstringparam.find("/") != std::string::npos) {
      member = ToType<int>(std::string(pdfstringparam, pdfstringparam.find("/") + 1));
      pdfstringparam = pdfstringparam.substr(0, pdfstringparam.find("/"));
    }
    firstmember = member;
    lastmember = member;
  }
  for (size_t j(firstmember); j <= lastmember; ++j) {
    pdfsandalphaslist.items.push_back(
        PDFs_And_AlphaS(pdfstringparam, j, beammask, alphasbeam));
  }
  return pdfsandalphaslist;
}


Variations::PDFs_And_AlphaS::PDFs_And_AlphaS():
    p_alphas(MODEL::as)
{
  // Workaround für C++03 (vector constructor is confused when fed
  // with NULL as the initial value for a pointer)
  // cf. https://gcc.gnu.org/ml/gcc-help/2013-02/msg00026.html
  PDF::PDF_Base *nullPtr = NULL;
  m_pdfs = std::vector<PDF::PDF_Base *>(2, nullPtr);
  m_pdfs[0] = rpa->gen.PDF(0);
  m_pdfs[1] = rpa->gen.PDF(1);
}


Variations::PDFs_And_AlphaS::PDFs_And_AlphaS(double alphasmz)
{
  // Workaround für C++03 (vector constructor is confused when fed
  // with NULL as the initial value for a pointer)
  // cf. https://gcc.gnu.org/ml/gcc-help/2013-02/msg00026.html
  PDF::PDF_Base *nullPtr = NULL;
  m_pdfs = std::vector<PDF::PDF_Base *>(2, nullPtr);
  m_pdfs[0] = rpa->gen.PDF(0);
  m_pdfs[1] = rpa->gen.PDF(1);

  Settings& s = Settings::GetMainSettings();
  int order_alphaS = s["ORDER_ALPHAS"].Get<int>();
  int th_alphaS = s["THRESHOLD_ALPHAS"].Get<int>();
  double MZ2 = sqr(Flavour(kf_Z).Mass());
  p_alphas = new MODEL::Running_AlphaS(alphasmz, MZ2, order_alphaS, th_alphaS);
  m_shoulddeletealphas = true;
}


Variations::PDFs_And_AlphaS::PDFs_And_AlphaS(
    std::string pdfname, int pdfmember, int beammask, int alphasbeam)
{
  // obtain PDFs
  for (int i(0); i < 2; ++i) {
    if ((1 << i) & beammask) {
      PDF::PDF_Arguments args{ rpa->gen.Bunch(i), i, pdfname, pdfmember };
      PDF::PDF_Base *pdf = PDF::PDF_Base::PDF_Getter_Function::GetObject(pdfname, args);
      if (pdf == NULL) THROW(fatal_error, "PDF set " + pdfname + " not available.");
      pdf->SetBounds();
      m_pdfs.push_back(pdf);
    } else if (rpa->gen.Bunch(i).IsHadron() || rpa->gen.Bunch(i).IsPhoton()) {
      m_pdfs.push_back(rpa->gen.PDF(i));
    } else {
      m_pdfs.push_back(NULL);
    }
  }
  m_shoulddeletepdfmask = beammask;

  // obtain AlphaS based on a loaded PDF or a new one (if none is found)
  PDF::PDF_Base *aspdf {m_pdfs[alphasbeam]};
  if (aspdf == NULL) {
    p_alphas = new MODEL::Running_AlphaS(pdfname, pdfmember);
    m_shoulddeletealphas = true;
  } else {
    p_alphas = new MODEL::Running_AlphaS(aspdf);
    m_shoulddeletealphas = (1 << alphasbeam) & beammask;
  }
  if (p_alphas == NULL) {
    THROW(fatal_error, "AlphaS for " + pdfname + " could not be initialised.");
  }
}


#if ENABLE_REWEIGHTING_FACTORS_HISTOGRAMS
ReweightingFactorHistogram::ReweightingFactorHistogram()
{
  // initialise histogram (many steps around 1.0, less far from 1.0)
  double binwidth(0.1);
  double distancefromone(0.0);
  while (distancefromone < 100.0) {
     Bin right(1.0 + distancefromone, EntryNumbers());
     Bin left(1.0 - distancefromone - binwidth, EntryNumbers());
     m_bins.insert(m_bins.begin(), left);
     m_bins.push_back(right);
     distancefromone += binwidth;
     binwidth *= 1.3;
  }
}


void ReweightingFactorHistogram::Fill(std::string name, double value)
{
  // find entry key index and ensure the key and the entres are present
  std::vector<EntryKey>::iterator findit =
  std::find(m_keys.begin(), m_keys.end(), name);
  size_t index(std::distance(m_keys.begin(), findit));
  if (findit == m_keys.end()) {
    m_keys.push_back(name);
    for (std::vector<Bin>::iterator it = m_bins.begin();
         it != m_bins.end(); ++it) {
      it->second.push_back(0);
    }
  }

  // increment entry number for appropriate bin
  for (std::vector<Bin>::iterator it = m_bins.begin();
       it != m_bins.end(); ++it) {
    if ((it+1) == m_bins.end() || value < (it+1)->first) {
      it->second[index] += + 1;
      break;
    }
  }
}


void ReweightingFactorHistogram::Write(std::string filenameaffix)
{
  Settings& s = Settings::GetMainSettings();
  // open file
  std::string outpath = s.Get<std::string>("ANALYSIS_OUTPUT");
  if (outpath.length() > 0) {
    size_t slashpos = outpath.rfind("/");
    if (slashpos != std::string::npos) {
      MakeDir(outpath.substr(0, slashpos));
    }
    if (slashpos != outpath.length() - 1) {
      outpath += ".";
    }
  }
  std::string filename = outpath + filenameaffix + ".RewFacHisto.dat";
  std::ofstream file;
  file.open(filename.c_str());

  // write head
  file << "# leftbinedge";
  for (std::vector<EntryKey>::const_iterator it(m_keys.begin());
       it != m_keys.end(); ++it) {
    file << " " << *it;
  }
  file << std::endl;

  // write bins
  for (std::vector<Bin>::const_iterator binit(m_bins.begin());
       binit != m_bins.end(); ++binit) {
    file << binit->first;
    EntryNumbers nums(binit->second);
    for (EntryNumbers::const_iterator numit(nums.begin());
         numit != nums.end(); ++numit) {
      file << " " << *numit;
    }
    file << std::endl;
  }
  file.close();
}
#endif


QCD_Variation_Params::~QCD_Variation_Params()
{
  if (m_deletepdfmask & 1) {
    if (p_pdf1) { delete p_pdf1; }
  }
  if (m_deletepdfmask & 2) {
    if (p_pdf2) { delete p_pdf2; }
  }
  if (m_deletealphas) {
    if (p_alphas) { delete p_alphas; }
  }
#if ENABLE_REWEIGHTING_FACTORS_HISTOGRAMS
  m_rewfachisto.Write(Name());
#endif
}

std::string QCD_Variation_Params::Name(Variations_Source source,
                                       Variations_Name_Type name_type) const {
  std::string divider;
  std::string name;
  std::string prefix;

  if (name_type == Variations_Name_Type::weight_name_convention) {
    divider = "__";
    if (source == Variations_Source::main)
      prefix = "ME.";
  }
  else {
    divider = " ";
    if (source == Variations_Source::main)
      name = "ME-only:" + divider;
    else
      name = "ME & PS:" + divider;
    if (IsTrivial())
      return name + "<no variation>";
  }

  // scale factors tags; note that we assume here that the shower scale factors
  // are either equal to the matrix-element scale factors, or 1.0
  if (name_type != Variations_Name_Type::human_readable || m_muR2fac != 1.0) {
    if (source == Variations_Source::main || !m_showermuR2enabled) {
      name += GenerateVariationNamePart(prefix + "MUR=", sqrt(m_muR2fac));
    } else {
      name += GenerateVariationNamePart("MUR=", sqrt(m_muR2fac));
    }
    name += divider;
  }
  if (name_type != Variations_Name_Type::human_readable || m_muF2fac != 1.0) {
    if (source == Variations_Source::main || !m_showermuF2enabled) {
      name += GenerateVariationNamePart(prefix + "MUF=", sqrt(m_muF2fac));
    } else {
      name += GenerateVariationNamePart("MUF=", sqrt(m_muF2fac));
    }
    name += divider;
  }

  // PDF tags
  if (p_pdf1 == NULL || p_pdf2 == NULL || p_pdf1->LHEFNumber() == p_pdf2->LHEFNumber()) {

    // there is only one relevant PDF ID (or none)
    PDF::PDF_Base* pdf {nullptr};
    if (p_pdf1 != NULL) {
      pdf = p_pdf1;
    } else if (p_pdf2 != NULL) {
      pdf = p_pdf2;
    } else if (p_alphas->GetAs()->PDF() != NULL) {
      pdf = p_alphas->GetAs()->PDF();
    }
    int pdfid(-1);
    if (pdf)
      pdfid = pdf->LHEFNumber();
    if (name_type != Variations_Name_Type::human_readable || pdf != rpa->gen.PDF(0)) {
      name += prefix + GenerateVariationNamePart("LHAPDF=", pdfid);
    }
  } else {
    // there are two relevant PDF IDs, quote both
    name += prefix +
            GenerateVariationNamePart("LHAPDF.BEAM1=", p_pdf1->LHEFNumber()) +
            divider + prefix +
            GenerateVariationNamePart("LHAPDF.BEAM2=", p_pdf2->LHEFNumber());
  }

  // append optional non-trival AlphaS(MZ) variation (which is not related to a
  // change in the PDF set)
  if (p_alphas != MODEL::as && p_alphas->GetAs()->PDF() != p_pdf1) {
    name += divider + prefix + GenerateVariationNamePart("ASMZ=", p_alphas->AsMZ());
  }

  return name;
}

bool QCD_Variation_Params::IsTrivial() const
{
  return IsQCDVariationTrivial(m_muR2fac, m_muF2fac, p_pdf1, p_pdf2, p_alphas);
}

std::string Qcut_Variation_Params::Name(Variations_Source source) const
{
  std::string level_prefix;
  if (source == Variations_Source::main)
    level_prefix = "ME.";
  return level_prefix + GenerateVariationNamePart("QCUT", m_scale_factor);
}

namespace ATOOLS {

  std::ostream& operator<<(std::ostream& o, const Variations_Type& t)
  {
    switch (t) {
      case Variations_Type::qcd:    return o << "QCD";
      case Variations_Type::qcut:   return o << "Qcut";
      case Variations_Type::custom: return o << "Custom";
    }
    return o;
  }

  std::ostream& operator<<(std::ostream& s, const Variations& v)
  {
    const Variations::Parameters_Vector * const paramsvec(v.GetParametersVector());
    s << "Named variations:" << std::endl;
    if (paramsvec->empty()) {
      return s << " None\n";
    }
    s << '\n';
    for (Variations::Parameters_Vector::const_iterator it(paramsvec->begin());
         it != paramsvec->end(); ++it) {
      s << (*it)->Name() << " (" << (*it)->m_deletepdfmask << ","
        << (*it)->m_deletealphas << ")" << '\n';
    }
    return s;
  }

  std::ostream& operator<<(std::ostream& s, const Variations_Source &c)
   {
    switch (c) {
      case Variations_Source::all:
        return s << "All";
      case Variations_Source::main:
        return s << "Main";
      case Variations_Source::sudakov:
        return s << "Sudakov";
    }
    return s;
  }

  bool IsQCDVariationTrivial(
      double muR2fac, double muF2fac,
      PDF::PDF_Base * const pdf1,
      PDF::PDF_Base * const pdf2,
      MODEL::Running_AlphaS * const alphas)
  {
    if (muR2fac != 1.0
        || muF2fac != 1.0
        || pdf1 != rpa->gen.PDF(0)
        || pdf2 != rpa->gen.PDF(1)
        || alphas != MODEL::as
       )
      return false;
    return true;
  }
}
