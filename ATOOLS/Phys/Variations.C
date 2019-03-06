#include "ATOOLS/Phys/Variations.H"

#include <iterator>
#include <numeric>
#include <algorithm>

#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Blob.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/PDF_Base.H"
#if defined USING__LHAPDF && defined USING__LHAPDF6
#include "LHAPDF/LHAPDF.h"
#endif

namespace ATOOLS {
  struct ScaleFactor_Pair: public std::pair<double, double> {
    double m_qcutfac;
    ScaleFactor_Pair(const double &a,const double &b,const double &c):
      std::pair<double, double>(a,b), m_qcutfac(c) {}
  };
}

using namespace ATOOLS;
using namespace SHERPA;

bool Variations::NeedsLHAPDF6Interface()
{
  // go through all parameters and return true if any non-double points to a
  // request for a PDF set/member
  std::vector<std::string> args = VariationArguments();
  if (args.empty()) {
    return false;
  }
  for (std::vector<std::string>::const_iterator paramsit(args.begin());
      paramsit != args.end(); ++paramsit) {
    std::vector<std::string> params(VariationArgumentParameters(*paramsit));
    for (std::vector<std::string>::const_iterator it(params.begin());
        it != params.end(); ++it) {
      if ((*it).empty()) continue;
      std::string stringparam(*it);
      // remove openening/closing brackets
      if (*stringparam.begin() == '[') {
        stringparam.erase(stringparam.begin(), stringparam.begin() + 1);
      }
      if (*(stringparam.end() - 1) == ']') {
        stringparam.erase(stringparam.end() - 1, stringparam.end());
      }
      // check if we can parse as a double
      double testdouble(1.0);
      MyStrStream ss(stringparam);
      ss >> testdouble;
      if (ss.fail()) {
        return true;
      }
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

Variations::Variations()
{
  ReadDefaults();
#if defined USING__LHAPDF && defined USING__LHAPDF6
  int lhapdfverbosity(0);
  const bool needslhapdf(NeedsLHAPDF6Interface());
  if (needslhapdf) {
    if (!s_loader->LibraryIsLoaded("LHAPDFSherpa")) {
      THROW(fatal_error, "LHAPDF interface is not initialised.");
    }
    lhapdfverbosity = LHAPDF::verbosity();
    LHAPDF::setVerbosity(0);
  }
#endif

  InitialiseParametersVector();

  if (!m_parameters_vector.empty()) {
    rpa->gen.AddCitation(1, "The Sherpa-internal reweighting is published in \\cite{Bothmann:2016nao}.");
  }

#if defined USING__LHAPDF && defined USING__LHAPDF6
  if (needslhapdf) {
    LHAPDF::setVerbosity(lhapdfverbosity);
  }
#endif
}


Variations::~Variations()
{
  for (Parameters_Vector::const_iterator it = m_parameters_vector.begin();
       it != m_parameters_vector.end();
       ++it) {
    delete *it;
  }
}


void Variations::ReadDefaults()
{
  Settings& s = Settings::GetMainSettings();
  m_includecentralvaluevariation =
    s["VARIATIONS_INCLUDE_CV"].SetDefault(true).Get<bool>();
  m_reweightsplittingalphasscales =
    s["REWEIGHT_SPLITTING_ALPHAS_SCALES"].SetDefault(false).Get<bool>();
  m_reweightsplittingpdfsscales =
    s["REWEIGHT_SPLITTING_PDF_SCALES"].SetDefault(false).Get<bool>();
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
          warningreps[warningsit->first].first  = (*paramsit)->m_name;
          warningreps[warningsit->first].second = warningsit->second;
        }
      }
    }
  }

  if (!m_warnings.empty() || !warningnums.empty()) {
    str << METHOD << "(): Reweighting statistics {" << std::endl;
  }
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
  if (!m_warnings.empty() || !warningnums.empty()) {
    str << "}" << std::endl;
  }
}


void Variations::InitialiseParametersVector()
{
  std::vector<std::string> args = VariationArguments();
  for (std::vector<std::string>::const_iterator it(args.begin());
      it != args.end(); ++it) {
    std::vector<std::string> params(VariationArgumentParameters(*it));
    AddParameters(params);
  }
}


std::vector<std::string> Variations::VariationArguments()
{
  Settings& s = Settings::GetMainSettings();
  auto args = s["VARIATIONS"].GetVector<std::string>();
  args.erase(std::remove(args.begin(), args.end(), "None"), args.end());
  return args;
}


std::vector<std::string> Variations::VariationArgumentParameters(std::string arg)
{
  return ToVector<std::string>(arg, ',');
}


void Variations::AddParameters(std::vector<std::string> stringparams)
{
  // up to 2 scale factor arguments and up to 1 PDF argument
  std::vector<std::string> scalestringparams;
  std::vector<std::string> pdfstringparams;

  ScaleFactorExpansions::code scalefactorexpansions(ScaleFactorExpansions::None);
  bool expandpdf(false);

  for (std::vector<std::string>::const_iterator it(stringparams.begin());
      it != stringparams.end(); it++) {

    if ((*it).empty()) continue;

    std::string stringparam(*it);

    // check and remove openening/closing brackets
    bool bracketopened(false);
    bool bracketclosed(false);
    if (*stringparam.begin() == '[') {
      bracketopened = true;
      stringparam.erase(stringparam.begin(), stringparam.begin() + 1);
    }
    if (*(stringparam.end() - 1) == ']') {
      bracketclosed = true;
      stringparam.erase(stringparam.end() - 1, stringparam.end());
    }

    // check if it's a scale factor or a PDF set/member and record the
    // consequence in terms of bracket expansion
    double sf(1.0);
    MyStrStream ss(stringparam);
    ss >> sf;
    if (ss.fail()) {
      // PDF set/member
      if (pdfstringparams.size() > 0) {
        THROW(not_implemented,
              "Only one PDF argument per variation is supported.");
      }
      if (bracketopened ^ bracketclosed) {
        THROW(fatal_error, "Unmatched bracket in PDF argument.");
      }
      expandpdf = bracketopened;
      pdfstringparams.push_back(stringparam);
    } else {
      // scale factor
      if (pdfstringparams.size() > 1) {
        THROW(not_implemented,
              "Only two scale factor arguments per variation are supported.");
      }
      if (bracketopened && bracketclosed) {
        if (scalestringparams.empty()) {
          scalefactorexpansions |= ScaleFactorExpansions::First;
        } else {
          scalefactorexpansions |= ScaleFactorExpansions::Second;
        }
      } else if (!bracketopened && bracketclosed) {
        scalefactorexpansions |= ScaleFactorExpansions::SevenPoint;
      }
      if ((bracketopened || bracketclosed) && sf == 1.0) {
        THROW(inconsistent_option,
              "When brackets are used to expand a scale factor, it must not be 1.0.");
      }
      scalestringparams.push_back(stringparam);
    }
  }

  // translate PDF string parameter into actual AlphaS and PDF objects
  std::vector<Variations::PDFs_And_AlphaS> pdfsandalphasvector;
  if (!pdfstringparams.empty()) {
    Settings& s = Settings::GetMainSettings();
    if (s["OVERRIDE_PDF_INFO"].Get<bool>()) {
      THROW(fatal_error,
            "`OVERRIDE_PDF_INFO: true` is incompatible with doing PDF/AlphaS variations.");
    }
    pdfsandalphasvector = PDFsAndAlphaSVector(pdfstringparams[0], expandpdf);
  }

  AddParameterExpandingScaleFactors(scalestringparams,
                                    scalefactorexpansions,
                                    pdfsandalphasvector);
}

void Variations::AddParameterExpandingScaleFactors(
    std::vector<std::string> scalestringparams, ScaleFactorExpansions::code expansions,
    std::vector<Variations::PDFs_And_AlphaS> pdfsandalphasvector)
{
  bool scalevariationrequested(!scalestringparams.empty());

  // read input scale factors
  double muR2fac(1.0);
  double muF2fac(1.0);
  double Qcutfac(1.0);
  if (scalevariationrequested) {
    if (scalestringparams.size() == 1) {
      muR2fac = muF2fac = ToType<double>(scalestringparams[0]);
    } else {
      muR2fac = ToType<double>(scalestringparams[0]);
      muF2fac = ToType<double>(scalestringparams[1]);
    }
    if (scalestringparams.size()>2)
      Qcutfac = ToType<double>(scalestringparams[2]);
  }

  // expand scale factors
  std::vector<ScaleFactor_Pair > scalefactorpairs;
  ScaleFactor_Pair defaultscalefactorpair(1.0, 1.0,1.0);

  if (!scalevariationrequested || expansions == ScaleFactorExpansions::None) {
    defaultscalefactorpair = ScaleFactor_Pair(muR2fac, muF2fac,Qcutfac);
    scalefactorpairs.push_back(defaultscalefactorpair);
  } else if (expansions & ScaleFactorExpansions::SevenPoint) {
    scalefactorpairs.push_back(ScaleFactor_Pair(muR2fac, 1.0,1.0));
    scalefactorpairs.push_back(ScaleFactor_Pair(1.0, muF2fac,1.0));
    scalefactorpairs.push_back(ScaleFactor_Pair(muR2fac, muF2fac,1.0));
    scalefactorpairs.push_back(ScaleFactor_Pair(1.0 / muR2fac, 1.0,1.0));
    scalefactorpairs.push_back(ScaleFactor_Pair(1.0, 1.0 / muF2fac,1.0));
    scalefactorpairs.push_back(ScaleFactor_Pair(1.0 / muR2fac, 1.0 / muF2fac,1.0));
  } else {
    if (expansions & ScaleFactorExpansions::First) {
      const double defaultmuF2fac = (expansions & ScaleFactorExpansions::Second) ? 1.0 : muF2fac;
      defaultscalefactorpair.second = defaultmuF2fac;
      scalefactorpairs.push_back(ScaleFactor_Pair(muR2fac, defaultmuF2fac,1.0));
      scalefactorpairs.push_back(ScaleFactor_Pair(1.0 / muR2fac, defaultmuF2fac,1.0));
    }
    if (expansions & ScaleFactorExpansions::Second) {
      const double defaultmuR2fac = (expansions & ScaleFactorExpansions::First) ? 1.0 : muR2fac;
      defaultscalefactorpair.first = defaultmuR2fac;
      scalefactorpairs.push_back(ScaleFactor_Pair(defaultmuR2fac, muF2fac,1.0));
      scalefactorpairs.push_back(ScaleFactor_Pair(defaultmuR2fac, 1.0 / muF2fac,1.0));
    }
  }

  // if there is no explicit PDF requested, we use the nominal one
  bool ownedpdfsandalphas(true);
  if (pdfsandalphasvector.empty()) {
    ownedpdfsandalphas = false;
    pdfsandalphasvector.push_back(PDFs_And_AlphaS());
  }

  for (std::vector<PDFs_And_AlphaS>::const_iterator pdfasit(pdfsandalphasvector.begin());
        pdfasit != pdfsandalphasvector.end(); pdfasit++) {

    bool assignedownershipofpdfsandalphas(!ownedpdfsandalphas);

    if (pdfasit != pdfsandalphasvector.begin()
        || (m_includecentralvaluevariation && expansions != ScaleFactorExpansions::None)) {
      AddParameters(defaultscalefactorpair.first, defaultscalefactorpair.second,
		    defaultscalefactorpair.m_qcutfac, pdfasit,
                    !assignedownershipofpdfsandalphas);
      assignedownershipofpdfsandalphas = true;
    }

    if (pdfasit == pdfsandalphasvector.begin()) {
      for (std::vector<ScaleFactor_Pair >::const_iterator scalefactorpairit(scalefactorpairs.begin());
          scalefactorpairit != scalefactorpairs.end(); scalefactorpairit++) {
        AddParameters(scalefactorpairit->first, scalefactorpairit->second,
		      scalefactorpairit->m_qcutfac, pdfasit,
                      !assignedownershipofpdfsandalphas);
        assignedownershipofpdfsandalphas = true;
      }
    }
  }
}


void Variations::AddParameters(double muR2fac, double muF2fac, double Qcutfac,
                               std::vector<PDFs_And_AlphaS>::const_iterator pdfsandalphas,
                               bool deletepdfsandalphas)
{
  const double showermuR2fac = (m_reweightsplittingalphasscales) ? muR2fac : 1.0;
  const double showermuF2fac = (m_reweightsplittingpdfsscales) ? muF2fac : 1.0;
  Variation_Parameters *params =
    new Variation_Parameters(
	muR2fac, muF2fac, showermuR2fac, showermuF2fac, Qcutfac,
        pdfsandalphas->m_pdfs[0], pdfsandalphas->m_pdfs[1],
        pdfsandalphas->p_alphas,
        deletepdfsandalphas);
  m_parameters_vector.push_back(params);
}


std::vector<Variations::PDFs_And_AlphaS> Variations::PDFsAndAlphaSVector(
    std::string pdfstringparam,
    bool expandpdf)
{
  // parse PDF member(s)
  std::vector<PDFs_And_AlphaS> pdfsandalphasvector;
  size_t firstmember(0);
  size_t lastmember(0);
  if (expandpdf) {
    // determine the number of set members to load
    bool lhapdflookupsuccessful(false);
    if (s_loader->LibraryIsLoaded("LHAPDFSherpa")) {
#if defined USING__LHAPDF && defined USING__LHAPDF6
      const std::vector<std::string>& availablepdfsets(LHAPDF::availablePDFSets());
      if (std::find(availablepdfsets.begin(), availablepdfsets.end(), pdfstringparam)
          != availablepdfsets.end()) {
        LHAPDF::PDFSet set(pdfstringparam);
        lastmember = set.size() - 1; // this assumes members are labelled 0...(set.size()-1)
        lhapdflookupsuccessful = true;
      }
#endif
    }
    if (!lhapdflookupsuccessful) {
      // the LHAPDF interface is not available or does not know about this set
      // provide a fallback at least for the default PDF set
      if (pdfstringparam == "NNPDF30NNLO") { 
        lastmember = 100;
      } else {
        THROW(not_implemented,
              "Full PDF set reweightings only work with LHAPDF6 sets or the internal default set."
              + std::string(" Otherwise specify explicitly."));
      }
    }
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
    pdfsandalphasvector.push_back(PDFs_And_AlphaS(pdfstringparam, j));
  }
  return pdfsandalphasvector;
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


Variations::PDFs_And_AlphaS::PDFs_And_AlphaS(
    std::string pdfname, int pdfmember)
{
  // obtain PDFs
  PDF::PDF_Base *aspdf(NULL);
  for (int i(0); i < 2; ++i) {
    if (rpa->gen.Bunch(i).IsHadron()) {
      PDF::PDF_Arguments args{ rpa->gen.Bunch(i), i, pdfname, pdfmember };
      PDF::PDF_Base *pdf = PDF::PDF_Base::PDF_Getter_Function::GetObject(pdfname, args);
      if (pdf == NULL) THROW(fatal_error, "PDF set " + pdfname + " not available.");
      pdf->SetBounds();
      m_pdfs.push_back(pdf);
      if (aspdf == NULL) {
        aspdf = pdf;
      }
    } else {
      m_pdfs.push_back(NULL);
    }
  }

  // obtain AlphaS based on a loaded PDF or a new one (if none is found)
  if (aspdf == NULL) {
    p_alphas = new MODEL::Running_AlphaS(pdfname, pdfmember);
  } else {
    p_alphas = new MODEL::Running_AlphaS(aspdf);
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


Variation_Parameters::~Variation_Parameters()
{
  if (m_deletepdfsandalphas) {
    if (p_pdf1) { delete p_pdf1; }
    if (p_pdf2) { delete p_pdf2; }
    if (p_alphas) { delete p_alphas; }
  }
#if ENABLE_REWEIGHTING_FACTORS_HISTOGRAMS
  m_rewfachisto.Write(m_name);
#endif
}


std::string Variation_Parameters::GenerateName() const
{
  const std::string divider("_");
  std::string name;
  if (p_pdf1 == NULL || p_pdf2 == NULL || p_pdf1->LHEFNumber() == p_pdf2->LHEFNumber()) {
    // there is only one relevant PDF ID
    int pdfid(-1);
    if (p_pdf1 != NULL) {
      pdfid = p_pdf1->LHEFNumber();
    } else if (p_pdf2 != NULL) {
      pdfid = p_pdf2->LHEFNumber();
    } else if (p_alphas->GetAs()->PDF() != NULL) {
      pdfid = p_alphas->GetAs()->PDF()->LHEFNumber();
    } else {
      // THROW(fatal_error, "Cannot obtain PDF IDF");
    }
    name = GenerateNamePart("MUR", sqrt(m_muR2fac)) + divider
           + GenerateNamePart("MUF", sqrt(m_muF2fac)) + divider
           + GenerateNamePart("PDF", pdfid);
  } else {
    // there are two relevant PDF IDs, quote both
    name = GenerateNamePart("MUR", sqrt(m_muR2fac)) + divider
           + GenerateNamePart("MUF", sqrt(m_muF2fac)) + divider
           + GenerateNamePart("PDF", p_pdf1->LHEFNumber()) + divider
           + GenerateNamePart("PDF", p_pdf2->LHEFNumber());
  }
  // append non-trivial shower scale factors
  if (m_showermuR2fac != 1.0 || m_showermuF2fac != 1.0) {
    name += divider + GenerateNamePart("PSMUR", sqrt(m_showermuR2fac));
    name += divider + GenerateNamePart("PSMUF", sqrt(m_showermuF2fac));
  }
  // append non-trivial shower scale factors
  if (m_Qcutfac != 1.0) {
    name += divider + GenerateNamePart("QCUT",m_Qcutfac);
  }
  return name;
}


template <typename U>
std::string Variation_Parameters::GenerateNamePart(std::string tag, U value) const
{
  return tag + ToString(value);
}


Subevent_Weights_Vector::Subevent_Weights_Vector():
  std::vector<double>()
{}


Subevent_Weights_Vector::Subevent_Weights_Vector(size_type count, const double& value):
  std::vector<double>(count, value)
{};


Subevent_Weights_Vector &
Subevent_Weights_Vector::operator*=(const double &scalefactor)
{
  for (iterator it(begin()); it != end(); ++it) {
    *it *= scalefactor;
  }
  return *this;
}

Subevent_Weights_Vector &
Subevent_Weights_Vector::operator*=(const Subevent_Weights_Vector &other)
{
  if (size() == other.size()) {
    for (size_t i(0); i < size(); ++i) {
      (*this)[i] *= other[i];
    }
  } else if (other.size() == 1) {
    *this *= other[0];
  }
  return *this;
}

Subevent_Weights_Vector &
Subevent_Weights_Vector::operator+=(const Subevent_Weights_Vector &other)
{
  if (size() != other.size()) {
    THROW(fatal_error, "Can not add subevent weights of different size.");
  }
  for (size_t i(0); i < size(); i++) {
    (*this)[i] *= other[i];
  }
  return *this;
}


void Variation_Weights::Reset()
{
  m_weights.clear();
}


Variation_Weights & Variation_Weights::operator*=(const double &scalefactor)
{
  if (GetNumberOfVariations() == 0) {
    return *this;
  }
  if (!AreWeightsInitialised()) {
    THROW(fatal_error, "Can not multiply uninitialised variation weights.");
  }
  for (auto& w : m_weights[Variations_Type::main])
    w *= scalefactor;
  return *this;
}


Variation_Weights & Variation_Weights::operator*=(const Variation_Weights &other)
{
  if (GetNumberOfVariations() == 0)
    return *this;
  if (!AreWeightsInitialised())
    THROW(fatal_error, "Can not multiply uninitialised variation weights.");
  if (other.m_weights.size() == 0)
    return *this;
  if (other.m_weights.size() > 1)
    THROW(fatal_error, "Refuse to multiply variation weights with more than one variation type.");
  for (Variations::Parameters_Vector::size_type i(0);
       i < GetNumberOfVariations();
       ++i) {
    this->m_weights[Variations_Type::main][i] *= other.GetVariationWeightAt(i);
  }
  return *this;
}

Variation_Weights & Variation_Weights::operator+=(const Variation_Weights &other)
{
  if (GetNumberOfVariations() == 0) {
    return *this;
  }
  if (!other.AreWeightsInitialised()) {
    return *this;
  }
  if (!AreWeightsInitialised()) {
    InitialiseWeights(Subevent_Weights_Vector(other.GetNumberOfSubevents(), 0.0), Variations_Type::main);
  } else if (GetNumberOfSubevents() != other.GetNumberOfSubevents()) {
    THROW(fatal_error, "Can not add variation weights with differing numbers of subevents.");
  }
  if (GetNumberOfVariations() != other.GetNumberOfVariations()) {
    THROW(fatal_error, "Can not add variation weights with differing numbers of variations.");
  }
  for (Variations::Parameters_Vector::size_type i(0);
       i < GetNumberOfVariations();
       ++i) {
    for (Subevent_Weights_Vector::size_type j(0);
	 j < GetNumberOfSubevents();
	 ++j) {
      this->m_weights[Variations_Type::main][i][j] += other.GetVariationWeightAt(i, Variations_Type::main, j);
    }
  }
  return *this;
}

void Variation_Weights::CombineSubeventWeights()
{
  if (!AreWeightsInitialised()) return;
  if (m_weights.size() > 1)
    THROW(fatal_error, std::string("Refuse to combine subevent weights when")
                       + " more than one variation type is stored.");
  auto& weights = m_weights.begin()->second;  // select the only type
  for (auto& subeventweights : weights) {
    const auto combinedweight =
      std::accumulate(subeventweights.begin(), subeventweights.end(), 0.0);
    subeventweights.clear();
    subeventweights.push_back(combinedweight);
  }
}

std::string Variation_Weights::GetVariationNameAt(Variations::Parameters_Vector::size_type i) const
{
  return p_variations->GetParametersVector()->at(i)->m_name;
}


double Variation_Weights::GetVariationWeightAt(
    Variations::Parameters_Vector::size_type paramidx,
    Variations_Type t,
    int subevtidx) const
{
  if (subevtidx < 0) {
    Subevent_Weights_Vector weights(GetNumberOfSubevents());
    for (const auto kv : m_weights)
      if (t == Variations_Type::all || t == kv.first)
        weights *= kv.second[paramidx];
    return std::accumulate(weights.begin(), weights.end(), 0.0);
  } else {
    double weight(1.0);
    for (const auto kv : m_weights) {
      if (t == Variations_Type::all || t == kv.first) {
        if (subevtidx > 0 && kv.second[paramidx].size() == 1) {
          if (kv.first == Variations_Type::main) {
            THROW(fatal_error,
                "The main variation weights do not have enough entries.");
          }
          weight *= kv.second[paramidx][0];
        } else {
          weight *= kv.second[paramidx][subevtidx];
        }
      }
    }
    return weight;
  }
}


Variations::Parameters_Vector::size_type Variation_Weights::CurrentParametersIndex() const
{
  if (!m_reweighting) THROW(fatal_error, "There is no ongoing reweighting.");
  return m_currentparametersindex;
}


size_t Variation_Weights::GetNumberOfVariations() const
{
  std::map<Variations_Type, std::vector<Subevent_Weights_Vector>>::const_iterator it{
    m_weights.find(Variations_Type::main) };
  if (it == m_weights.end())
    return 0;
  return it->second.size();
}


size_t Variation_Weights::GetNumberOfSubevents() const
{
  std::map<Variations_Type, std::vector<Subevent_Weights_Vector>>::const_iterator it{
    m_weights.find(Variations_Type::main) };
  if (it == m_weights.end())
    return 0;
  return it->second[0].size();
}


void Variation_Weights::InitialiseWeights(const Subevent_Weights_Vector & subweights,
                                          const Variations_Type t)
{
  const size_t size(p_variations->GetParametersVector()->size());
  m_weights[t].clear();
  m_weights[t].reserve(size);
  for (size_t i(0); i < size; ++i) {
    m_weights[t].push_back(subweights);
  }
}


bool Variation_Weights::AreWeightsInitialised(
    const Variations_Type t) const
{
  return (m_weights.find(t) != m_weights.end());
}


namespace ATOOLS {

  std::ostream& operator<<(std::ostream& s, const Variations& v)
  {
    const Variations::Parameters_Vector * const paramsvec(v.GetParametersVector());
    s << "Named variations:" << std::endl;
    for (Variations::Parameters_Vector::const_iterator it(paramsvec->begin());
         it != paramsvec->end(); ++it) {
      s << (*it)->m_name << " (" << (*it)->m_deletepdfsandalphas << ")" << std::endl;
    }
    return s;
  }

  std::ostream& operator<<(std::ostream& s, const Subevent_Weights_Vector& v)
  {
    if (v.size() == 1) {
      s << v[0];
    } else {
      s << "(";
      for (size_t j{ 0 }; j < v.size(); ++j) {
        if (j != 0)
          s << ", ";
        s << v[j];
      }
      s << ")";
    }
    return s;
  }

  std::ostream& operator<<(std::ostream& s, const Variations_Type &c)
   {
    switch (c) {
      case Variations_Type::all:
        return s << "All";
      case Variations_Type::main:
        return s << "Main";
      case Variations_Type::sudakov:
        return s << "Sudakov";
    }
  }

  std::ostream& operator<<(std::ostream& s, const Variation_Weights& weights)
  {
    const Variations::Parameters_Vector * const paramsvec(
        weights.p_variations->GetParametersVector());
    s << "Variation weights: {" << std::endl;
    for (Variations::Parameters_Vector::size_type i(0);
         i < paramsvec->size(); ++i) {
      s << "    " << (*paramsvec)[i]->m_name << ": ";
      for (auto& kv : weights.m_weights)
        s << kv.first << "=" << kv.second[i] << " ";
      s << std::endl;
    }
    s << "}";
    return s;
  }

  // Explicit template instantiations
  template <> Blob_Data<Variation_Weights>::~Blob_Data() {}
  template class Blob_Data<Variation_Weights>;

}