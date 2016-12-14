#include "ATOOLS/Phys/Variations.H"

#include <iterator>
#include <numeric>

#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Smart_Pointer.C"
#include "ATOOLS/Phys/Blob.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/PDF_Base.H"
#if defined USING__LHAPDF && defined USING__LHAPDF6
#include "LHAPDF/LHAPDF.h"
#endif

using namespace ATOOLS;
using namespace SHERPA;

bool Variations::NeedsLHAPDF6Interface(std::string inputpath)
{
  // set up data reader
  Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.AddWordSeparator("\t");
  reader.SetInputPath(inputpath);
  // go through all parameters and return true if any non-double points to a
  // request for a PDF set/member
  std::vector<std::string> args = VariationArguments(&reader);
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

Variations::Variations(Data_Reader * const reader):
  m_includecentralvaluevariation(reader->GetValue<int>("VARIATIONS_INCLUDE_CV", 1)),
  m_reweightsplittingalphasscales(reader->GetValue<int>("REWEIGHT_SPLITTING_ALPHAS_SCALES", 0)),
  m_reweightsplittingpdfsscales(reader->GetValue<int>("REWEIGHT_SPLITTING_PDF_SCALES", 0))
{
#if defined USING__LHAPDF && defined USING__LHAPDF6
  int lhapdfverbosity(0);
  const bool needslhapdf(NeedsLHAPDF6Interface(reader->InputPath()));
  if (needslhapdf) {
    if (!s_loader->LibraryIsLoaded("LHAPDFSherpa")) {
      THROW(fatal_error, "LHAPDF interface is not initialised.");
    }
    lhapdfverbosity = LHAPDF::verbosity();
    LHAPDF::setVerbosity(0);
  }
#endif

  InitialiseParametersVector(reader);

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


void Variations::InitialiseParametersVector(Data_Reader * const reader)
{
  std::vector<std::string> args = VariationArguments(reader);
  for (std::vector<std::string>::const_iterator it(args.begin());
      it != args.end(); ++it) {
    std::vector<std::string> params(VariationArgumentParameters(*it));
    AddParameters(params, reader);
  }
}


std::vector<std::string> Variations::VariationArguments(Data_Reader * const reader)
{
  std::vector<std::string> args;
  reader->VectorFromFile(args, "VARIATIONS");
  if (args.size() == 1 && args[0] == "None") {
    args.clear();
  }
  return args;
}


std::vector<std::string> Variations::VariationArgumentParameters(std::string arg)
{
  const std::string delimiter = ",";
  std::vector<std::string> params;
  size_t pos = 0;
  while (true) {
    pos = arg.find(delimiter);
    params.push_back(arg.substr(0, pos));
    if (pos == std::string::npos) {
      break;
    }
    arg.erase(0, pos + delimiter.length());
  }
  return params;
}


void Variations::AddParameters(std::vector<std::string> stringparams,
    Data_Reader * const reader)
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
    if (reader->GetValue<int>("OVERRIDE_PDF_INFO",0)==1) {
      THROW(fatal_error, "OVERRIDE_PDF_INFO=1 is incompatible with doing PDF/AlphaS variations.");
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
  if (scalevariationrequested) {
    if (scalestringparams.size() == 1) {
      muR2fac = muF2fac = ToType<double>(scalestringparams[0]);
    } else {
      muR2fac = ToType<double>(scalestringparams[0]);
      muF2fac = ToType<double>(scalestringparams[1]);
    }
  }

  // expand scale factors
  std::vector<std::pair<double, double> > scalefactorpairs;
  std::pair<double, double> defaultscalefactorpair(1.0, 1.0);

  if (!scalevariationrequested || expansions == ScaleFactorExpansions::None) {
    defaultscalefactorpair = std::make_pair(muR2fac, muF2fac);
    scalefactorpairs.push_back(defaultscalefactorpair);
  } else if (expansions & ScaleFactorExpansions::SevenPoint) {
    scalefactorpairs.push_back(std::pair<double, double>(muR2fac, 1.0));
    scalefactorpairs.push_back(std::pair<double, double>(1.0, muF2fac));
    scalefactorpairs.push_back(std::pair<double, double>(muR2fac, muF2fac));
    scalefactorpairs.push_back(std::pair<double, double>(1.0 / muR2fac, 1.0));
    scalefactorpairs.push_back(std::pair<double, double>(1.0, 1.0 / muF2fac));
    scalefactorpairs.push_back(std::pair<double, double>(1.0 / muR2fac, 1.0 / muF2fac));
  } else {
    if (expansions & ScaleFactorExpansions::First) {
      const double defaultmuF2fac = (expansions & ScaleFactorExpansions::Second) ? 1.0 : muF2fac;
      defaultscalefactorpair.second = defaultmuF2fac;
      scalefactorpairs.push_back(std::pair<double, double>(muR2fac, defaultmuF2fac));
      scalefactorpairs.push_back(std::pair<double, double>(1.0 / muR2fac, defaultmuF2fac));
    }
    if (expansions & ScaleFactorExpansions::Second) {
      const double defaultmuR2fac = (expansions & ScaleFactorExpansions::First) ? 1.0 : muR2fac;
      defaultscalefactorpair.first = defaultmuR2fac;
      scalefactorpairs.push_back(std::pair<double, double>(defaultmuR2fac, muF2fac));
      scalefactorpairs.push_back(std::pair<double, double>(defaultmuR2fac, 1.0 / muF2fac));
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
      AddParameters(defaultscalefactorpair.first, defaultscalefactorpair.second, pdfasit,
                    !assignedownershipofpdfsandalphas);
      assignedownershipofpdfsandalphas = true;
    }

    if (pdfasit == pdfsandalphasvector.begin()) {
      for (std::vector<std::pair<double, double> >::const_iterator scalefactorpairit(scalefactorpairs.begin());
          scalefactorpairit != scalefactorpairs.end(); scalefactorpairit++) {
        AddParameters(scalefactorpairit->first, scalefactorpairit->second, pdfasit,
                      !assignedownershipofpdfsandalphas);
        assignedownershipofpdfsandalphas = true;
      }
    }
  }
}


void Variations::AddParameters(double muR2fac, double muF2fac,
                               std::vector<PDFs_And_AlphaS>::const_iterator pdfsandalphas,
                               bool deletepdfsandalphas)
{
  const double showermuR2fac = (m_reweightsplittingalphasscales) ? muR2fac : 1.0;
  const double showermuF2fac = (m_reweightsplittingpdfsscales) ? muF2fac : 1.0;
  Variation_Parameters *params =
    new Variation_Parameters(
        muR2fac, muF2fac, showermuR2fac, showermuF2fac,
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


Variations::PDFs_And_AlphaS::PDFs_And_AlphaS(std::string pdfname, size_t pdfmember)
{
  Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.AddWordSeparator("\t");
  // obtain PDFs
  PDF::PDF_Base *aspdf(NULL);
  for (int i(0); i < 2; ++i) {
    if (rpa->gen.Bunch(i).IsHadron()) {
      PDF::PDF_Arguments args(rpa->gen.Bunch(i), &reader, i, pdfname, pdfmember);
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
  // open file
  Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.AddWordSeparator("\t");
  std::string outpath = reader.GetValue<std::string>("ANALYSIS_OUTPUT","Analysis/");
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
      THROW(fatal_error, "Cannot obtain PDF IDF");
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


void Variation_Weights::Reset()
{
  m_weights.clear();
  m_initialised = false;
}


Variation_Weights & Variation_Weights::operator*=(const double &scalefactor)
{
  if (GetNumberOfVariations() == 0) {
    return *this;
  }
  if (!m_initialised) {
    THROW(fatal_error, "Can not multiply uninitialised variation weights.");
  }
  for (std::vector<Subevent_Weights_Vector>::iterator it(m_weights.begin());
       it != m_weights.end(); ++it) {
    *it *= scalefactor;
  }
  return *this;
}


Variation_Weights & Variation_Weights::operator*=(const Variation_Weights &other)
{
  if (GetNumberOfVariations() == 0) {
    return *this;
  }
  if (!m_initialised) {
    THROW(fatal_error, "Can not multiply uninitialised variation weights.");
  }
  if (!other.m_initialised) {
    return *this;
  }
  for (Variations::Parameters_Vector::size_type i(0);
       i < GetNumberOfVariations();
       ++i) {
    this->m_weights[i] *= other.GetVariationWeightAt(i);
  }
  return *this;
}


std::string Variation_Weights::GetVariationNameAt(Variations::Parameters_Vector::size_type i) const
{
  return p_variations->GetParametersVector()->at(i)->m_name;
}


double Variation_Weights::GetVariationWeightAt(Variations::Parameters_Vector::size_type paramidx,
                                               int subevtidx) const
{
  if (subevtidx < 0) {
    return std::accumulate(m_weights[paramidx].begin(), m_weights[paramidx].end(), 0.0);
  } else { 
    return m_weights[paramidx][subevtidx];
  }
}


Variations::Parameters_Vector::size_type Variation_Weights::CurrentParametersIndex() const
{
  if (!m_reweighting) THROW(fatal_error, "There is no ongoing reweighting.");
  return m_currentparametersindex;
}


void Variation_Weights::InitialiseWeights(const Subevent_Weights_Vector & subweights) {
  const size_t size(p_variations->GetParametersVector()->size());
  m_weights.clear();
  m_weights.reserve(size);
  for (size_t i(0); i < size; ++i) {
    m_weights.push_back(subweights);
  }
  m_initialised = true;
}


namespace ATOOLS {

  std::ostream& operator<<(std::ostream &s, const Variations &v)
  {
    const Variations::Parameters_Vector * const paramsvec(v.GetParametersVector());
    s << "Named variations:" << std::endl;
    for (Variations::Parameters_Vector::const_iterator it(paramsvec->begin());
         it != paramsvec->end(); ++it) {
      s << (*it)->m_name << " (" << (*it)->m_deletepdfsandalphas << ")" << std::endl;
    }
    return s;
  }


  std::ostream & operator<<(std::ostream & s, const Variation_Weights & weights)
  {
    const Variations::Parameters_Vector * const paramsvec(weights.p_variations->GetParametersVector());
    s << "Variation weights: {" << std::endl;
    for (Variations::Parameters_Vector::size_type i(0);
         i < paramsvec->size(); ++i) {
      s << "    " << (*paramsvec)[i]->m_name << ": ";
      if (!weights.m_initialised) {
        s << "not initialised";
      } else {
        s << weights.m_weights[i];
      }
      s << std::endl;
    }
    s << "}" << std::endl;
    return s;
  }

  // Explicit template instantiations
  template <> Blob_Data<Variation_Weights>::~Blob_Data() {}
  template class Blob_Data<Variation_Weights>;
  template class SP(Variation_Weights);

}
