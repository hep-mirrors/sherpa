#include "SHERPA/Tools/Variations.H"

#include <iterator>

#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Smart_Pointer.C"
#include "ATOOLS/Phys/Blob.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "PDF/Main/PDF_Base.H"
#if defined USING__LHAPDF && defined USING__LHAPDF6
#include "LHAPDF/LHAPDF.h"
#endif

using namespace ATOOLS;
using namespace SHERPA;

typedef void (*PDF_Init_Function)();

Variations::Variations(Data_Reader * const reader)
{
#if defined USING__LHAPDF && defined USING__LHAPDF6
  LoadLHAPDFInterfaceIfNecessary(reader);
  const int lhapdfverbosity(LHAPDF::verbosity());
  LHAPDF::setVerbosity(0);
#endif

  // read settings that are relevant for parsing the prompted variations
  m_reweightsplittingalphasscales = reader->GetValue<int>("REWEIGHT_SPLITTING_ALPHAS_SCALES", 0);
  m_reweightsplittingpdfsscales = reader->GetValue<int>("REWEIGHT_SPLITTING_PDF_SCALES", 0);
  if (m_reweightsplittingpdfsscales) {
    THROW(not_implemented, "PDF scale factors in shower splittings is not implemented yet.");
  }

  InitialiseParametersVector(reader);

  if (!m_parameters_vector.empty()) {
    rpa->gen.AddCitation(1, "The Sherpa-internal reweighting is published in \\cite{Bothmann:2016nao}.");
  }

#if defined USING__LHAPDF && defined USING__LHAPDF6
  LHAPDF::setVerbosity(lhapdfverbosity);
#endif
}


Variations::~Variations()
{
#if VARIATION_WARNINGS
  CollectAndOutputWarnings();
#endif

  for (Parameters_Vector::const_iterator it = m_parameters_vector.begin();
       it != m_parameters_vector.end();
       ++it) {
    delete *it;
  }
}


void Variations::CollectAndOutputWarnings()
{
  PRINT_FUNC("Variation warning statistics");

  // initialise with own warnings, but note that parameter-specific warnings
  // will also be taken into account
  bool emptystats(m_warnings.size() == 0);

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
        emptystats = false;
        warningnums[warningsit->first] += warningsit->second;
        if (warningreps[warningsit->first].second < warningsit->second) {
          warningreps[warningsit->first].first  = (*paramsit)->m_name;
          warningreps[warningsit->first].second = warningsit->second;
        }
      }
    }
  }

  if (emptystats) {
    msg_Info() << "No warnings reported." << std::endl;
  } else {
    msg_Info() << "Total number of warnings:" << std::endl;
    // output own warnings (warning -> number of reports)
    for (std::map<std::string, unsigned long>::const_iterator
         it(m_warnings.begin()); it != m_warnings.end(); ++it) {
      msg_Info() << it->first << ": " << it->second << std::endl;
    }
    msg_Info() << "Total number of warnings (parameter-specific):" << std::endl;
    // output parameter-specific warnings (warning -> number of reports)
    for (std::map<std::string, unsigned long>::const_iterator it(warningnums.begin());
         it != warningnums.end(); ++it) {
      msg_Info() << it->first << ": " << it->second << std::endl;
    }
    // output parameter-specific warnings
    // (warning -> (reporter that reported largest number of reports))
    msg_Info() << "Parameter variation that reported most warnings:" << std::endl;
    for (std::map<std::string, std::pair<std::string, unsigned long> >::const_iterator
         it(warningreps.begin());
         it != warningreps.end(); ++it) {
      msg_Info() << it->first << ": " << it->second.first;
      msg_Info() << " (" << it->second.second << ")" << std::endl;
    }
  }
}

void Variations::IncrementOrInitialiseWarningCounter(const std::string name)
{
#if VARIATION_WARNINGS
  m_warnings[name]++;
#endif
}

#if defined USING__LHAPDF && defined USING__LHAPDF6
void Variations::LoadLHAPDFInterfaceIfNecessary(Data_Reader * const reader)
{
  // check whether LHAPDF is already loaded, if not load and init interface
  if (!s_loader->LibraryIsLoaded("LHAPDFSherpa")) {
    s_loader->AddPath(std::string(LHAPDF_PATH)+"/lib");
    s_loader->LoadLibrary("LHAPDF");
    void *init(s_loader->GetLibraryFunction("LHAPDFSherpa","InitPDFLib"));
    if (init==NULL) THROW(fatal_error,"Cannot load PDF library LHAPDFSherpa");
    ((PDF_Init_Function)init)();
  }
  std::string path;
  if (reader->ReadFromFile(path,"LHAPDF_GRID_PATH")) LHAPDF::setPaths(path);
}
#endif


void Variations::InitialiseParametersVector(Data_Reader * const reader)
{
  std::vector<std::string> args = VariationArguments(reader);
  PRINT_FUNC(args.size() << " variations");
  for (std::vector<std::string>::const_iterator it(args.begin());
      it != args.end(); ++it) {
    std::vector<std::string> params(VariationArgumentParameters(*it));
    AddParameters(params, reader);
  }
  msg_Info() << *this;
}


std::vector<std::string> Variations::VariationArguments(Data_Reader * const reader)
{
  std::vector<std::string> varargs, pdfvarargs;
  varargs = VariationArguments(reader, "SCALE_VARIATIONS");
  // the PDF_VARIATIONS has a more specific syntax that must be transformed to
  // the more general syntax of SCALE_VARIATIONS
  pdfvarargs = VariationArguments(reader, "PDF_VARIATIONS");
  for (size_t i(0); i < pdfvarargs.size(); ++i) {
    varargs.push_back("1.,1.," + pdfvarargs[i]);
  }
  return varargs;
}


std::vector<std::string> Variations::VariationArguments(Data_Reader * const reader,
                                                        std::string tag)
{
  std::vector<std::string> args;
  reader->VectorFromFile(args, tag);
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
  // parse ME and PS scale factors
  const double muR2fac(ToType<double>(stringparams[0]));
  const double muF2fac(ToType<double>(stringparams[1]));
  const double showermuR2fac = (m_reweightsplittingalphasscales) ? muR2fac : 1.0;

  // parse PDF member(s)
  std::vector<PDFs_And_AlphaS> pdfsandalphasvector;
  bool ownedpdfsandalphas(false);
  if (stringparams.size() > 2) {
    // PDF variation requested
    ownedpdfsandalphas = true;
    std::string pdfname(stringparams[2]);
    if (pdfname.find("[all]") == std::string::npos) {
      // single PDF member: "Set/i" or just "Set"
      int member(0);
      if (pdfname.find("/") != std::string::npos) {
        member = ToType<int>(std::string(pdfname, pdfname.find("/") + 1));
        pdfname = pdfname.substr(0, pdfname.find("/"));
      }
      pdfsandalphasvector.push_back(PDFs_And_AlphaS(reader, pdfname, member));
    } else {
      // all PDF members: "Set[all]"
#if defined USING__LHAPDF && defined USING__LHAPDF6
      // check whether interface is loaded
      if (!s_loader->LibraryIsLoaded("LHAPDFSherpa"))
        THROW(fatal_error, "LHAPDF interface not initialised."
              + std::string(" Add LHAPDFSherpa to PDF_LIBRARY"));
      // assume members are labeled 0..n
      pdfname=pdfname.substr(0, pdfname.find("[all]"));
      const std::vector<std::string>& availablepdfsets(LHAPDF::availablePDFSets());
      if (std::find(availablepdfsets.begin(), availablepdfsets.end(), pdfname)
          == availablepdfsets.end()) {
        THROW(fatal_error, "PDF set " + pdfname + " not available.");
      }
      LHAPDF::PDFSet set(pdfname);
      for (size_t j(0); j < set.size(); ++j) {
        pdfsandalphasvector.push_back(PDFs_And_AlphaS(reader, pdfname, j));
      }
#else
      THROW(not_implemented,"Full set reweightings only work with LHAPDF6."
          + std::string(" Otherwise specify separately."));
#endif
    }
  } else {
    pdfsandalphasvector.push_back(PDFs_And_AlphaS());
  }

  for (std::vector<PDFs_And_AlphaS>::const_iterator pdfasit(pdfsandalphasvector.begin());
        pdfasit != pdfsandalphasvector.end(); pdfasit++) {
    Variation_Parameters *params =
      new Variation_Parameters(muR2fac, muF2fac,
          showermuR2fac,
          pdfasit->m_pdfs[0],
          pdfasit->m_pdfs[1],
          pdfasit->p_alphas,
          ownedpdfsandalphas);
    m_parameters_vector.push_back(params);
  }
}


Variations::PDFs_And_AlphaS::PDFs_And_AlphaS():
    p_alphas(MODEL::as->GetAs(PDF::isr::hard_process))
{
  // Workaround für C++03 (vector constructor is confused when fed
  // with NULL as the initial value for a pointer)
  // cf. https://gcc.gnu.org/ml/gcc-help/2013-02/msg00026.html
  PDF::PDF_Base *nullPtr = NULL;
  m_pdfs = std::vector<PDF::PDF_Base *>(2, nullPtr);
  m_pdfs[0] = rpa->gen.PDF(0);
  m_pdfs[1] = rpa->gen.PDF(1);
}


Variations::PDFs_And_AlphaS::PDFs_And_AlphaS(Data_Reader * const reader,
    std::string pdfname, size_t pdfmember)
{
  // obtain PDFs
  PDF::PDF_Base *aspdf(NULL);
  for (int i(0); i < 2; ++i) {
    if (rpa->gen.Bunch(i).IsHadron()) {
      PDF::PDF_Arguments args(rpa->gen.Bunch(i), reader, i, pdfname, pdfmember);
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
    p_alphas = new MODEL::One_Running_AlphaS(pdfname, pdfmember);
  } else {
    p_alphas = new MODEL::One_Running_AlphaS(aspdf);
  }
  if (p_alphas == NULL) {
    THROW(fatal_error, "AlphaS for " + pdfname + " could not be initialised.");
  }
}


Variation_Parameters::~Variation_Parameters()
{
  if (m_deletepdfsandalphas) {
    if (p_pdf1) { delete p_pdf1; }
    if (p_pdf2) { delete p_pdf2; }
    if (p_alphas) { delete p_alphas; }
  }
}

void Variation_Parameters::IncrementOrInitialiseWarningCounter(const std::string name)
{
#if VARIATION_WARNINGS
  m_warnings[name]++;
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
    } else if (p_alphas->PDF() != NULL) {
      pdfid = p_alphas->PDF()->LHEFNumber();
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
  if (m_showermuR2fac != 1.0) {
    name += divider + GenerateNamePart("PSMUR", sqrt(m_showermuR2fac));
    name += divider + GenerateNamePart("PSMUF", 1.0);
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


namespace SHERPA {

  std::ostream& operator<<(std::ostream &s, const Variations &v)
  {
    const Variations::Parameters_Vector * const paramsvec(v.GetParametersVector());
    s << "Named variations:" << std::endl;
    for (Variations::Parameters_Vector::const_iterator it(paramsvec->begin());
         it != paramsvec->end(); ++it) {
      s << (*it)->m_name << std::endl;
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
}


namespace ATOOLS {

  // Explicit template instantiations
  template <> Blob_Data<Variation_Weights>::~Blob_Data() {}
  template class Blob_Data<Variation_Weights>;
  template class SP(Variation_Weights);

}
