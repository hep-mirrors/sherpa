#include "AddOns/HepMC/Output_HepMC3.H"
#include "HepMC3/GenEvent.h"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include "HepMC3/Writer.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3/WriterHEPEVT.h"
#ifdef USING__HEPMC3__ROOT
#include "HepMC3/WriterRootTree.h"
#include "HepMC3/WriterRoot.h"
#endif
/** If Sherpa globaly uses GZIP, enable HepMC3 compression with ZLIB support */
#if defined(USING__GZIP)
#define HEPMC3_USE_COMPRESSION 1
#define HEPMC3_Z_SUPPORT 1
#endif
#include "HepMC3/WriterGZ.h"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

template <class T> void set_writer_precision(T& writer, const int precision) { writer->set_precision(precision);}
template <> void set_writer_precision(std::shared_ptr<HepMC3::WriterHEPEVT>& writer, const int precision) { }
template <class T>
std::shared_ptr<HepMC3::Writer> create_writer(std::ofstream& outstream, const std::string &use_compression, const int precision,  const std::string& basename,  std::string & ext) {
#if HEPMC3_USE_COMPRESSION
#if HEPMC3_Z_SUPPORT
    if (use_compression == "GZ" )   {
        ext += ".gz";
        outstream.open((basename + ext).c_str());
        if (!outstream.good()) THROW(fatal_error, "Could not open event file " + basename + ext + ".");
        auto X = std::make_shared< HepMC3::WriterGZ<T,HepMC3::Compression::z> >(outstream);
        set_writer_precision(X->writer(), precision);
        return X;
    }
#endif
#if HEPMC3_LZMA_SUPPORT
    if (use_compression == "LZMA" ) {
        ext += ".lz";
        outstream.open((basename + ext).c_str());
        if (!outstream.good()) THROW(fatal_error, "Could not open event file " + basename + ext + ".");
        auto X = std::make_shared< HepMC3::WriterGZ<T,HepMC3::Compression::lzma> >(outstream);
        set_writer_precision(X->writer(), precision);
        return X;
    }
#endif
#if HEPMC3_BZ2_SUPPORT
    if (use_compression == "BZ2" )  {
        ext += ".bz2";
        outstream.open((basename + ext).c_str());
        if (!outstream.good()) THROW(fatal_error, "Could not open event file " + basename + ext + ".");
        auto X = std::make_shared< HepMC3::WriterGZ<T,HepMC3::Compression::bz2> >(outstream);
        set_writer_precision(X->writer(), precision);
        return X;
    }
#endif
#endif
    if (use_compression == "" )  {
       outstream.open((basename + ext).c_str());
       if (!outstream.good()) THROW(fatal_error, "Could not open event file " + basename + ext + ".");
       auto X = std::make_shared<T>(outstream);
       set_writer_precision(X, precision);
       return X;
    }
    THROW(fatal_error, "Usupported HEPMC3_COMPRESSION type " + use_compression + ".");
    return nullptr;
}

Output_HepMC3::Output_HepMC3(const Output_Arguments &args) :
  Output_Base("HepMC3")
{
  m_basename = args.m_outpath + "/" + args.m_outfile;
/** If compression is enabled and ZLIB type is enabled, that should be the default compression */
#if HEPMC3_USE_COMPRESSION && HEPMC3_Z_SUPPORT
  m_compression = Settings::GetMainSettings()["HEPMC3_COMPRESSION"].SetDefault("GZ").Get<std::string>();
#else
  m_compression = Settings::GetMainSettings()["HEPMC3_COMPRESSION"].SetDefault("").Get<std::string>();
#endif
  m_iotype = Settings::GetMainSettings()["HEPMC3_IO_TYPE"].SetDefault(0).Get<int>();
  int precision = Settings::GetMainSettings()["HEPMC3_OUTPUT_PRECISION"].SetDefault(12).Get<int>();
  m_short = Settings::GetMainSettings()["HEPMC3_SHORT"].SetDefault(false).Get<bool>();
#ifdef USING__MPI
  if (mpi->Size()>1) {
    m_basename += "_"+rpa->gen.Variable("RNG_SEED");
  }
#endif

switch (m_iotype)
    {
    case 0:
        p_writer = create_writer<HepMC3::WriterAscii>(m_outstream, m_compression, precision, m_basename, m_ext);
        break;
    case 1:
        p_writer = create_writer<HepMC3::WriterHEPEVT>(m_outstream, m_compression, precision, m_basename, m_ext);
        break;
    case 2:
        p_writer = create_writer<HepMC3::WriterAsciiHepMC2>(m_outstream, m_compression, precision, m_basename, m_ext);
        break;
    case 3:
#ifdef USING__HEPMC3__ROOT
        p_writer = std::make_shared<HepMC3::WriterRoot>(m_basename);
#else
        THROW(fatal_error,"Asked for Root output, but Sherpa/HepMC3 was compiled without Root output support.");
#endif
        break;
    case 4:
#ifdef USING__HEPMC3__ROOT
        p_writer = std::make_shared<HepMC3::WriterRootTree>(m_basename);
#else
        THROW(fatal_error,"Asked for RootTree output, but Sherpa/HepMC3 was compiled without RootTree output support.");
#endif
        break;
    default:
        THROW(fatal_error, "Output format HEPMC3_IO_TYPE is undefined.");
        break;
    }
  m_run_info= std::make_shared<HepMC::GenRunInfo>();
  HepMC::GenRunInfo::ToolInfo tool;
  tool.name = std::string("SHERPA-MC");
  tool.version = std::string(SHERPA_VERSION)+"."+std::string(SHERPA_SUBVERSION);
  tool.description = std::string(SHERPA_NAME);
  m_run_info->tools().push_back(tool);
}

Output_HepMC3::~Output_HepMC3()
{
  p_writer->close();
  m_outstream.close();
}

void Output_HepMC3::SetXS(const ATOOLS::Weights_Map& xs,
			        const ATOOLS::Weights_Map& xserr)
{
  // Only copy for now, we have to wait until the event weights have been
  // added (when Output()), otherwise HepMC3::GenCrossSection will not be
  // initialised correctly.
  m_xs = xs;
  m_err = xserr;
}

void Output_HepMC3::Output(Blob_List* blobs)
{
  if (m_short) {
    m_hepmc3.Sherpa2ShortHepMC(blobs, m_run_info);
  } else {
    m_hepmc3.Sherpa2HepMC(blobs, m_run_info);
  }
  HepMC::GenEvent* q=m_hepmc3.GenEvent();
  if (q)  m_hepmc3.AddCrossSection(*q, m_xs, m_err);
  std::vector<HepMC::GenEvent*> subevents(m_hepmc3.GenSubEventList());
  for (size_t i = 0; i<subevents.size(); ++i) {
    m_hepmc3.AddCrossSection(*subevents[i], m_xs, m_err);
  }
  if (subevents.size()) {
    for (size_t i = 0; i<subevents.size(); ++i) {
      if (p_writer)    p_writer->write_event(*subevents[i]);
    }
  }
  else if (q) {
    if (p_writer)  p_writer->write_event(*(q));
  }
}

void Output_HepMC3::ChangeFile()
{
  /*This should be implemented in HepMC3 library.*/
}

DECLARE_GETTER(Output_HepMC3,"HepMC3", Output_Base,Output_Arguments);

Output_Base *ATOOLS::Getter<Output_Base,Output_Arguments,Output_HepMC3>::
operator()(const Output_Arguments &args) const
{
  return new Output_HepMC3(args);
}

void ATOOLS::Getter<Output_Base,Output_Arguments,Output_HepMC3>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"HepMC3 output";
}

DECLARE_GETTER(Output_HepMC3_GenEvent,"HepMC3_GenEvent", Output_Base,Output_Arguments);

Output_Base *ATOOLS::Getter<Output_Base,Output_Arguments,Output_HepMC3_GenEvent>::
operator()(const Output_Arguments &args) const
{
  return new Output_HepMC3_GenEvent(args);
}

void ATOOLS::Getter<Output_Base,Output_Arguments,Output_HepMC3_GenEvent>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"HepMC3 GenEvent output";
}

DECLARE_GETTER(Output_HepMC3_Short,"HepMC3_Short", Output_Base,Output_Arguments);

Output_Base *ATOOLS::Getter<Output_Base,Output_Arguments,Output_HepMC3_Short>::
operator()(const Output_Arguments &args) const
{
  return new Output_HepMC3_Short(args);
}

void ATOOLS::Getter<Output_Base,Output_Arguments,Output_HepMC3_Short>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"HepMC3 Short output";
}

