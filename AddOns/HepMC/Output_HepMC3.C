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

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Output_HepMC3::Output_HepMC3(const Output_Arguments &args) :
  Output_Base("HepMC3")
{
  m_basename = args.m_outpath + "/" + args.m_outfile;
  m_iotype = Settings::GetMainSettings()["HEPMC3_IO_TYPE"].SetDefault(0).Get<int>();
  int precision = Settings::GetMainSettings()["HEPMC3_OUTPUT_PRECISION"].SetDefault(12).Get<int>();
  m_short = Settings::GetMainSettings()["HEPMC3_SHORT"].SetDefault(false).Get<bool>();
#ifdef USING__GZIP
  m_ext += ".gz";
#endif
#ifdef USING__MPI
  if (mpi->Size()>1) {
    m_basename += "_"+rpa->gen.Variable("RNG_SEED");
  }
#endif

switch (m_iotype)
    {
    case 0:
    {
        m_outstream.open((m_basename+m_ext).c_str());
        if (!m_outstream.good())THROW(fatal_error, "Could not open event file "+m_basename+m_ext+".");
        HepMC::WriterAscii* t_writer=new HepMC::WriterAscii(m_outstream);
        t_writer->set_precision(precision);
        p_writer=t_writer;
    }
    break;
    case 1:
        m_outstream.open((m_basename+m_ext).c_str());
        if (!m_outstream.good())THROW(fatal_error, "Could not open event file "+m_basename+m_ext+".");
        p_writer=new HepMC::WriterHEPEVT(m_outstream);
        break;
    case 2:
    {
        m_outstream.open((m_basename+m_ext).c_str());
        if (!m_outstream.good())THROW(fatal_error, "Could not open event file "+m_basename+m_ext+".");
        HepMC::WriterAsciiHepMC2* t_writer=new HepMC::WriterAsciiHepMC2(m_outstream);
        t_writer->set_precision(precision);
        p_writer=t_writer;
    }
    break;
    case 3:
#ifdef USING__HEPMC3__ROOT
        p_writer=new HepMC::WriterRoot(m_basename);
#else
        THROW(fatal_error,"Asked for Root output, but Sherpa/HepMC3 was compiled without Root output support.");
#endif
        break;
    case 4:
#ifdef USING__HEPMC3__ROOT
        p_writer=new HepMC::WriterRootTree(m_basename);
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

