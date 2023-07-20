#include "AddOns/HepMC/Output_HepMC3_Short.H"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenCrossSection.h"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"

#include "HepMC3/Writer.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3/WriterHEPEVT.h"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Output_HepMC3_Short::Output_HepMC3_Short(const Output_Arguments &args) :
  Output_Base("HepMC3S")
{
  m_basename=args.m_outpath+"/"+args.m_outfile;
  m_iotype=args.p_reader->GetValue<int>("HEPMC3_IO_TYPE",0);
  int precision= args.p_reader->GetValue<int>("HEPMC3_OUTPUT_PRECISION",12);
  #ifdef USING__MPI
  if (mpi->Size()>1) {
    m_basename+="_"+rpa->gen.Variable("RNG_SEED");
  }
  #endif
  switch (m_iotype) {
    case 0:
    {
        HepMC::WriterAscii* t_writer=new HepMC::WriterAscii( m_basename);
        t_writer->set_precision(precision);
        p_writer=t_writer;
    }
    break;
    case 1:
        p_writer=new HepMC::WriterHEPEVT( m_basename);
        break;
    case 2:
    {
        HepMC::WriterAsciiHepMC2* t_writer=new HepMC::WriterAsciiHepMC2( m_basename);
        t_writer->set_precision(precision);
        p_writer=t_writer;
    }
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

Output_HepMC3_Short::~Output_HepMC3_Short()
{
  p_writer->close();
}

void Output_HepMC3_Short::SetXS(const double& xs, const double& xserr)
{
  // Only copy for now, we have to wait until the event weights have been
  // added (when Output()), otherwise HepMC3::GenCrossSection will not be
  // initialised correctly.
  m_xs = xs;
  m_err = xserr;
}

void Output_HepMC3_Short::Output(Blob_List* blobs, const double weight)
{
  m_hepmc3.Sherpa2ShortHepMC(blobs, m_run_info);
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
    if (p_writer)    p_writer->write_event(*(q));
  }

}

void Output_HepMC3_Short::ChangeFile()
{
  /*This should be implemented in HepMC3 library.*/
}

DECLARE_GETTER(Output_HepMC3_Short,"HepMC3_Short",
	       Output_Base,Output_Arguments);

Output_Base *ATOOLS::Getter<Output_Base,Output_Arguments,Output_HepMC3_Short>::
operator()(const Output_Arguments &args) const
{
  return new Output_HepMC3_Short(args);
}

void ATOOLS::Getter<Output_Base,Output_Arguments,Output_HepMC3_Short>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"HepMC3 short output";
}

