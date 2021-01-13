#include "AddOns/HepMC/Output_HepMC3_Genevent.H"
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
#ifdef USING__HEPMC3__ROOT
#include "HepMC3/WriterRootTree.h"
#include "HepMC3/WriterRoot.h"
#endif

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Output_HepMC3_Genevent::Output_HepMC3_Genevent(const Output_Arguments &args) :
  Output_Base("HepMC3")
{
  m_basename=args.m_outpath+"/"+args.m_outfile;
  m_iotype=args.p_reader->GetValue<int>("HEPMC3_IO_TYPE",0);
  int precision= args.p_reader->GetValue<int>("HEPMC3_OUTPUT_PRECISION",12);
#ifdef USING__MPI
  if (mpi->Size()>1) {
    m_basename+="_"+rpa->gen.Variable("RNG_SEED");
  }
#endif

switch (m_iotype)
    {
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

 p_xs= std::make_shared<HepMC::GenCrossSection>();
 m_run_info= std::make_shared<HepMC::GenRunInfo>();
 HepMC::GenRunInfo::ToolInfo tool;
 tool.name = std::string("SHERPA-MC");
 tool.version = std::string(SHERPA_VERSION)+"."+std::string(SHERPA_SUBVERSION);
 tool.description = std::string(SHERPA_NAME);
 m_run_info->tools().push_back(tool);
}

Output_HepMC3_Genevent::~Output_HepMC3_Genevent()
{
  p_writer->close();

}

void Output_HepMC3_Genevent::SetXS(const double& xs, const double& xserr)
{
  p_xs->set_cross_section(xs, xserr);
}

void Output_HepMC3_Genevent::Output(Blob_List* blobs, const double weight) 
{
  m_hepmc3.Sherpa2HepMC(blobs, m_run_info);
  HepMC::GenEvent* q=m_hepmc3.GenEvent();
  if (q) 
  {
  q->set_cross_section(p_xs);
  if (p_writer)    p_writer->write_event(*(q));
  }
}

void Output_HepMC3_Genevent::ChangeFile()
{
  /*This should be implemented in HepMC3 library.*/
}

DECLARE_GETTER(Output_HepMC3_Genevent,"HepMC3_GenEvent",
	       Output_Base,Output_Arguments);

Output_Base *ATOOLS::Getter<Output_Base,Output_Arguments,Output_HepMC3_Genevent>::
operator()(const Output_Arguments &args) const
{
  return new Output_HepMC3_Genevent(args);
}

void ATOOLS::Getter<Output_Base,Output_Arguments,Output_HepMC3_Genevent>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"HepMC3 GenEvent output";
}
