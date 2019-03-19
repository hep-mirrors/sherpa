#include "AddOns/HepMC/Output_HepMC3_Genevent.H"
#include "HepMC/GenEvent.h"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"


#include "HepMC/Writer.h"
#include "HepMC/WriterAscii.h"
#include "HepMC/WriterHEPEVT.h"
#ifdef USING__HEPMC3__WRITERROOTTREE
#include "HepMC/WriterRootTree.h"
#endif
#ifdef USING__HEPMC3__WRITERROOT
#include "HepMC/WriterRoot.h"
#endif

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Output_HepMC3_Genevent::Output_HepMC3_Genevent(const Output_Arguments &args) :
  Output_Base("HepMC3")
{
  m_basename=args.m_outpath+"/"+args.m_outfile;
  m_ext=".hepmc3g";
  m_iotype=args.p_reader->GetValue<int>("HEPMC3_IO_TYPE",0);
  
#ifdef USING__MPI
  if (mpi->Size()>1) {
    m_basename+="_"+rpa->gen.Variable("RNG_SEED");
  }
#endif

   switch (m_iotype)
                {
                case 0:
                    p_writer=new HepMC::WriterAscii( m_basename);
                    break;
                case 1:
                    p_writer=new HepMC::WriterHEPEVT(m_basename);
                    break;
                case 3:
#ifdef USING__HEPMC3__WRITERROOT
                    p_writer=new HepMC::WriterRoot(m_basename);
#else 
                    THROW(fatal_error,"Asked for Root output, but Sherpa was compiled without Root output support.");                    
#endif
                    break;
                case 4:
#ifdef USING__HEPMC3__WRITERROOTTREE
                    p_writer=new HepMC::WriterRootTree(m_basename);
#else 
                    THROW(fatal_error,"Asked for RootTree output, but Sherpa was compiled without RootTree output support.");                    
#endif
                    break;
                default:
                    THROW(fatal_error, "Output format HEPMC3_IO_TYPE= is undefined.");
                    break;
                }
}

Output_HepMC3_Genevent::~Output_HepMC3_Genevent()
{
  p_writer->close();
}



void Output_HepMC3_Genevent::Output(Blob_List* blobs, const double weight) 
{
  m_hepmc3.Sherpa2HepMC(blobs, weight);
  HepMC::GenEvent* q=m_hepmc3.GenEvent();
  if (q) 
  {
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
