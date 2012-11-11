#include "AddOns/HepMC/Output_HepMC2_Genevent.H"
#include "HepMC/GenEvent.h"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

#ifdef USING__HEPMC2__IOGENEVENT
#include "HepMC/IO_GenEvent.h"
#endif

#ifdef USING__HEPMC2__DEFS
#include "HepMC/HepMCDefs.h"
#ifdef HEPMC_HAS_CROSS_SECTION
#include "HepMC/GenCrossSection.h"
#endif
#endif

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Output_HepMC2_Genevent::Output_HepMC2_Genevent(const Output_Arguments &args) :
  Output_Base("HepMC2")
{
  m_basename=args.m_outpath+"/"+args.m_outfile;
  m_ext=".hepmc2g";
  int precision       = args.p_reader->GetValue<int>("OUTPUT_PRECISION",12);
#ifdef USING__HEPMC2__IOGENEVENT
  p_iogenevent = new HepMC::IO_GenEvent(m_outstream);
#ifdef HEPMC_HAS_CROSS_SECTION
  p_iogenevent->precision(precision);
#endif
#else
  THROW(fatal_error,"HepMC::IO_GenEvent asked for, but HepMC version too old.");
#endif
#ifdef HEPMC_HAS_CROSS_SECTION
  p_xs=new HepMC::GenCrossSection();
#endif
#ifdef USING__GZIP
  m_ext += ".gz";
#endif
  m_outstream.open((m_basename+m_ext).c_str());
  if (!m_outstream.good())
    THROW(fatal_error, "Could not open event file "+m_basename+m_ext+".");
  m_outstream.precision(precision);
}

Output_HepMC2_Genevent::~Output_HepMC2_Genevent()
{
  m_outstream.close();
#ifdef USING__HEPMC2__IOGENEVENT
  delete p_iogenevent;
#endif
#ifdef HEPMC_HAS_CROSS_SECTION
  delete p_xs;
#endif
}

void Output_HepMC2_Genevent::SetXS(const double& xs, const double& xserr)
{
#ifdef HEPMC_HAS_CROSS_SECTION
  p_xs->set_cross_section(xs, xserr);
#endif
}

void Output_HepMC2_Genevent::Output(Blob_List* blobs, const double weight) 
{
#ifdef USING__HEPMC2__IOGENEVENT
  m_hepmc2.Sherpa2HepMC(blobs, weight);
#ifdef HEPMC_HAS_CROSS_SECTION
  m_hepmc2.GenEvent()->set_cross_section(*p_xs);
#endif
  p_iogenevent->write_event(m_hepmc2.GenEvent());
#endif
}

void Output_HepMC2_Genevent::ChangeFile()
{
#ifdef USING__HEPMC2__IOGENEVENT
  delete p_iogenevent;
  m_outstream.close();
  std::string newname(m_basename+m_ext);
  for (size_t i(0);FileExists(newname);
       newname=m_basename+"."+ToString(++i)+m_ext);
  m_outstream.open(newname.c_str());
  if (!m_outstream.good())
    THROW(fatal_error, "Could not open event file "+newname+".")
  p_iogenevent = new HepMC::IO_GenEvent(m_outstream);
#endif
}

DECLARE_GETTER(HepMC_GenEvent_Output_Getter,"HepMC_GenEvent",
	       Output_Base,Output_Arguments);

Output_Base *HepMC_GenEvent_Output_Getter::operator()
(const Output_Arguments &args) const
{
  return new Output_HepMC2_Genevent(args);
}

void HepMC_GenEvent_Output_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"HepMC GenEvent output";
}

