#include "AddOns/HepMC/Output_HepMC2_Short.H"
#include "HepMC/GenEvent.h"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include "HepMC/IO_GenEvent.h"
#include "HepMC/HepMCDefs.h"
#include "HepMC/GenCrossSection.h"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Output_HepMC2_Short::Output_HepMC2_Short(const Output_Arguments &args) :
  Output_Base{ "HepMC2S" }
{
  m_basename=args.m_outpath+"/"+args.m_outfile;
  m_ext=".hepmc";
  const int precision{
    Settings::GetMainSettings()["EVENT_OUTPUT_PRECISION"].Get<int>() };
  p_iogenevent = new HepMC::IO_GenEvent(m_outstream);
  p_iogenevent->precision(precision);
  p_xs=new HepMC::GenCrossSection();
  p_event=new HepMC::GenEvent();
#ifdef USING__GZIP
  m_ext += ".gz";
#endif
#ifdef USING__MPI
  if (mpi->Size()>1) {
    m_basename+="_"+rpa->gen.Variable("RNG_SEED");
  }
#endif
  m_outstream.open((m_basename+m_ext).c_str());
  if (!m_outstream.good())
    THROW(fatal_error, "Could not open event file "+m_basename+m_ext+".");
  m_outstream.precision(precision);
}

Output_HepMC2_Short::~Output_HepMC2_Short()
{
  delete p_iogenevent;
  delete p_xs;
  m_outstream.close();
  delete p_event;
}

void Output_HepMC2_Short::SetXS(const ATOOLS::Weights_Map& xs,
			        const ATOOLS::Weights_Map& xserr)
{
  p_xs->set_cross_section(xs.Nominal(), xserr.Nominal());
}

void Output_HepMC2_Short::Output(Blob_List* blobs)
{
  p_event->clear();
  m_hepmc2.Sherpa2ShortHepMC(blobs, *p_event);
  std::vector<HepMC::GenEvent*> subevents(m_hepmc2.GenSubEventList());
  p_event->set_cross_section(*p_xs);
  for (size_t i(0); i<subevents.size(); ++i) {
    subevents[i]->set_cross_section(*p_xs);
  }
  if (subevents.size()) {
    for (size_t i(0); i<subevents.size(); ++i) {
      p_iogenevent->write_event(subevents[i]);
    }
    m_hepmc2.DeleteGenSubEventList();
  } else {
    p_iogenevent->write_event(p_event);
  }
}

void Output_HepMC2_Short::ChangeFile()
{
  delete p_iogenevent;
  m_outstream.close();
  std::string newname(m_basename+m_ext);
  for (size_t i(0);FileExists(newname);
       newname=m_basename+"."+ToString(++i)+m_ext);
  m_outstream.open(newname.c_str());
  if (!m_outstream.good())
    THROW(fatal_error, "Could not open event file "+newname+".")
  p_iogenevent = new HepMC::IO_GenEvent(m_outstream);
}

DECLARE_GETTER(Output_HepMC2_Short,"HepMC_Short",
	       Output_Base,Output_Arguments);

Output_Base *ATOOLS::Getter<Output_Base,Output_Arguments,Output_HepMC2_Short>::
operator()(const Output_Arguments &args) const
{
  return new Output_HepMC2_Short(args);
}

void ATOOLS::Getter<Output_Base,Output_Arguments,Output_HepMC2_Short>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"HepMC short output";
}

