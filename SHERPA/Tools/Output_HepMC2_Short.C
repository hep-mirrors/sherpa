#include "SHERPA/Tools/Output_HepMC2_Short.H"
#include "HepMC/GenEvent.h"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
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

Output_HepMC2_Short::Output_HepMC2_Short(std::string basename,
                                         std::string ext,int precision) :
  Output_Base(basename, ext, precision)
{
#ifdef USING__HEPMC2__IOGENEVENT
  p_iogenevent = new HepMC::IO_GenEvent(m_outstream);
  p_iogenevent->precision(precision);
#else
  THROW(fatal_error,"HepMC::IO_GenEvent asked for, but HepMC version too old.");
#endif
#ifdef HEPMC_HAS_CROSS_SECTION
  p_xs=new HepMC::GenCrossSection();
#endif
  p_event=new HepMC::GenEvent();
}

Output_HepMC2_Short::~Output_HepMC2_Short()
{
#ifdef USING__HEPMC2__IOGENEVENT
  delete p_iogenevent;
#endif
#ifdef HEPMC_HAS_CROSS_SECTION
  delete p_xs;
#endif
  delete p_event;
}

void Output_HepMC2_Short::SetXS(const double& xs, const double& xserr)
{
#ifdef HEPMC_HAS_CROSS_SECTION
  p_xs->set_cross_section(xs, xserr);
#endif
}

void Output_HepMC2_Short::Output(Blob_List* blobs, const double weight)
{
#ifdef USING__HEPMC2__IOGENEVENT
  p_event->clear();
  m_hepmc2.Sherpa2ShortHepMC(blobs, *p_event, weight);
#ifdef HEPMC_HAS_CROSS_SECTION
  p_event->set_cross_section(*p_xs);
#endif
  p_iogenevent->write_event(p_event);
#endif
}

void Output_HepMC2_Short::ChangeFile(string number)
{
#ifdef USING__HEPMC2__IOGENEVENT
  delete p_iogenevent;
  Output_Base::ChangeFile(number);
  p_iogenevent = new HepMC::IO_GenEvent(m_outstream);
#endif
}
