#include "Output_HepMC2_Genevent.H"
#include "HepMC/GenEvent.h"
#include "CXXFLAGS_PACKAGES.H"
#include "Exception.H"
#ifdef USING__HEPMC2__IOGENEVENT
#include "HepMC/IO_GenEvent.h"
#endif

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Output_HepMC2_Genevent::Output_HepMC2_Genevent(std::string basename,
                                               std::string ext,int precision) :
  Output_Base(basename, ext, precision)
{
#ifdef USING__HEPMC2__IOGENEVENT
  p_iogenevent = new HepMC::IO_GenEvent(m_outstream);
#else
  THROW(fatal_error,"HepMC::IO_GenEvent asked for, but HepMC version too old.");
#endif
}

Output_HepMC2_Genevent::~Output_HepMC2_Genevent()
{
#ifdef USING__HEPMC2__IOGENEVENT
  delete p_iogenevent;
#endif
}

void Output_HepMC2_Genevent::Output(Blob_List* blobs, const double weight) 
{
#ifdef USING__HEPMC2__IOGENEVENT
  m_hepmc2.Sherpa2HepMC(blobs, weight);
  p_iogenevent->write_event(m_hepmc2.GenEvent());
#endif
}

void Output_HepMC2_Genevent::ChangeFile(string number)
{
#ifdef USING__HEPMC2__IOGENEVENT
  delete p_iogenevent;
  Output_Base::ChangeFile(number);
  p_iogenevent = new HepMC::IO_GenEvent(m_outstream);
#endif
}
