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

void Output_HepMC2_Genevent::Output(Blob_List* blobs, const double weight) 
{
  m_hepmc2.Sherpa2HepMC(blobs);
#ifdef USING__HEPMC2__IOGENEVENT
  HepMC::IO_GenEvent io(m_outstream);
  io.write_event(m_hepmc2.GenEvent());
#else
  THROW(fatal_error,"HepMC::IO_GenEvent asked for, but HepMC version too old.");
#endif
}
