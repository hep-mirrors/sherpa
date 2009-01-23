#include "Output_HepMC2.H"
#include "HepMC/GenEvent.h"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

void Output_HepMC2::Output(Blob_List* blobs, const double weight) 
{
  m_hepmc2.Sherpa2HepMC(blobs, weight);
  m_hepmc2.GenEvent()->print(m_outstream);
}
