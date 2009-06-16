#include "SHERPA/Tools/Output_HepEvt.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

void Output_HepEvt::Output(Blob_List* blobs, const double weight) 
{
  m_hepevt.Sherpa2HepEvt(blobs);
  m_hepevt.SetWeight(weight);
  m_hepevt.WriteFullHepEvt(m_outstream,m_hepevt.Nhep());
}
