#include "Hadron_Decay_Handler.H"
#include "Hadrons.H"
#include "Message.H"
#include "Random.H"
#include "Vector.H"
#include "Data_Read.H"

using namespace SHERPA;
using namespace HADRONS;
using namespace ATOOLS;
using namespace std;


Hadron_Decay_Handler::Hadron_Decay_Handler(Hadrons * _hadrons) :
  m_decmodel(string("Hadrons")), m_mode(1),
  p_lund(NULL), p_hadrons(_hadrons)
{
  SwitchOfLundDecays();
}

Hadron_Decay_Handler::Hadron_Decay_Handler(Lund_Interface * _lund) :
  m_decmodel(string("Lund")), m_mode(0),
  p_lund(_lund), p_hadrons(NULL)
{
}

void Hadron_Decay_Handler::EraseTreated(std::set<int> * hadrons)
{
  if (m_mode==0) hadrons->clear();
  if (m_mode==1) {
    map<kf::code,Decay_Table *> * cans = p_hadrons->GetDecayMap();
    for (map<kf::code,Decay_Table *>::iterator citer=cans->begin();citer!=cans->end();citer++) {
      msg.Debugging()<<"Killing flavours: "<<citer->first<<" ("<<cans->size()<<" ) "<<hadrons->size()<<endl;
      hadrons->erase(int(citer->first));
      msg.Debugging()<<"                  "<<citer->first<<" ("<<cans->size()<<" ) "<<hadrons->size()<<endl;
    }
  }
}


Hadron_Decay_Handler::~Hadron_Decay_Handler() 
{
}

void Hadron_Decay_Handler::DeletePointers()
{
}

bool Hadron_Decay_Handler::FillHadronDecayBlobs(Particle *part,
						Blob_List *blob_list,
						Particle_List *part_list )
{
  msg_Tracking()<<"Hadron_Decay_Handler::FillHadronDecayBlobs "<<part->Flav()<<endl;
  msg.Debugging()<<"Momentum: "<<part->Momentum()<<endl;

  // perform decay 
  switch( m_mode ) {
  case 1: p_hadrons->PerformDecay( part, blob_list, part_list );
    break;
  case 0: p_lund->PerformDecay( part, blob_list, part_list );
    break;
  }
  
  return 1;
}

void Hadron_Decay_Handler::SwitchOfLundDecays()
{
  std::map<ATOOLS::kf::code,ATOOLS::Decay_Table *>::iterator dtiter;
  for (dtiter=p_hadrons->GetDecayMap()->begin();
       dtiter!=p_hadrons->GetDecayMap()->end();dtiter++) {
    p_lund->SwitchOfDecays(dtiter->first);
  }
}
