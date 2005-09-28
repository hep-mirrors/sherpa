#include "Hadron_Decay_Handler.H"
#include "Message.H"
#include "Random.H"
#include "Vector.H"
#include "Data_Read.H"


using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

#ifdef USING__Hadrons
using namespace HADRONS;
#endif


#ifdef USING__Hadrons
Hadron_Decay_Handler::Hadron_Decay_Handler(Hadrons * _hadrons) :
  m_decmodel(string("Hadrons")), m_mode(1),
  p_hadrons(_hadrons),
  p_lund(NULL)
{
  //std::cout<<"Hadron_Decay_Handler::Hadron_Decay_Handler(Lund_Interface * _hadrons) :"<<std::endl
  //	   <<"   "<<_hadrons<<" -> "<<p_hadrons<<" "<<std::endl;
  SwitchOfLundDecays();
}
#endif


Hadron_Decay_Handler::Hadron_Decay_Handler(Lund_Interface * _lund) :
  m_decmodel(string("Lund")), m_mode(0),
#ifdef USING__Hadrons
  p_hadrons(NULL), 
#endif
  p_lund(_lund)
{
  //std::cout<<"Hadron_Decay_Handler::Hadron_Decay_Handler(Lund_Interface * _lund) :"<<std::endl
  //	   <<"   "<<_lund<<" -> "<<p_lund<<" "<<std::endl;
}

void Hadron_Decay_Handler::EraseTreated(std::set<int> * hadrons)
{
  if (m_mode==0) hadrons->clear();
#ifdef USING__Hadrons
  if (m_mode==1) {
    map<kf::code,Decay_Table *> * cans = p_hadrons->GetDecayMap();
    for (map<kf::code,Decay_Table *>::iterator citer=cans->begin();citer!=cans->end();citer++) {
      //msg.Debugging()<<"Killing flavours: "<<citer->first<<" ("<<cans->size()<<" ) "<<hadrons->size()<<endl;
      hadrons->erase(int(citer->first));
      //msg.Debugging()<<"                  "<<citer->first<<" ("<<cans->size()<<" ) "<<hadrons->size()<<endl;
    }
  }
#endif
}


Hadron_Decay_Handler::~Hadron_Decay_Handler() 
{
}

void Hadron_Decay_Handler::DeletePointers()
{
}

void Hadron_Decay_Handler::PrepareDecays(Blob * blob) {
  //msg_Tracking()<<"Hadron_Decay_Handler::PrepareDecays "<<endl<<(*blob)<<endl;
  switch( m_mode ) {
#ifdef USING__Hadrons
  case 1: 
    break;
#endif
  case 0: 
    p_lund->PerformAllDecays(blob);
    break;
  }
}

bool Hadron_Decay_Handler::FillHadronDecayBlobs(Particle *part,
						Blob_List *blob_list,
						Particle_List *part_list )
{
  msg_Tracking()<<"Hadron_Decay_Handler::FillHadronDecayBlobs "<<part->Flav()<<endl
		<<"Momentum: "<<part->Momentum()<<endl;

  // perform decay 
  switch( m_mode ) {
#ifdef USING__Hadrons
  case 1: p_hadrons->PerformDecay( part, blob_list, part_list );
    break;
#endif
  case 0: 
    p_lund->PerformDecay( part, blob_list, part_list );
    break;
  }
  
  return 1;
}

void Hadron_Decay_Handler::SwitchOfLundDecays()
{
#ifdef USING__Hadrons
  std::map<ATOOLS::kf::code,ATOOLS::Decay_Table *>::iterator dtiter;
  for (dtiter=p_hadrons->GetDecayMap()->begin();
       dtiter!=p_hadrons->GetDecayMap()->end();dtiter++) {
    p_lund->SwitchOfDecays(dtiter->first);
  }
#endif
}
