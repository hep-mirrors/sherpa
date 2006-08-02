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
    map<ATOOLS::Flavour,Decay_Table *> * cans = p_hadrons->GetDecayMap();
    for (map<ATOOLS::Flavour,Decay_Table *>::iterator citer=cans->begin();citer!=cans->end();citer++) {
      hadrons->erase(int(citer->first.Kfcode()));
      Spin_Correlation_Tensor::AddPossibleParticle( citer->first.Kfcode() );
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

// with spin correlation 
bool Hadron_Decay_Handler::FillHadronDecayBlobs(Particle *part,
						Blob_List *blob_list,
                        Spin_Density_Matrix * sigma, 
                        Spin_Density_Matrix * decmatr, 
						Particle_List *part_list )
{
  msg_Tracking()<<"Hadron_Decay_Handler::FillHadronDecayBlobs "<<part->Flav()<<endl
    <<"     with spin correlations."<<endl
    <<"     Momentum: "<<part->Momentum()<<endl;

  // perform decay 
  switch( m_mode ) {
#ifdef USING__Hadrons
    case 1: {
      // create first decay blobs to start the recursion with a head start for mass smearing
      Blob* blob=p_hadrons->CreateDecayBlobSkeleton(part,blob_list,part_list);
      
      if(blob) {
        if( decmatr ) (*decmatr) = p_hadrons->PerformDecay( blob, blob_list, part_list, sigma );
        else p_hadrons->PerformDecay( blob, blob_list, part_list, NULL );
      }
      break;
    }
#endif
    case 0: {
      p_lund->PerformDecay( part, blob_list, part_list );
      break;
    }
  }
  
  return 1;
}

void Hadron_Decay_Handler::SwitchOfLundDecays()
{
#ifdef USING__Hadrons
  std::map<ATOOLS::Flavour,ATOOLS::Decay_Table *>::iterator dtiter;
  for (dtiter=p_hadrons->GetDecayMap()->begin();
       dtiter!=p_hadrons->GetDecayMap()->end();dtiter++) {
    p_lund->SwitchOfDecays(dtiter->first.Kfcode());
  }
#endif
}
