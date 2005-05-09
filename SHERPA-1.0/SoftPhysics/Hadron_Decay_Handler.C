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

bool Hadron_Decay_Handler::FillHadronDecayBlobs(Particle *part,ATOOLS::Blob_List *blob_list,
												Particle_List *part_list )
{
  msg_Tracking()<<"Hadron_Decay_Handler::FillHadronDecayBlobs "<<part->Flav()<<endl;
  msg.Debugging()<<"Momentum: "<<part->Momentum()<<endl;

  // perform decay with resting decayer (CMS system)
  bool anti = part->Flav().IsAnti();					// antiparticle ?
  FlavourSet daughters; 
  vector <Vec4D> v_mom;
  vector <Vec4D> v_pos;
  string type("");
  bool transform( false );								// Lorentztrafo momenta
  int status( 0 );
  switch( m_mode ) {
	case 1: daughters = p_hadrons->PerformDecay( v_mom, v_pos );
			type = "Sherpa"; 
			transform =  true;
			status = 1;									// no need for second check
			break;
	case 0: daughters = p_lund->PerformDecay( part, v_mom, v_pos );
			type = "Pythia_v6.214";
			transform =  false;
			status = 0;									// see if sherpa can treat remaining part
			break;
  }

  // transform momentum into Lab System and create blob
  Poincare lambda(part->Momentum());
  lambda.Invert();
  Blob * blob;											// decay blob
  blob = new Blob();
  blob->SetStatus(status);
  blob->SetType( btp::Hadron_Decay );
  blob->SetTypeSpec(type);
  blob->SetId();
  msg.Debugging()<<"created new blob: #"<<blob->Id()<<" with status "<<status<<endl;
  blob->AddToInParticles( part );
  if( part->Info() == 'P' ) part->SetInfo('p');
  if( part->Info() == 'D' ) part->SetInfo('d');
  blob_list->push_back( blob );
  Particle * particle;									// daughter part.
  Vec4D momentum;										// daughter mom.
  Flavour flav;											// daughter flav.
  FlSetConstIter dit;									// Iterator
  int i(0);

  // treat every daughter
  for( dit = daughters.begin(); dit != daughters.end(); dit++ ) {
	i++;
	v_mom[i] = ( transform )? lambda*v_mom[i] : v_mom[i];		// Lorentz trafo (if appl.)
	momentum = v_mom[i];
	flav = (*dit);
	if( anti ) flav = flav.Bar();
	msg.Debugging()<<"Daughters: "<<i<<"   "<<flav<<"  "<<momentum<<endl;
	particle = new Particle( -1, flav, momentum );
	if( part_list ) particle->SetNumber( part_list->size() );
	else particle->SetNumber( 0 );
	particle->SetStatus(1);
	particle->SetInfo('D');
	blob->SetPosition( v_pos[i] );
	if( part_list ) part_list->push_back( particle ); 
	blob->AddToOutParticles( particle );
	// check if daughter can be treated as well
	bool can_sherpa(false);
	switch( m_mode ) {
	  case 1 : can_sherpa = p_hadrons->FindDecay(particle->RefFlav()); 
			   break;
	  case 0:  can_sherpa = p_lund->FindDecay( particle );		   
			   break;
	}
	if( can_sherpa ) {
	  blob->SetStatus(1);
	  FillHadronDecayBlobs( particle, blob_list, part_list );
	}
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
