#include "Hadrons.H"
#include "Data_Reader.H"
#include "Flavour.H"
#include "Message.H"
#include "Random.H"
#include <time.h>
#include "Hadron_Decay_Channel.H"
#include <fstream>

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

Hadrons::Hadrons(string _path,string _file) : 
  m_path(_path), m_file(_file) 
{ 
  msg_Tracking()<<"In Hadrons: ("<<_path<<") "<<_file<<std::endl;
  ReadInDecayTables();
}

bool Hadrons::FindDecay(const ATOOLS::Flavour & Decayer)
{
  if (p_decaymap->find(Decayer.Kfcode())==p_decaymap->end()) return false;
  p_table = (*p_decaymap)[Decayer.Kfcode()];
  return true;
}

Hadron_Decay_Channel * Hadrons::ChooseDecayChannel()
{
  const int nchan = p_table->NumberOfDecayChannels();
  Decay_Channel * dec_channel;
  double TotalWidth = p_table->TotalWidth();
  bool channel_choosen (0);
  int k (0);

  // dice decay channel acc. to BR
  while( !channel_choosen ) {
	k = rand() % nchan;									// dice decay channel
	double r = 1.*rand()/RAND_MAX;						// take it acc. to BR
	if( r < p_table->Width(k) / TotalWidth ) {
	  dec_channel = p_table->GetDecayChannel(k);
	  channel_choosen = 1;
	}
  }	

//  k = 4;
//  dec_channel = p_table->GetDecayChannel(k);
  msg_Tracking()<<"     choosen channel: "<<k<<endl;
  if( p_channelmap->find(dec_channel ) == p_channelmap->end() ) {
	msg.Error()<<"Error in Hadrons::ChooseDecayChannel() \n"
	  		   <<"      Couldn't find appropriate channel pointer for "<<dec_channel->GetDecaying()
			   <<" decay. \n"
			   <<"      Don't know what to do, will abort."<<endl;
	abort();
  }
  return (*p_channelmap)[dec_channel];
}

void Hadrons::ChooseDecayKinematics( Vec4D * _p, Vec4D * _v, Hadron_Decay_Channel * _hdc )
{
  double value(0.);
  const double max = _hdc->GetPS()->Maximum();
  int trials(0);											// number of trials
  do {
	value = _hdc->Differential();							// current val. of |M|^2
	trials++;
  } while( ran.Get() > value/max );
//  ofstream f("trials.out",ios::app|ios::out );
//  f<<trials<<endl;
  msg_Tracking()<<"Hadrons::ChooseDecayKinematics:  # Trials "<<trials
  	        <<"   <=>  "<<100./trials<<" %"<<endl;
  for( int i=0; i < _hdc->DecayChannel()->NumberOfDecayProducts()+1; i++ ) {
	_p[i] = _hdc->Momenta(i);
	_v[i] = Vec4D( 1., 1., 0., 0. );
  }
}

void Hadrons::PerformDecay( Particle * part, Blob_List * blob_list, Particle_List * part_list )
{
  msg_Tracking()<<"Hadrons::PerformDecay() for "<<part->Flav()<<endl;
  msg.Debugging()<<"Momentum: "<<part->Momentum()<<endl;
  if( part->Flav().Kfcode() != p_table->Flav().Kfcode() ) {
	msg.Error()<<"ERROR in Hadrons::PerformDecay() : Particle in selected decay table is not \n"
	           <<"        identical to the decayer.\n"
			   <<"        >>>   "<<p_table->Flav()<<" <-> "<<part->Flav()<<"   <<<"<<endl
			   <<"        Don't know which to take, will abort."<<endl;
	abort();
  }

  bool anti = part->Flav().IsAnti();					// antiparticle ?

  // choose decay channel acc. to BR
  Hadron_Decay_Channel * hdc = ChooseDecayChannel();
  FlavourSet daughters = hdc->DecayChannel()->GetDecayProducts();
  const int n = hdc->DecayChannel()->NumberOfDecayProducts()+1;

  // choose a kinematics that corresponds to the ME kinematic distribution
  Vec4D mom[n];
  Vec4D pos[n];
  ChooseDecayKinematics( mom, pos, hdc ); 

  // transform momentum into Lab System and create blob
  Poincare lambda(part->Momentum());
  lambda.Invert();
  Blob * blob;											// decay blob
  blob = new Blob();
  blob->SetStatus(1);
  blob->SetType( btp::Hadron_Decay );
  blob->SetTypeSpec( "Sherpa" );
  blob->SetId();
  msg.Debugging()<<"created new blob: #"<<blob->Id()<<" with status 1"<<endl;
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
	momentum = lambda * mom[i];							// Lorentz trafo
	flav = (*dit);
	if( anti ) flav = flav.Bar();
	msg.Debugging()<<"Daughters: "<<i<<"   "<<flav<<"  "<<momentum<<endl;
	particle = new Particle( -1, flav, momentum );
	if( part_list ) particle->SetNumber( part_list->size() );
	else particle->SetNumber( 0 );
	particle->SetStatus(1);
	particle->SetInfo('D');
	blob->SetPosition( pos[i] );
	if( part_list ) part_list->push_back( particle ); 
	blob->AddToOutParticles( particle );
	// check if daughter can be treated as well
	if( FindDecay(particle->RefFlav()) ) {
	  PerformDecay( particle, blob_list, part_list );
	}
  }
}
  
void Hadrons::ReadInDecayTables()
{
  Data_Reader reader = Data_Reader(string("->"),string(";"),string("!"));
  reader.SetAddCommandLine(false);
  reader.AddComment("#");
  reader.AddComment("//");
  reader.SetInputPath(m_path);
  reader.SetInputFile(m_file);

  vector<vector<string> > Decayers;
  reader.SetMatrixType(reader.MTransposed);
  if(!reader.MatrixFromFile(Decayers)) {
    msg.Error()<<"ERROR in Hadrons::ReadInDecayTables() :\n"
	       <<"   Read in failure, will abort."<<endl;
    abort();
  }

  p_decaymap = new map<ATOOLS::kf::code,ATOOLS::Decay_Table *>;
  p_channelmap = new map< Decay_Channel*, Hadron_Decay_Channel* >;
  Decay_Table * dt;
  Flavour fl;
  for (size_t i=0;i<Decayers.size();++i) {
    fl = Flavour(kf::code(atoi((Decayers[i][0]).c_str())));
    if (p_decaymap->find(fl.Kfcode())!=p_decaymap->end()) {
      msg.Error()<<"ERROR in Hadrons::ReadInDecayTables() :"<<endl
		 <<"   Flavour "<<fl
		 <<" already in map. Don't know what to do, will abort."<<endl;
      abort();
    }
    dt = InitialiseOneDecayTable(Decayers[i]);
    (*p_decaymap)[fl.Kfcode()] = dt;
  }
}

Decay_Table * Hadrons::InitialiseOneDecayTable(vector<string> line)
{
  Decay_Table * dt              = new Decay_Table(Flavour(kf::code(atoi((line[0]).c_str()))));
  Decay_Table_Reader * dtreader = new Decay_Table_Reader(m_path,line[1]);
  if (dtreader->FillDecayTable(dt)>0) {
    msg.Out()<<"Found "<<dt->NumberOfDecayChannels()<<" decay channels for "<<dt->Flav()<<endl;
    dtreader->FillInMatrixElementsAndPS(dt,p_channelmap);
    msg.Out()<<"Initialised a new decay table : "<<endl;dt->Output();msg.Out()<<endl;
  }
  else { 
    msg.Error()<<"WARNING in Hadrons::InitialiseOneDecayTable : "<<endl
	       <<"   No decay channels found for "<<dt->Flav()<<" in file "<<line[1]<<endl
	       <<"   Will continue and hope for the best."<<endl;
    delete dt; 
    dt = NULL;
  }
  delete dtreader;
  return dt;
}
