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

FlavourSet Hadrons::PerformDecay( vector<Vec4D> &_mom, vector<Vec4D> &_pos )
{
  const int nchan = p_table->NumberOfDecayChannels();
  Decay_Channel * dec_channel;
  msg_Tracking()<<"Hadrons::PerformDecay(): Do Decay for "<<p_table->Flav()<<"."<<endl;
  double TotalWidth = p_table->TotalWidth();
  bool channel_choosen(0);
  int k(0);

  // dice decay channel acc. to BR
  while( !channel_choosen ) {
	k = rand() % nchan;									// dice decay channel
	double r = 1.*rand()/RAND_MAX;						// take it acc. to BR
	if( r < p_table->Width(k) / TotalWidth ) {
	  dec_channel = p_table->GetDecayChannel(k);
	  channel_choosen = 1;
	}
  }	// now a channel is choosen
//  k = 4;
//  dec_channel = p_table->GetDecayChannel(k);
  msg_Tracking()<<"     choosen channel: "<<k<<endl;
  if( p_channelmap->find(dec_channel ) == p_channelmap->end() ) {
	msg.Error()<<"Error in Hadrons::PerformDecay() \n"
	  		   <<"      Couldn't find appropriate channel pointer for "<<dec_channel->GetDecaying()
			   <<" decay. \n"
			   <<"      Don't know what to do, will abort."<<endl;
	abort();
  }
  Hadron_Decay_Channel * hdc = (*p_channelmap)[dec_channel];

  // choose a kinematics that corresponds to the ME kinematic distribution
  double value(0.);
  const double max = hdc->GetPS()->Maximum();
  int trials(0);											// number of trials
  do {
	value = hdc->Differential();							// current val. of |M|^2
	trials++;
  } while( ran.Get() > value/max );
//  ofstream f("trials.out",ios::app|ios::out );
//  f<<trials<<endl;
  msg.Info()<<"Hadrons::PerformDecay:  # Trials "<<trials
	        <<"   <=>  "<<100./trials<<" %"<<endl;
  Vec4D * mom_pointer = hdc->Momenta();									// diced momemta
  const int n = dec_channel->NumberOfDecayProducts()+1;
  _pos.clear(); _mom.clear();
  for( int i=0; i<n; i++ ) {
	_pos.push_back( Vec4D( 1., 1., 0., 0. ) );
	_mom.push_back( mom_pointer[i] );
  }

  return dec_channel->GetDecayProducts();
}
  
void Hadrons::ReadInDecayTables()
{
  Data_Reader reader = Data_Reader(string("->"),string(";"),string("!"));
  reader.AddComment("#");
  reader.AddComment("//");
  reader.SetInputPath(m_path);
  reader.SetInputFile(m_file);

  vector<vector<string> > Decayers;
  reader.SetMatrixType(reader.MTransposed);
  reader.MatrixFromFile(Decayers);

  p_decaymap = new map<ATOOLS::kf::code,ATOOLS::Decay_Table *>;
  p_channelmap = new map< Decay_Channel*, Hadron_Decay_Channel* >;
  Decay_Table * dt;
  Flavour fl;
  for (size_t i=0;i<Decayers.size();++i) {
	fl = Flavour(kf::code(atoi((Decayers[i][0]).c_str())));
	if (p_decaymap->find(fl.Kfcode())!=p_decaymap->end()) {
	  msg.Error()<<"ERROR in Hadrons::ReadInDecayTables() :"<<endl
				 <<"   Flavour "<<fl<<" already in map. Don't know what to do, will abort."<<endl;
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
