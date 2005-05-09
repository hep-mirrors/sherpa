#include "Decay_Table_Reader.H"
#include "Hadron_Decay_Channel.H"
#include "Decay_Table.H"
#include "Data_Reader.H"
#include "Message.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

Decay_Table_Reader::Decay_Table_Reader(string path,string file)
{
  if (file==string("")) {
    msg.Error()<<"Error in Decay_Table_Reader::Decay_Table_Reader("
	       <<path<<","<<file<<") : "<<endl
	       <<"   No file specified, will return and hope for the best."<<endl;
    return;
  }
  msg_Info()<<"New Decay_Table_Reader("<<path+file<<")"<<endl;
  Data_Reader reader = Data_Reader(string("|"),string(";"),string("!"));
  reader.AddComment("#");
  reader.AddComment("//");
  reader.SetInputPath(path);
  reader.SetInputFile(file);
  reader.SetMatrixType(reader.MTransposed);
  reader.MatrixFromFile(m_helpsvv,"");
}

int Decay_Table_Reader::FillDecayTable(Decay_Table * dt)
{
  int              nchannels = 0;									// number of channels
  vector<int>      helpkfc;											// aux. kf-code used to extract flavours
  double           BR;												// branching ratio
  Flavour          flav;											// involved flavors
  Decay_Channel  * dc;												// pointer on decay channel

  for (int i=0;i<m_helpsvv.size();i++) {							// for each element
    if (ExtractFlavours(helpkfc,m_helpsvv[i][0])) {					//   extract flavours
      BR = double(atof(m_helpsvv[i][1].c_str()));					//   get branching ratio
      if (BR>1.e-6) {												
		dc = new Decay_Channel(dt->Flav());							//   create new decay channel
		for (int j=0;j<helpkfc.size();++j) {
		  flav = Flavour(kf::code(abs(helpkfc[j])));
		  if (helpkfc[j]<0) flav = flav.Bar();
		  dc->AddDecayProduct(flav);				
		}
		dc->SetWidth(BR*dt->Flav().Width());						//   set decay width
		dc->SetProcessName();										//   set process name
		if (m_helpsvv[i].size()>2) dc->SetMEType(m_helpsvv[i][2]);
		if (m_helpsvv[i].size()>3) dc->SetPSFile(m_helpsvv[i][3]);
		dt->AddDecayChannel(dc);									//   add to ex. decay channels
		nchannels++;								
      } // if BR>...
    }
  } // for
  return nchannels;												// return number of channels
}

void Decay_Table_Reader::FillInMatrixElementsAndPS( Decay_Table * dt, Channel_Map * chmap )
{
  Decay_Channel        * dc;								// pointer on decay channel
  Hadron_Decay_Channel * hdc;								// pointer on hadronic decay channel
  for (int i=0;i<dt->NumberOfDecayChannels();i++) {			// for each channel in decay table
    dc  = dt->GetDecayChannel(i);							//   get channel
    hdc = new Hadron_Decay_Channel(dc);						//   create new hadronic channel
    hdc->InitialisePhaseSpace(m_helpsvv[i]);				//   initialise phase space
	(*chmap)[ dc ] = hdc;									//   map channel to hadronic channel
  }
}

bool Decay_Table_Reader::ExtractFlavours(vector<int> & helpkfc,string help)
{
  helpkfc.clear();													
  size_t pos = help.find("{");
  bool             hit;
  if (pos!=string::npos) help = help.substr(pos+1);
  else {
    msg.Error()<<"WARNING in Decay_Table_Reader:: : "<<endl
	       <<"   Something wrong with final state of decay (Bracket missing) :"<<help<<endl
	       <<"   Will skip it."<<endl;
    return false;
  }
  pos    = help.find("}");
  if (pos!=string::npos) help = help.substr(0,pos);
  else {
    msg.Error()<<"WARNING in Decay_Table_Reader:: : "<<endl
	       <<"   Something wrong with final state of decay (Bracket missing) :"<<help<<endl
	       <<"   Will skip it."<<endl;
    return false;
  }
  hit    = true;
  while (hit) {
    pos      = help.find(",");
    if (pos!=string::npos) {
      helpkfc.push_back(atoi((help.substr(0,pos)).c_str()));
      help  = help.substr(pos+1);
    }
    else {
      helpkfc.push_back(atoi(help.c_str()));
      hit = false;
    }
  }
  if (helpkfc.size()<2) {
    msg.Error()<<"WARNING in Decay_Table_Reader:: : "<<endl
	       <<"   Something wrong with final state of decay (Too little particles) :";
    for (int j=0;j<helpkfc.size();j++) msg.Error()<<helpkfc[j]<<" ";
    msg.Error()<<endl<<"   Will skip it."<<endl;
    return false;
  } 
  return true;
}
