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
  std::cout<<"New Decay_Table_Reader("<<path+file<<")"<<std::endl;
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
  int              nchannels = 0;
  vector<int>      helpkfc;
  double           BR;
  Flavour          flav;
  Decay_Channel  * dc;

  for (int i=0;i<m_helpsvv.size();i++) {
    if (ExtractFlavours(helpkfc,m_helpsvv[i][0])) {
      BR = double(atof(m_helpsvv[i][1].c_str()));
      if (BR>1.e-6) {
	dc = new Decay_Channel(dt->Flav());
	for (int j=0;j<helpkfc.size();++j) {
	  flav = Flavour(kf::code(abs(helpkfc[j])));
	  if (helpkfc[j]<0) flav = flav.Bar();
	  dc->AddDecayProduct(flav);
	}
	dc->SetWidth(BR*dt->Flav().Width());
	dc->SetProcessName();
	if (m_helpsvv[i].size()>2) dc->SetMEType(m_helpsvv[i][2]);
	if (m_helpsvv[i].size()>3) dc->SetPSFile(m_helpsvv[i][3]);
	dt->AddDecayChannel(dc);
	nchannels++;
      }
    }
  }
  return nchannels;
}

void Decay_Table_Reader::FillInMatrixElementsAndPS(Decay_Table * dt)
{
  Decay_Channel        * dc;
  Hadron_Decay_Channel * hdc;
  for (int i=0;i<dt->NumberOfDecayChannels();i++) {
    dc  = dt->GetDecayChannel(i);
    hdc = new Hadron_Decay_Channel(dc);
    hdc->InitialisePhaseSpace(m_helpsvv[i]);
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
