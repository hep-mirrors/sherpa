#include "Hadron_Decays.H"
#include "Data_Reader.H"
#include "Flavour.H"
#include "Message.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

Hadron_Decays::Hadron_Decays(string _path,string _file) : 
  m_path(_path), m_file(_file) 
{ 
  msg.Out()<<"In Hadron_Decays: |"<<_path<<"|"<<_file<<std::endl;
  ReadInDecayTables();
  abort();
}

void Hadron_Decays::ReadInDecayTables()
{
  Data_Reader reader = Data_Reader(string("->"),string(";"),string("!"));
  reader.AddComment("#");
  reader.AddComment("//");
  reader.SetInputPath(m_path);
  reader.SetInputFile(m_file);

  vector<vector<string> > Decayers;
  reader.SetMatrixType(reader.MTransposed);
  reader.MatrixFromFile(Decayers);

  Decay_Table * dt;
  for (size_t i=0;i<Decayers.size();++i) {
    dt = InitialiseOneDecayTable(Decayers[i]);
  }
}

Decay_Table * Hadron_Decays::InitialiseOneDecayTable(vector<string> line)
{
  Decay_Table * dt              = new Decay_Table(Flavour(kf::code(atoi((line[0]).c_str()))));
  Decay_Table_Reader * dtreader = new Decay_Table_Reader(m_path,line[1]);
  if (dtreader->FillDecayTable(dt)>0) {
    cout<<"Found "<<dt->NumberOfDecayChannels()<<" decay channels for "<<dt->Flav()<<endl;
    dtreader->FillInMatrixElementsAndPS(dt);
    msg.Out()<<"Initialised a new decay table : "<<endl;dt->Output();msg.Out()<<endl;
  }
  else { 
    msg.Error()<<"WARNING in Hadron_Decays::InitialiseOneDecayTable : "<<endl
	       <<"   No decay channels found for "<<dt->Flav()<<" in file "<<line[1]<<endl
	       <<"   Will continue and hope for the best."<<endl;
    delete dt; 
    dt = NULL;
  }
  delete dtreader;
  return dt;
}
