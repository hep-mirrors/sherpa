#include "Selector.H"
#include "Message.H"

#include "Algebra_Interpreter.H"
#include "Data_Reader.H"
#include "Exception.H"
#include "Run_Parameter.H"
#include "MyStrStream.H"

using namespace ATOOLS;
using namespace std;


Selector_Base::~Selector_Base() 
{ 
  if (m_sel_log!=NULL) delete m_sel_log;
}

void Selector_Base::Output() { 
  if (!(msg_LevelIsTracking())) return;
  if(m_sel_log) {
    m_sel_log->Output();
    msg_Out()<<m_name<<"  total number of rejections: "<<m_sel_log->Rejections()<<std::endl;
  }
}

void Selector_Base::Add(Selector_Base *) {
  msg_Error()<<"Selector_Base::Add : Virtual method."<<std::endl;
}

void Selector_Base::BuildCuts(Cut_Data *) { 
  msg_Error()<<"Selector_Base::BuildCuts : Virtual method."<<std::endl;
}

void Selector_Base::UpdateCuts(double,double,Cut_Data *) { 
  msg_Error()<<"Selector_Base::BuildCuts : Virtual method."<<std::endl;
}

void Selector_Base::SetRange(std::vector<Flavour>,double,double) { 
  msg_Error()<<"Selector_Base::SetRange : Virtual method."<<std::endl;
}

void Selector_Base::SetRange(std::vector<Flavour>,int,double,double) { 
  msg_Error()<<"Selector_Base::SetRange : Virtual method."<<std::endl;
}

void Selector_Base::SetRange(std::vector<Flavour>,
			     std::vector<std::pair<double,double> > &)
{
  msg_Error()<<"Selector_Base::SetRange : Virtual method."<<std::endl;
}

bool Selector_Base::GetValue(const std::string &name,double &value)
{
  msg_Error()<<"Selector_Base::GetValue("<<name<<",..): "
	     <<"Virtual method called."<<std::endl;
  return false;
}

double Selector_Base::ActualValue() const
{ 
  return 2.; 
}

int    Selector_Base::NeedUpdate()                         { return 0; }
int    Selector_Base::IsConditional()                      { return 0; }
void   Selector_Base::SetSRange(double _smin,double _smax) { m_smin = _smin; m_smax = _smax; }
void   Selector_Base::SetName(std::string _name)           { m_name = _name; }
//double Selector_Base::Smin()                               { return m_smin; }
//double Selector_Base::Smax()                               { return m_smax; }
std::string Selector_Base::Name()                          { return m_name; }
void   Selector_Base::SetProcessName(const std::string &name) { m_procname=name; }
std::string Selector_Base::ProcessName() const { return m_procname; }

/*-----------------------------------------------------------------------------------

  Selector_Data

  -----------------------------------------------------------------------------------*/


Selector_Data::Selector_Data() {}

Selector_Data::Selector_Data(std::string path, std::string filename) {
  if (!ReadInData(path, filename)) {
    msg_Error()<<"Error in Selector_Data::Selector_Data("<<path<<")."<<endl
			  <<"Cannot initialise any selector. Abort."<<endl;
    abort();
  }
  ControlOutput();
}

bool Selector_Data::ReadInData(std::string path, std::string filename) 
{
  Data_Reader reader(" ","\\;","!");
  reader.AddWordSeparator("\t");
  reader.AddComment("#");
  reader.AddComment("//");
  reader.SetAddCommandLine(false);
  reader.SetInputPath(path);
  reader.SetInputFile(filename);
  reader.SetMatrixType(mtc::transposed);

  vector<vector<string> > svv;
  reader.MatrixFromFile(svv);

  Algebra_Interpreter * ip = reader.Interpreter();
  
  std::string keyword;
  Mom_Data    dat;
  int         crit1,crit2;
  Flavour     flav,flav2;
  for (size_t i=0;i<svv.size();++i) {
    dat.type=-1;
    dat.flavs.clear();
    dat.bounds.resize(1);
    keyword=svv[i][0];
    if (keyword == string("JetFinder")) {
      dat.type = 1;
      rpa.gen.SetVariable("Y_CUT",svv[i][1]);
      rpa.gen.SetVariable("DELTA_R",svv[i][2]);
      dat.bounds.front().first=2.0;
      dat.bounds.front().second=1.0;
    }
    else if (keyword == string("ConeFinder")) {
      dat.type = 2;
      dat.bounds.front().first=ToType<double>(ip->Interprete(svv[i][1]));
    }
    else if (keyword == string("DipoleFinder")) {
      dat.type = 3;
      dat.bounds.front().first=ToType<double>(ip->Interprete(svv[i][1]));
      dat.bounds.front().second=ToType<double>(ip->Interprete(svv[i][2]));
    }
    else if (keyword == string("Energy")) {
      dat.type = 11;
      crit1=ToType<int>(ip->Interprete(svv[i][1]));
      dat.bounds.front().first=ToType<double>(ip->Interprete(svv[i][2]));
      dat.bounds.front().second=ToType<double>(ip->Interprete(svv[i][3]));
      flav = Flavour((kf_code)abs(crit1));
      if (crit1<0) flav = flav.Bar();
      dat.flavs.push_back(flav);
    }
    else if (keyword == string("PT")) {
      dat.type = 12;
      crit1=ToType<int>(ip->Interprete(svv[i][1]));
      dat.bounds.front().first=ToType<double>(ip->Interprete(svv[i][2]));
      dat.bounds.front().second=ToType<double>(ip->Interprete(svv[i][3]));
      flav = Flavour((kf_code)abs(crit1));
      if (crit1<0) flav = flav.Bar();
      dat.flavs.push_back(flav);
    }
    else if (keyword == string("Rapidity")) {
      dat.type = 13;
      crit1=ToType<int>(ip->Interprete(svv[i][1]));
      dat.bounds.front().first=ToType<double>(ip->Interprete(svv[i][2]));
      dat.bounds.front().second=ToType<double>(ip->Interprete(svv[i][3]));
      flav = Flavour((kf_code)abs(crit1));
      if (crit1<0) flav = flav.Bar();
      dat.flavs.push_back(flav);
    }
    else if (keyword == string("BeamAngle")) {
      dat.type = 14;
      crit1=ToType<int>(ip->Interprete(svv[i][1]));
      crit2=ToType<int>(ip->Interprete(svv[i][2]));
      dat.bounds.front().first=ToType<double>(ip->Interprete(svv[i][3]));
      dat.bounds.front().second=ToType<double>(ip->Interprete(svv[i][4]));
      flav = Flavour((kf_code)abs(crit1));
      if (crit1<0) flav = flav.Bar();
      dat.flavs.push_back(flav);
      dat.help = crit2;
    }  
    else if (keyword == string("ET")) {
      dat.type = 15;
      crit1=ToType<int>(ip->Interprete(svv[i][1]));
      dat.bounds.front().first=ToType<double>(ip->Interprete(svv[i][2]));
      dat.bounds.front().second=ToType<double>(ip->Interprete(svv[i][3]));
      flav = Flavour((kf_code)abs(crit1));
      if (crit1<0) flav = flav.Bar();
      dat.flavs.push_back(flav);
    }
    else if (keyword == string("PseudoRapidity")) {
      dat.type = 16;
      crit1=ToType<int>(ip->Interprete(svv[i][1]));
      dat.bounds.front().first=ToType<double>(ip->Interprete(svv[i][2]));
      dat.bounds.front().second=ToType<double>(ip->Interprete(svv[i][3]));
      flav = Flavour((kf_code)abs(crit1));
      if (crit1<0) flav = flav.Bar();
      dat.flavs.push_back(flav);
    }
    else if (keyword == string("Mass")) {
      dat.type = 21;
      crit1=ToType<int>(ip->Interprete(svv[i][1]));
      crit2=ToType<int>(ip->Interprete(svv[i][2]));
      dat.bounds.front().first=ToType<double>(ip->Interprete(svv[i][3]));
      dat.bounds.front().second=ToType<double>(ip->Interprete(svv[i][4]));
      flav = Flavour((kf_code)abs(crit1));
      if (crit1<0) flav = flav.Bar();
      (dat.flavs).push_back(flav);
      flav = Flavour((kf_code)abs(crit2));
      if (crit2<0) flav = flav.Bar();
      dat.flavs.push_back(flav);
    }
    else if (keyword == string("Angle")) {
      dat.type = 22;
      crit1=ToType<int>(ip->Interprete(svv[i][1]));
      crit2=ToType<int>(ip->Interprete(svv[i][2]));
      dat.bounds.front().first=ToType<double>(ip->Interprete(svv[i][3]));
      dat.bounds.front().second=ToType<double>(ip->Interprete(svv[i][4]));
      flav = Flavour((kf_code)abs(crit1));
      if (crit1<0) flav = flav.Bar();
      (dat.flavs).push_back(flav);
      flav = Flavour((kf_code)abs(crit2));
      if (crit2<0) flav = flav.Bar();
      dat.flavs.push_back(flav);
    }  
    else if (keyword == string("X")) {
      dat.type = 24;
      dat.bounds.front().first=ToType<double>(ip->Interprete(svv[i][1]));
      dat.bounds.front().second=ToType<double>(ip->Interprete(svv[i][2]));
    }
    else if (keyword == string("BFKL_PT")) {
      dat.type = 25;
      crit1=ToType<int>(ip->Interprete(svv[i][1]));
      dat.bounds.front().first=ToType<double>(ip->Interprete(svv[i][2]));
      dat.bounds.front().second=ToType<double>(ip->Interprete(svv[i][3]));
      flav = Flavour((kf_code)abs(crit1));
      if (crit1<0) flav = flav.Bar();
      dat.flavs.push_back(flav);
    }
    else if (keyword == string("DeltaEta")) {
      dat.type = 26;
      crit1=ToType<int>(ip->Interprete(svv[i][1]));
      crit2=ToType<int>(ip->Interprete(svv[i][2]));
      dat.bounds.front().first=ToType<double>(ip->Interprete(svv[i][3]));
      dat.bounds.front().second=ToType<double>(ip->Interprete(svv[i][4]));
      flav = Flavour((kf_code)abs(crit1));
      if (crit1<0) flav = flav.Bar();
      (dat.flavs).push_back(flav);
      flav = Flavour((kf_code)abs(crit2));
      if (crit2<0) flav = flav.Bar();
      dat.flavs.push_back(flav);
    }
        else if (keyword == string("DeltaPhi")) {
      dat.type = 27;
      crit1=ToType<int>(ip->Interprete(svv[i][1]));
      crit2=ToType<int>(ip->Interprete(svv[i][2]));
      dat.bounds.front().first=ToType<double>(ip->Interprete(svv[i][3]));
      dat.bounds.front().second=ToType<double>(ip->Interprete(svv[i][4]));
      flav = Flavour((kf_code)abs(crit1));
      if (crit1<0) flav = flav.Bar();
      (dat.flavs).push_back(flav);
      flav = Flavour((kf_code)abs(crit2));
      if (crit2<0) flav = flav.Bar();
      dat.flavs.push_back(flav);
    }
    else if (keyword == string("DeltaR")) {
      dat.type = 28;
      crit1=ToType<int>(ip->Interprete(svv[i][1]));
      crit2=ToType<int>(ip->Interprete(svv[i][2]));
      dat.bounds.front().first=ToType<double>(ip->Interprete(svv[i][3]));
      dat.bounds.front().second=ToType<double>(ip->Interprete(svv[i][4]));
      flav = Flavour((kf_code)abs(crit1));
      if (crit1<0) flav = flav.Bar();
      (dat.flavs).push_back(flav);
      flav = Flavour((kf_code)abs(crit2));
      if (crit2<0) flav = flav.Bar();
      dat.flavs.push_back(flav);
    }
    else if (keyword.find('"')==0 && 
	     keyword[keyword.length()-1]=='"') {
      dat.type=126;
      dat.helps=keyword.substr(1);
      dat.helps.erase(dat.helps.length()-1,1);
      dat.helps+="|"+(svv[i].size()>3?svv[i][3]:"");
      Data_Reader reader(",",";","!","=");
      reader.SetString(svv[i][1]);
      std::vector<int> flavs;
      if (!reader.VectorFromString(flavs,""))
 	THROW(critical_error,"Invalid Syntax in Selector.dat: '"+svv[i][1]+"'");
      dat.flavs.clear();
      for (size_t j(0);j<flavs.size();++j) {
 	flav=Flavour((kf_code)abs(flavs[j]));
 	if (flavs[j]<0) flav=flav.Bar();
 	dat.flavs.push_back(flav);
      }
      reader.SetString(svv[i][2]);
      std::vector<std::vector<double> > crits;
      if (!reader.MatrixFromString(crits,""))
	THROW(critical_error,"Invalid Syntax in Selector.dat: '"+svv[i][2]+"'");
      dat.bounds.clear();
      for (size_t j(0);j<crits.size();++j) {
	if (crits[j].size()<2) 
	  THROW(critical_error,
		"Invalid Syntax in Selector.dat: '"+svv[i][2]+"'");
	dat.bounds.push_back(std::pair<double,double>
			     (crits[j][0],crits[j][1]));
      }
    }
    else if (keyword.find("Bias")!=std::string::npos) {
      dat.type=0;
      if (keyword=="ET_Bias") dat.type=101;
      if (keyword=="PT_Bias") dat.type=102;
      if (keyword=="Eta_Bias") dat.type=105;
      if (keyword=="Delta_Eta_Bias") dat.type=113;
      if (keyword=="Delta_Phi_Bias") dat.type=114;
      if (keyword=="Delta_R_Bias") dat.type=115;
      if (keyword=="Mass_Bias") dat.type=116;
      if (dat.type>110) {
        crit1=ToType<int>(ip->Interprete(svv[i][1]));
        crit2=ToType<int>(ip->Interprete(svv[i][2]));
        std::string values=svv[i][3];
        dat.helps=svv[i][4];
	flav = Flavour((kf_code)abs(crit1));
	if (crit1<0) flav  = flav.Bar();
	flav2 = Flavour((kf_code)abs(crit2));
	if (crit2<0) flav2 = flav2.Bar();
	dat.flavs.push_back(flav);
	dat.flavs.push_back(flav2);
	Data_Reader reader(",",";","!","=");
	reader.SetString(values);
	std::vector<std::vector<double> > crits;
	if (!reader.MatrixFromString(crits,""))
	  THROW(critical_error,"Invalid Syntax in Selector.dat: '"+values+"'");
	dat.bounds.clear();
	for (size_t i(0);i<crits.size();++i) {
	  if (crits[i].size()<2) 
	    THROW(critical_error,"Invalid Syntax in Selector.dat: '"
		  +values+"'");
	  dat.bounds.push_back(std::pair<double,double>
			       (crits[i][0],crits[i][1]));
	}
      }
      else if (dat.type>0) {
        crit1=ToType<int>(ip->Interprete(svv[i][1]));
        std::string values=svv[i][2];
        dat.helps=svv[i][3];
	flav = Flavour((kf_code)abs(crit1));
	if (crit1<0) flav = flav.Bar();
	dat.flavs.push_back(flav);
	Data_Reader reader(",",";","!","=");
	reader.SetString(values);
	std::vector<std::vector<double> > crits;
	if (!reader.MatrixFromString(crits,""))
	  THROW(critical_error,"Invalid Syntax in Selector.dat: '"+values+"'");
	dat.bounds.clear();
	for (size_t i(0);i<crits.size();++i) {
	  if (crits[i].size()<1) 
	    THROW(critical_error,"Invalid Syntax in Selector.dat: '"
		  +values+"'");
	  dat.bounds.push_back(std::pair<double,double>
			       (crits[i].front(),crits[i].size()>1?
				crits[i][1]:rpa.gen.Ecms()));
	}
      }
    }
    if (dat.type>0) data.push_back(dat);
  }
  return 1;
}

void Selector_Data::ControlOutput() {
  if (data.size()<=0) {
    msg_Tracking()<<"Selector_Data empty."<<endl;
    return;
  }
  msg_Debugging()<<"Selector_Data : "<<endl;
  for (size_t i=0;i<data.size();i++) {
    switch (data[i].type) {
    case 1:  msg_Debugging()<<"Jet_Finder :      "; break;
    case 2:  msg_Debugging()<<"Cone_Finder:      "; break;
    case 3:  msg_Debugging()<<"Dipole_Jetfinder: "; break;
    case 11: msg_Debugging()<<"Energies   :      "; break;
    case 12: msg_Debugging()<<"PTs        :      "; break;
    case 13: msg_Debugging()<<"Rapidities :      "; break;
    case 14: msg_Debugging()<<"BeamAngles :      "; break;
    case 15: msg_Debugging()<<"ETs        :      "; break;
    case 16: msg_Debugging()<<"PseudoRaps :      "; break;  
    case 21: msg_Debugging()<<"Masses     :      "; break;
    case 22: msg_Debugging()<<"Angles     :      "; break;
    case 26: msg_Debugging()<<"DeltaEtaR  :      "; break;
    case 27: msg_Debugging()<<"DeltaPhi   :      "; break;
    case 28: msg_Debugging()<<"DeltaR     :      "; break;
    case 101: msg_Debugging()<<"ET_Bias    :      "; break;
    case 102: msg_Debugging()<<"PT_Bias    :      "; break;
    case 103: msg_Debugging()<<"Eta_Bias   :      "; break;
    case 113: msg_Debugging()<<"DEta_Bias  :      "; break;
    case 114: msg_Debugging()<<"DPhi_Bias  :      "; break;
    case 115: msg_Debugging()<<"DR_Bias    :      "; break;
    case 116: msg_Debugging()<<"Mass_Bias  :      "; break;
    case 126: msg_Debugging()<<"Variable   :      "; break;
    } 
    for (size_t j(0);j<data[i].bounds.size();++j) 
      msg_Debugging()<<"{"<<data[i].bounds[j].first<<","<<data[i].bounds[j].second<<"}";
    msg_Debugging()<<" : ";
    for (size_t j=0;j<(data[i].flavs).size();j++) msg_Debugging()<<(data[i]).flavs[j]<<" ";
    if (data[i].type == 14) msg_Debugging()<<" with "<<data[i].help;
    if (data[i].type==126) msg_Debugging()<<" -> "<<data[i].helps;
    msg_Debugging()<<endl;
  }
}

void Selector_Data::FillData(int i,int & type,std::vector<Flavour> & flavs,
			     std::vector<std::pair<double,double> > &bounds,int & help) {
  if ( (i<0) || (i>(int)data.size()) ) {
    msg_Error()<<"Error in Selector_Data::Data("<<i<<"). "
			  <<"Delimiter out of bounds."<<endl
			  <<"   Size : "<<data.size()<<endl;
    abort();
  }
  type  = data[i].type;
  flavs = data[i].flavs;
  bounds= data[i].bounds;
  help  = data[i].help;
}

void Selector_Data::SetData(int  _type,const std::vector<Flavour> & _flavs,
			    std::vector<std::pair<double,double> > &_bounds,int  _help) 
{
  RemoveData(_type);
  AddData(_type,_flavs,_bounds,_help);
}

void Selector_Data::AddData(int _type,const std::vector<Flavour> & _flavs,
			    std::vector<std::pair<double,double> > &_bounds,int _help) 
{
  Mom_Data dat;
  dat.type  = _type;
  dat.flavs = _flavs;
  dat.help  = _help;
  dat.bounds= _bounds;
  data.push_back(dat);
}

Mom_Data Selector_Data::RemoveData(int & _type) 
{
  Mom_Data last;
  for (std::vector<Mom_Data>::iterator it=data.begin();it!=data.end();++it) {
    if (it->type==_type) {
      last=*it;
      data.erase(it--);
    }
  }
  return last;
}


void Selector_Log::Output() { 
  msg_Info()<<"  Selector "<<m_name<<" rejection quota  : "
	    <<double(m_rejected)/double(m_rejected+m_passed)
	    <<"  ("<<m_rejected<<" / "<<m_passed+m_rejected<<")"<<std::endl;
}
