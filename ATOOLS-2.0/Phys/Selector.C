#include "Selector.H"
#include "Message.H"

#include "Algebra_Interpreter.H"
#include "MyStrStream.H"

using namespace ATOOLS;
using namespace std;


Selector_Base::~Selector_Base() { }

void Selector_Base::Output() { 
  if (!(msg.LevelIsTracking())) return;
  if(m_sel_log) {
    m_sel_log->Output();
    msg.Out()<<m_name<<"  total number of rejections: "<<m_sel_log->Rejections()<<std::endl;
  }
}

void Selector_Base::Add(Selector_Base *) {
  msg.Error()<<"Selector_Base::Add : Virtual method."<<std::endl;
}

double * Selector_Base::ActualValue() {
  msg.Error()<<"Selector_Base::ActualValue :"
	     <<m_name<<" Virtual method."<<std::endl;
  return 0;
}

void Selector_Base::BuildCuts(Cut_Data *) { 
  msg.Error()<<"Selector_Base::BuildCuts : Virtual method."<<std::endl;
}

void Selector_Base::UpdateCuts(double,double,Cut_Data *) { 
  msg.Error()<<"Selector_Base::BuildCuts : Virtual method."<<std::endl;
}

void Selector_Base::SetRange(std::vector<Flavour>,double,double) { 
  msg.Error()<<"Selector_Base::SetRange : Virtual method."<<std::endl;
}

void Selector_Base::SetRange(std::vector<Flavour>,int,double,double) { 
  msg.Error()<<"Selector_Base::SetRange : Virtual method."<<std::endl;
}

bool Selector_Base::GetValue(const std::string &name,double &value)
{
  msg.Error()<<"Selector_Base::GetValue("<<name<<",..): "
	     <<"Virtual method called."<<std::endl;
  return false;
}

int    Selector_Base::NeedUpdate()                         { return 0; }
int    Selector_Base::IsConditional()                      { return 0; }
void   Selector_Base::SetSRange(double _smin,double _smax) { m_smin = _smin; m_smax = _smax; }
void   Selector_Base::SetName(std::string _name)           { m_name = _name; }
//double Selector_Base::Smin()                               { return m_smin; }
//double Selector_Base::Smax()                               { return m_smax; }
std::string Selector_Base::Name()                          { return m_name; }


/*-----------------------------------------------------------------------------------

  Selector_Data

  -----------------------------------------------------------------------------------*/


Selector_Data::Selector_Data() {}

Selector_Data::Selector_Data(std::string path) {
  if (!ReadInData(path)) {
    msg.Error()<<"Error in Selector_Data::Selector_Data("<<path<<")."<<endl
			  <<"Cannot initialise any selector. Abort."<<endl;
    abort();
  }
  ControlOutput();
}

bool Selector_Data::ReadInData(std::string filename) {
  ifstream from(filename.c_str());
  if (!from) {
    msg.Error()<<"Error in Selector_Data::ReadInData("<<filename<<"). "
			  <<"File does not exist."<<endl;
    return 0;
  }
  
  std::string keyword;
  Mom_Data    dat;
  int         crit1,crit2;
  Flavour     flav;
  for(;from;) {
    dat.flavs.clear();
    keyword=string("");
    from>>keyword;
    if (keyword == string("JetFinder")) {
      dat.type = 1;
      std::string dmin, dmax;
      from>>dmin>>dmax;
      Algebra_Interpreter inter;
      dat.min=ToType<double>(inter.Interprete(dmin));
      dat.max=ToType<double>(inter.Interprete(dmax));
      data.push_back(dat);
    }
    if (keyword == string("ConeFinder")) {
      dat.type = 2;
      from>>dat.min;
      data.push_back(dat);
    }
    if (keyword == string("Energy")) {
      dat.type = 11;
      from>>crit1>>dat.min>>dat.max;
      flav = Flavour(kf::code(abs(crit1)));
      if (crit1<0) flav = flav.Bar();
      (dat.flavs).push_back(flav);
      data.push_back(dat);
    }
    if (keyword == string("PT")) {
      dat.type = 12;
      from>>crit1>>dat.min>>dat.max;
      flav = Flavour(kf::code(abs(crit1)));
      if (crit1<0) flav = flav.Bar();
      (dat.flavs).push_back(flav);
      data.push_back(dat);
    }
    if (keyword == string("BFKL_PT")) {
      dat.type = 25;
      std::string dmin, dmax;
      from>>crit1>>dmin>>dmax;
      Algebra_Interpreter inter;
      dat.min=ToType<double>(inter.Interprete(dmin));
      dat.max=ToType<double>(inter.Interprete(dmax));
      flav = Flavour(kf::code(abs(crit1)));
      if (crit1<0) flav = flav.Bar();
      (dat.flavs).push_back(flav);
      data.push_back(dat);
    }
    if (keyword == string("X")) {
      dat.type = 24;
      std::string dmin, dmax;
      from>>dmin>>dmax;
      Algebra_Interpreter inter;
      dat.min=ToType<double>(inter.Interprete(dmin));
      dat.max=ToType<double>(inter.Interprete(dmax));
      data.push_back(dat);
    }
    if (keyword == string("PT2")) {
      dat.type = 23;
      from>>crit1>>crit2>>dat.min>>dat.max;
      Flavour flav1 = Flavour(kf::code(abs(crit1)));
      Flavour flav2 = Flavour(kf::code(abs(crit2)));
      if (crit1<0) flav1 = flav1.Bar();
      if (crit2<0) flav2 = flav2.Bar();
      (dat.flavs).push_back(flav1);
      (dat.flavs).push_back(flav2);
      data.push_back(dat);
    }
    if (keyword == string("Rapidity")) {
      dat.type = 13;
      from>>crit1>>dat.min>>dat.max;
      flav = Flavour(kf::code(abs(crit1)));
      if (crit1<0) flav = flav.Bar();
      (dat.flavs).push_back(flav);
      data.push_back(dat);
    }
    if (keyword == string("BeamAngle")) {
      dat.type = 14;
      from>>crit1>>crit2>>dat.min>>dat.max;
      flav = Flavour(kf::code(abs(crit1)));
      if (crit1<0) flav = flav.Bar();
      (dat.flavs).push_back(flav);
      dat.help = crit2;
      data.push_back(dat);
    }  
    if (keyword == string("ET")) {
      dat.type = 15;
      from>>crit1>>dat.min>>dat.max;
      flav = Flavour(kf::code(abs(crit1)));
      if (crit1<0) flav = flav.Bar();
      (dat.flavs).push_back(flav);
      data.push_back(dat);
    }
    if (keyword == string("PseudoRapidity")) {
      dat.type = 16;
      from>>crit1>>dat.min>>dat.max;
      flav = Flavour(kf::code(abs(crit1)));
      if (crit1<0) flav = flav.Bar();
      (dat.flavs).push_back(flav);
      data.push_back(dat);
    }
    if (keyword == string("Mass")) {
      dat.type = 21;
      from>>crit1>>crit2>>dat.min>>dat.max;
      flav = Flavour(kf::code(abs(crit1)));
      if (crit1<0) flav = flav.Bar();
      (dat.flavs).push_back(flav);
      flav = Flavour(kf::code(abs(crit2)));
      if (crit2<0) flav = flav.Bar();
      (dat.flavs).push_back(flav);
      data.push_back(dat);
    }
    if (keyword == string("Angle")) {
      dat.type = 22;
      from>>crit1>>crit2>>dat.min>>dat.max;
      flav = Flavour(kf::code(abs(crit1)));
      if (crit1<0) flav = flav.Bar();
      (dat.flavs).push_back(flav);
      flav = Flavour(kf::code(abs(crit2)));
      if (crit2<0) flav = flav.Bar();
      (dat.flavs).push_back(flav);
      data.push_back(dat);
    }  
    if (keyword == string("SummedPT")) {
      dat.type = 32;
      from>>crit1>>dat.min>>dat.max;
      flav = Flavour(kf::code(abs(crit1)));
      if (crit1<0) flav = flav.Bar();
      (dat.flavs).push_back(flav);
      data.push_back(dat);
    }
  }
  from.close();
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
    case 1:  msg_Debugging()<<"Jet_Finder : "; break;
    case 2:  msg_Debugging()<<"Cone_Finder: "; break;
    case 11: msg_Debugging()<<"Energies   : "; break;
    case 12: msg_Debugging()<<"PTs        : "; break;
    case 13: msg_Debugging()<<"Rapidities : "; break;
    case 14: msg_Debugging()<<"BeamAngles : "; break;
    case 15: msg_Debugging()<<"ETs        : "; break;
    case 16: msg_Debugging()<<"PseudoRaps : "; break;  
    case 21: msg_Debugging()<<"Masses     : "; break;
    case 22: msg_Debugging()<<"Angles     : "; break;
    case 32: msg_Debugging()<<"SummedPT   : "; break;
    } 
    msg_Debugging()<<data[i].min<<" ... "<<data[i].max<<" : ";
    for (size_t j=0;j<(data[i].flavs).size();j++) msg_Debugging()<<(data[i]).flavs[j]<<" ";
    if (data[i].type == 14) msg_Debugging()<<" with "<<data[i].help;
    msg_Debugging()<<endl;
  }
}

void Selector_Data::Data(int i,int & type,std::vector<Flavour> & flavs,
			 int & help,double & min,double & max) {
  if ( (i<0) || (i>(int)data.size()) ) {
    msg.Error()<<"Error in Selector_Data::Data("<<i<<"). "
			  <<"Delimiter out of bounds."<<endl
			  <<"   Size : "<<data.size()<<endl;
    abort();
  }
  type  = data[i].type;
  flavs = data[i].flavs;
  min   = data[i].min;
  max   = data[i].max;
  help  = data[i].help;
}

void Selector_Data::SetData(int  _type,const std::vector<Flavour> & _flavs,
			    int  _help,double  _min,double  _max) 
{
  RemoveData(_type);
  AddData(_type,_flavs,_help,_min,_max);
}

void Selector_Data::AddData(int _type,const std::vector<Flavour> & _flavs,
			    int _help,double _min,double _max) 
{
  Mom_Data dat;
  dat.type  = _type;
  dat.flavs = _flavs;
  dat.help  = _help;
  dat.min   = _min;
  dat.max   = _max;
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
