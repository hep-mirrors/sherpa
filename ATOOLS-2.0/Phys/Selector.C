#include "Selector.H"
#include "Message.H"
#include "Run_Parameter.H"

using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;


Selector_Base::~Selector_Base() { }

void Selector_Base::Output() { 
  if (!(AORGTOOLS::rpa.gen.Debugging())) return;
  if(m_sel_log) {
    m_sel_log->Output();
    AORGTOOLS::msg.Out()<<m_name<<"  total number of rejections: "<<m_sel_log->Rejections()<<std::endl;
  }
}

/*-----------------------------------------------------------------------------------

  Selector_Data

  -----------------------------------------------------------------------------------*/


Selector_Data::Selector_Data(std::string path) {
  if (!ReadInData(path)) {
    AORGTOOLS::msg.Error()<<"Error in Selector_Data::Selector_Data("<<path<<")."<<endl
			  <<"Cannot initialise any selector. Abort."<<endl;
    abort();
  }
  ControlOutput();
}

bool Selector_Data::ReadInData(std::string filename) {
  ifstream from(filename.c_str());
  if (!from) {
    AORGTOOLS::msg.Error()<<"Error in Selector_Data::ReadInData("<<filename<<"). "
			  <<"File does not exist."<<endl;
    return 0;
  }
  
  std::string keyword;
  Mom_Data    dat;
  int         crit1,crit2;
  Flavour     flav;
  for(;from;) {
    dat.flavs.erase(dat.flavs.begin(),dat.flavs.end());
    keyword=string("");
    from>>keyword;
    if (keyword == string("JetFinder")) {
      dat.type = 1;
      from>>dat.min>>dat.max;
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
    AORGTOOLS::msg.Debugging()<<"Selector_Data empty."<<endl;
    return;
  }
  AORGTOOLS::msg.Debugging()<<"Selector_Data : "<<endl;
  for (int i=0;i<data.size();i++) {
    switch (data[i].type) {
    case 1:  AORGTOOLS::msg.Debugging()<<"Jet_Finder : "; break;
    case 2:  AORGTOOLS::msg.Debugging()<<"Cone_Finder: "; break;
    case 11: AORGTOOLS::msg.Debugging()<<"Energies   : "; break;
    case 12: AORGTOOLS::msg.Debugging()<<"PTs        : "; break;
    case 13: AORGTOOLS::msg.Debugging()<<"Rapidities : "; break;
    case 14: AORGTOOLS::msg.Debugging()<<"BeamAngles : "; break;
    case 15: AORGTOOLS::msg.Debugging()<<"ETs        : "; break;
    case 16: AORGTOOLS::msg.Debugging()<<"PseudoRaps : "; break;  
    case 21: AORGTOOLS::msg.Debugging()<<"Masses     : "; break;
    case 22: AORGTOOLS::msg.Debugging()<<"Angles     : "; break;
    case 32: AORGTOOLS::msg.Debugging()<<"SummedPT   : "; break;
    } 
    AORGTOOLS::msg.Debugging()<<data[i].min<<" ... "<<data[i].max<<" : ";
    for (int j=0;j<(data[i].flavs).size();j++) AORGTOOLS::msg.Debugging()<<(data[i]).flavs[j]<<" ";
    if (data[i].type == 14) AORGTOOLS::msg.Debugging()<<" with "<<data[i].help;
    AORGTOOLS::msg.Debugging()<<endl;
  }
}

void Selector_Data::Data(int i,int & type,std::vector<APHYTOOLS::Flavour> & flavs,
			 int & help,double & min,double & max) {
  if ( (i<0) || (i>data.size()) ) {
    AORGTOOLS::msg.Error()<<"Error in Selector_Data::Data("<<i<<"). "
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
