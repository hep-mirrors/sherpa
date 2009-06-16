#include "SHERPA/HerwigTools/MCatNLO_Interface.H"

#include "SHERPA/HerwigTools/MCatNLO_Wrapper.H"
#include "SHERPA/HerwigTools/Herwig_Wrapper.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Math/Scaling.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

std::ostream& SHERPA::operator<<(std::ostream& ostr,const ptp::code ptpc)
{
  switch (ptpc) {
  case (ptp::Unspecified) : return ostr<<"Unspecified. ";
  case (ptp::Single_Boson): return ostr<<"Single_Boson.";
  case (ptp::Lepton_Pair) : return ostr<<"Lepton_Pair  ";
  case (ptp::Higgs)       : return ostr<<"Higgs.       ";
  case (ptp::Heavy_Quark) : return ostr<<"Heavy_Quark. ";
  case (ptp::Boson_Pair)  : return ostr<<"Boson_Pair.  ";
  default: return ostr<<"Not implemented yet!"<<std::endl;
  }
}

MCatNLO_Interface::MCatNLO_Interface(std::string _m_path,std::string _m_file,bool sherpa) :
  m_path(_m_path),m_file(_m_file), m_proc(ptp::Unspecified)
{ 
  p_herwig = new Herwig_Interface(m_path,m_file,sherpa);
  if (!ReadInTheParameters()) {
    msg_Error()<<"Error in MCatNLO_Interface::MCatNLO_Interface: "<<std::endl
	       <<"   Did not succeed to initialize the Interface, abort the run."<<std::endl;
    abort();
  }
}

MCatNLO_Interface::~MCatNLO_Interface()
{
  if (p_herwig) { delete p_herwig; p_herwig = NULL; }
}

bool MCatNLO_Interface::ReadInTheParameters()
{
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader(" ",";","!","=");
  reader->AddWordSeparator("\t");
  reader->SetInputPath(m_path);
  reader->SetInputFile(m_file);
  reader->AddIgnore("(");
  reader->AddIgnore(")");
  reader->AddIgnore(",");

  std::string nlo_events;
  if (!reader->ReadFromFile(nlo_events,"NLO_EVENTFILE"))  {
    msg_Error()<<"Error in MCatNLO_Interface::ReadInTheParameters("<<m_path<<"/"<<m_file<<"):"<<std::endl
	       <<"   No value found for token NLO_EVENTFILE. Return false."<<std::endl;
    return false;
  }
  
  ifstream from(nlo_events.c_str());
  if (!from) {
    msg_Error()<<"Error in MCatNLO_Interface::ReadInTheParameters("<<m_path<<"/"<<m_file<<"):"<<std::endl
	       <<"   No NLO event-file found : "<<nlo_events<<" not in run-directory."<<std::endl
	       <<"   Return false."<<std::endl;
    return false;
  }
  char   buffer[200];
  string buf, act, pdf;
  int    pos, proc, help, numb, ipdf;
  double ecm, fac[4];

  // c.m. energy and scale choices.
  from.getline(buffer,200);
  buf  = string(buffer);
  buf  = buf.substr(0,buf.find(string("-->")));
  Shorten(buf);
  pos  = buf.find(string(" "));
  act  = buf.substr(0,pos);
  Shorten(act);
  ecm = atof(act.c_str());
  if (!IsEqual(ecm,rpa.gen.Ecms())) {
    msg_Error()<<"Error in MCatNLO_Interface::ReadInTheParameters("<<m_path<<"/"<<m_file<<"):"<<std::endl
	       <<"   C.m. energies of file and run do not coincide : "
	       <<ecm<<" vs. "<<rpa.gen.Ecms()<<std::endl
	       <<"   Return false."<<std::endl;
    return false;
  }
  for (int i=0;i<4;i++) {
    buf = buf.substr(pos);
    Shorten(buf);
    pos = buf.find(string(" "));
    act = buf.substr(0,pos);
    Shorten(act);
    fac[i] = atof(act.c_str());
  }
  // process code.
  from.getline(buffer,200);
  buf  = string(buffer);
  buf  = buf.substr(0,buf.find(string("-->")));
  Shorten(buf);
  pos  = buf.find(string(" "));
  act  = buf.substr(0,pos);
  Shorten(act);
  proc = atoi(act.c_str());
  help = proc/10;
  switch (help) {
  case 135:
  case 136:
  case 137:
  case 146:
  case 147:
    hwproc.iproc = -proc;
    m_proc       = ptp::Lepton_Pair; 
    break;
  case  285:
  case 1285: 
  case  286:
  case  287:
  case  288:
    hwproc.iproc = -proc;
    m_proc       = ptp::Boson_Pair; 
    break;
  default: 
    msg_Error()<<"Error in MCatNLO_Interface::ReadInTheParameters("<<m_path<<"/"<<m_file<<"):"<<std::endl
	       <<"   Process number not understood : "<<proc<<std::endl
	       <<"   Return false."<<std::endl;
    return false;    
  }
  // masses etc. - ignore them.
  from.getline(buffer,200);
  from.getline(buffer,200);
  // beams
  from.getline(buffer,200);
  buf  = string(buffer);
  buf  = buf.substr(0,buf.find(string("-->")));
  Shorten(buf);
  pos  = buf.find(string(" "));
  act  = buf.substr(0,pos);
  Shorten(act);
  buf  = buf.substr(pos);
  Shorten(buf);
  if (!(act==string("P") || act==string("PBAR")) ||
      !(buf==string("P") || buf==string("PBAR")) ) {
    msg_Error()<<"Error in MCatNLO_Interface::ReadInTheParameters("<<m_path<<"/"<<m_file<<"):"<<std::endl
	       <<"   Funny beam particles (1,2), not understood : ("<<act<<","<<buf<<")"<<std::endl
	       <<"   Return false."<<std::endl;
    return false;    
  } 
  // pdf 
  from.getline(buffer,200);
  buf  = string(buffer);
  buf  = buf.substr(0,buf.find(string("-->")));
  Shorten(buf);
  pos  = buf.find(string(" "));
  pdf  = buf.substr(0,pos);
  Shorten(pdf);
  buf  = buf.substr(pos);
  Shorten(buf);
  ipdf = atoi(buf.c_str());
  MakeFortranString(hwprch.autpdf[0],pdf,20);
  MakeFortranString(hwprch.autpdf[1],pdf,20);
  hwpram.modpdf[0] = hwpram.modpdf[1] = ipdf;
  // alpha s choice and scheme, ignore them
  from.getline(buffer,200);
  ctmplam.tmplam = -1.;
  // number of events
  from.getline(buffer,200);
  buf  = string(buffer);
  buf  = buf.substr(0,buf.find(string("-->")));
  Shorten(buf);
  numb = atoi(buf.c_str());
  if (numb<rpa.gen.NumberOfEvents()) {
    msg_Error()<<"Error in MCatNLO_Interface::ReadInTheParameters("<<m_path<<"/"<<m_file<<"):"<<std::endl
	       <<"   Not enough NLO events in event file : "<<numb<<", "
	       <<rpa.gen.NumberOfEvents()<<" wanted."<<std::endl
	       <<"   Continue anyhow."<<std::endl;
  }
  MakeFortranString(vvjin.qqin,nlo_events,50);
  from.close();
  delete reader;
  return true;
}

bool MCatNLO_Interface::Initialize() {
  mccinit(1);
  return true;
}

void MCatNLO_Interface::Shorten(std::string& str) {
  //kill initial spaces
  for (;;) {    
    if (int(str[0])==32 || int(str[0])==9) str = str.substr(1);
    else break;
  }
  //kill final spaces
  for (;;) {    
    if (int(str[str.length()-1])==32 ||
	//Tabulator
	int(str[str.length()-1])==9) str = str.substr(0,str.length()-1);
    else break;
  }
}

bool MCatNLO_Interface::OneEvent(ATOOLS::Blob_List *const bloblist,double &weight)
{
  bool test = p_herwig->OneEvent(bloblist,weight);
  return test;
}

void MCatNLO_Interface::Terminate()
{
  mcctrm();   
}
 
