#include "SimpleXSecs.H"

#include "QCD_Processes.H"
#include "QED_Processes.H"

#include "Run_Parameter.H"
#include "Random.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace std;


SimpleXSecs::SimpleXSecs(std::string _path,std::string _file,
			 MODEL::Model_Base * _model) :
  m_path(_path), m_file(_file), XS_Group(0,0,NULL)
{
  m_atoms    = 1;
  p_dataread = new Data_Read(m_path+m_file);
}

SimpleXSecs::~SimpleXSecs() 
{
  if (p_dataread) { delete p_dataread; p_dataread = NULL; }
  if (p_seldata)  { delete p_seldata;  p_seldata  = NULL; }
  for(int i=m_xsecs.size();i>0;i--) {
    if (m_xsecs[i-1]) delete m_xsecs[i-1];
  }
}

bool SimpleXSecs::InitializeProcesses(BEAM::Beam_Spectra_Handler * _beam,
				      PDF::ISR_Handler * _isr) {
  p_beam         = _beam; 
  p_isr          = _isr;
  
  string xsfile  = p_dataread->GetValue<string>("XSFILE",string("XS.dat"));
  string selfile = p_dataread->GetValue<string>("XS-SELECTORFILE",string("XSSelector.dat"));
  p_seldata      = new Selector_Data(m_path+selfile);

  int    _scale_scheme   = p_dataread->GetValue<int>("SCALE SCHEME",0);
  int    _kfactor_scheme = p_dataread->GetValue<int>("KFACTOR SCHEME",0);
  double _scale_factor   = p_dataread->GetValue<double>("SCALE FACTOR",1.);
  ifstream from((m_path+xsfile).c_str());
  if (!from) {
    msg.Error()<<"Error in EXTRAXS::InitializeProcesses : "<<endl
	       <<"   File : "<<(m_path+xsfile).c_str()<<" not found ! Abort program execution."<<endl;
    abort();
  }

  string    buf;
  char      buffer[100];
  int       position;
  Flavour * flavs = new Flavour[4];
  for(;from;) {
    from.getline(buffer,100);
    if (buffer[0] != '%' && strlen(buffer)>0) {
      buf      = string(buffer);
      position = buf.find(string("All QCD 2->2")); 
      if (position>-1 && position<buf.length()) {
	flavs[0] = flavs[1] = flavs[2] = flavs[3] = Flavour(kf::jet);
	m_xsecs.push_back(new QCD_Processes(p_isr,p_beam,flavs,p_seldata,
					    _scale_scheme,_kfactor_scheme,_scale_factor));
	m_maxjet = 4;
	continue;
      }
      position = buf.find(string("All QED ee->qq")); 
      if (position>-1 && position<buf.length()) {
	if ((p_isr->Flav(0)== Flavour(kf::e)) &&
	    (p_isr->Flav(1)== Flavour(kf::e).Bar())) {
	  flavs[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::e);
	  flavs[1] = APHYTOOLS::Flavour(APHYTOOLS::kf::e).Bar();
	}
	else if ((p_isr->Flav(0)== Flavour(kf::e).Bar()) &&
		 (p_isr->Flav(1)== Flavour(kf::e))) {
	  flavs[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::e).Bar();
	  flavs[1] = APHYTOOLS::Flavour(APHYTOOLS::kf::e);
	}
	else {
	  msg.Error()<<"ERROR in EXTRAXS::QED_Processes() : "<<std::endl
		     <<"   Mismatch of flavours."<<endl;
	}
	flavs[2] = flavs[3] = Flavour(kf::jet);
	m_xsecs.push_back(new QED_Processes(p_isr,p_beam,flavs,p_seldata,
					    _scale_scheme,_kfactor_scheme,_scale_factor));
	m_maxjet = 2;
	continue;
      }
    }
  }
  delete [] flavs;
  from.close();

  if (m_xsecs.size()>0) return 1;
  msg.Error()<<"Error in SimpleXSecs::InitializeProcesses."<<endl
	     <<"   Did not find any process in "<<m_path+xsfile<<endl;
  return 0;
}

bool SimpleXSecs::CalculateTotalXSec() 
{
  bool okay = 1;
  for (int i=0;i<m_xsecs.size();i++) {
    okay = okay && m_xsecs[i]->CalculateTotalXSec();
    m_totalxs += m_xsecs[i]->Total();
  }
  msg.Events()<<"In SimpleXSecs::CalculateTotalXSec() = "
	      <<m_totalxs*AORGTOOLS::rpa.Picobarn()<<" pb."<<endl;
  return okay;
}

void  SimpleXSecs::SelectOne() {
  DeSelect();
  if (m_totalxs==0) p_selected = m_xsecs[int(ran.Get()*m_xsecs.size())];
  else {
    double disc = m_totalxs * ran.Get(); 
    for (int i=0;i<m_xsecs.size();i++) {
      disc -= m_xsecs[i]->Total();
      if (disc<0.) {
	p_selected = m_xsecs[i];
	p_selected->DeSelect();
	return;
      }
    }
    if (disc>0.) { 
      msg.Error()<<"Error in Process_Group::SelectOne() : "<<"Total xsec = "<<m_totalxs<<std::endl;
      return;
    }
  }
}

bool SimpleXSecs::UnweightedEvent(int mode)
{
  SelectOne();
  if (p_selected->OneEvent()) return 1;
  return 0;
}

double SimpleXSecs::WeightedEvent(int mode)
{
  SelectOne();
  double weight=p_selected->OneEvent();
  if (weight>0.) return weight;
  return 0.;
}

bool SimpleXSecs::PrepareXSecTables() {}
bool SimpleXSecs::LookUpXSec(double,bool,std::string) {}
void SimpleXSecs::SingleEvents() {}
