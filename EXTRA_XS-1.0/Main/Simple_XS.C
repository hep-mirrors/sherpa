#include "Simple_XS.H"

#include "QCD_Processes.H"
#include "QED_Processes.H"

#include "Run_Parameter.H"
#include "Random.H"
#include "MyStrStream.H"
#include "Data_Read.H"
#include "Model_Base.H"
#include "XS_Selector.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;


Simple_XS::Simple_XS(const std::string &path,const std::string &file,
			 MODEL::Model_Base *const model):
  XS_Group(0,0,NULL),
  m_path(path), 
  m_file(file),
  m_maxjet(2)
{
  m_atoms=1;
  p_dataread = new Data_Read(m_path+m_file);
}

Simple_XS::~Simple_XS() 
{
  if (p_selectordata) delete p_selectordata;
  delete p_dataread; 
}

bool Simple_XS::InitializeProcesses(BEAM::Beam_Spectra_Handler *const beamhandler,
				    PDF::ISR_Handler *const isrhandler,const bool construct) 
{
  p_beamhandler=beamhandler; 
  p_isrhandler=isrhandler;
  std::string processfile=p_dataread->GetValue<std::string>("PROCESS_FILE",
							    std::string("Processes.dat"));
  std::string selectorfile=p_dataread->GetValue<std::string>("SELECTOR_FILE",
							       std::string("Selector.dat"));
  p_selectordata = new Selector_Data(m_path+selectorfile);
  m_scalescheme=p_dataread->GetValue<int>("SCALE_SCHEME",0);
  m_kfactorscheme=p_dataread->GetValue<int>("KFACTOR_SCHEME",0);
  m_scalefactor=p_dataread->GetValue<double>("SCALE_FACTOR",1.);
  if (!construct) return true;
  ifstream from((m_path+processfile).c_str());
  if (!from) {
    msg.Error()<<"Error in Simple_XS::InitializeProcesses : "<<endl
	       <<"   File : "<<(m_path+processfile).c_str()<<" not found ! Abort program execution."<<endl;
    exit(145);
  }
  char buffer[100];
  size_t position;
  int         flag;
  string      buf,ini,fin;
  int         nIS,   nFS;
  Flavour   * IS,  * FS, * flavs;
  int         order_ew,order_strong; 
  while(from) {
    from.getline(buffer,100);
    if (buffer[0] != '%' && strlen(buffer)>0) {
      buf        = string(buffer);
      position   = buf.find(string("Process :")); 
      flag       = 0;
      if (position!=std::string::npos && position<buf.length()) {
	flag     = 1;
	buf      = buf.substr(position+9);
	position = buf.find(string("->"));
	if (position > 0) {
	  ini    = buf.substr(0,position);
	  fin    = buf.substr(position+2);
	  nIS    = ExtractFlavours(IS,ini);
	  nFS    = ExtractFlavours(FS,fin);
	  if (!(p_isrhandler->CheckConsistency(IS))) {
	    msg.Error()<<"Error in initialising ISR_Handler."<<endl
		       <<" "<<p_isrhandler->Flav(0)<<" -> "<<IS[0]<<", "
		       <<" "<<p_isrhandler->Flav(1)<<" -> "<<IS[1]<<endl
		       <<"  Delete it and ignore the process."<<endl;
	    flag = 0;
	  }
	  if (flag==0) {
	    delete [] IS;
	    delete [] FS;
	  }
	  else {
	    do {
	      from.getline(buffer,100);
	      if (buffer[0] != '%' && strlen(buffer)>0) {
		buf      = string(buffer);
		position = buf.find(string("End process"));
		if (!from) {
		  msg.Error()<<"Error in Amegic::InitializeProcesses("<<m_path+processfile<<")."<<endl
			     <<"   End of file reached without 'End process'-tag."<<endl
			     <<"   Continue and hope for the best."<<endl;
		  position     = 0;
		}
	      }
	    }
	    while (position==std::string::npos);
	    msg.Debugging()<<"Read in process :";
	    for (short int i=0;i<nIS;i++) msg.Debugging()<<" "<<IS[i].Name();
	    msg.Debugging()<<" -> ";
	    for (short int i=0;i<nFS;i++) msg.Debugging()<<FS[i].Name()<<" ";
	    msg.Debugging()<<" EW("<<order_ew<<"), QCD("<<order_strong<<")"<<endl;
	    if (nIS+nFS>(int)m_nmax) m_nmax = nIS+nFS;
	    flavs              = new Flavour[nIS+nFS];
	    for (int i=0;i<nIS;i++) flavs[i]     = IS[i]; 
	    for (int i=0;i<nFS;i++) flavs[i+nIS] = FS[i]; 
	    bool single        = 1;
	    for (int i=0;i<nIS+nFS;i++) {
	      if (flavs[i].Size()>1) { single = 0; break; }
	    } 
	    double summass = 0.;
	    for (int i=0;i<nFS;i++) summass += flavs[i+nIS].Mass();
	    if (summass<rpa.gen.Ecms()) {
	      if (single) {
 		XS_Group *group=FindGroup(nIS,nFS,flavs,m_scalescheme,m_kfactorscheme,m_scalefactor);
		group->Add(group->XSSelector()->GetXS(nIS,nFS,flavs,p_isrhandler->KMROn()));
		p_selected=group;
	      }
	    }
	    else {
	      msg.Out()<<"Kicked Process: ";
	      for (short int i=0;i<nIS;i++) msg.Out()<<" "<<IS[i].Name();
	      msg.Out()<<" -> ";
	      for (short int i=0;i<nFS;i++) msg.Out()<<FS[i].Name()<<" ";
	      msg.Out()<<", kinematically not allowed."<<endl;
	    }	    
	    delete [] flavs;
	  }
	}
      }
    }
  }
  if (m_xsecs.size()>0) return true;
  msg.Error()<<"Simple_XS::InitializeProcesses("<<beamhandler<<","<<isrhandler<<"):"<<std::endl
	     <<"   Did not find any process in '"<<m_path+processfile<<"' !"<<std::endl;
  return false;
}

XS_Group *Simple_XS::FindGroup(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours,
			       const int scalescheme,const int kfactorscheme, 
			       const double scalefactor)
{
  ATOOLS::Flavour *dummy = new ATOOLS::Flavour[nin+nout];
  for (size_t i=0;i<nin+nout;++i) dummy[i]=flavours[i];
  for (size_t i=0;i<m_xsecs.size();++i) {
    if (nin==m_xsecs[i]->NIn() && nout==m_xsecs[i]->NOut()) {
      const ATOOLS::Flavour *test=m_xsecs[i]->Flavours();
      bool hit=true;
      for (size_t j=0; j<nin+nout;++j) if (test[j]!=dummy[j]) hit=false;
      if (hit) {
	delete [] dummy;
	return dynamic_cast<XS_Group*>(m_xsecs[i]);
      }
    }
  }
  XS_Group *newgroup = new XS_Group(nin,nout,dummy,scalescheme,kfactorscheme,scalefactor,
				    p_beamhandler,p_isrhandler,p_selectordata);
  m_xsecs.push_back(newgroup);
  delete [] dummy;
  return newgroup;
}

int Simple_XS::ExtractFlavours(ATOOLS::Flavour *&flavours,std::string buffer)
{
  int ii[20];
  char pc[20],pp[20];
  double pd[20],angle[20];
  int count = 0;
  buffer += string(" "); 
  short int i;
  for (i=0;i<20;i++) { pp[i]=' '; pd[i]=0.; pc[i]=' '; angle[i]=0.;}
  while(buffer.length()>1) {
    if (buffer[0]==' ') buffer = buffer.substr(1);
    else {
      int next = buffer.find(string(" "));
      string pn = buffer.substr(0,next);
      size_t nxt = pn.find(string("("));
      string number;
      if (nxt==string::npos) number = pn.substr(0,pn.length());
      else {
	number = pn.substr(0,nxt);
	pn.erase(0,nxt);
	if(pn[pn.length()-1]==')'){
	  pn.erase(0,1);
	  pn.erase(pn.length()-1,1);
	  size_t lh=pn.find(string("l"));
	  if(lh!=string::npos) {
	    pc[count]='l';
	    pp[count]='+';
	    string ha = pn.substr(lh+1,pn.length()-lh);
	    MyStrStream astream;
	    astream<<ha;
	    astream>>angle[count];
	  }
	  else {
	    pp[count]=pn[pn.length()-1];
	    pc[count]=pn[pn.length()-2];
	    pn.erase(pn.length()-2,2);
	  }
	  MyStrStream pstream;
	  pstream<<pn;
	  pstream>>pd[count];	  
	}
      }
      MyStrStream sstream;
      sstream<<number;
      sstream>>ii[count];
      count++;
      buffer = buffer.substr(next);
    }
  }
  flavours = new Flavour[count];
  for (i=0;i<count;i++) {
    flavours[i] = Flavour(kf::code(iabs(ii[i])));
    if (ii[i]<0) flavours[i] = flavours[i].Bar();
  }
  return count;
}

bool Simple_XS::CalculateTotalXSec(const std::string &resultpath) 
{
  Reset();
  bool okay = 1;
  for (size_t i=0;i<m_xsecs.size();++i) {
    okay = okay && m_xsecs[i]->CalculateTotalXSec(resultpath);
    m_totalxs += m_xsecs[i]->TotalXS();
  }
  msg.Events()<<"In Simple_XS::CalculateTotalXSec() = "
	      <<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb."<<endl;
  return okay;
}

void  Simple_XS::SelectOne() {
  DeSelect();
  if (m_totalxs==0) p_selected = m_xsecs[int(ran.Get()*m_xsecs.size())];
  else {
    double disc = m_totalxs * ran.Get(); 
    for (size_t i=0;i<m_xsecs.size();++i) {
      disc -= m_xsecs[i]->TotalXS();
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

bool Simple_XS::OneEvent(const int mode)
{
  SelectOne();
  if (p_selected->OneEvent()) return 1;
  return 0;
}

ATOOLS::Blob_Data_Base *Simple_XS::WeightedEvent(const int mode)
{
  SelectOne();
  return p_selected->WeightedEvent();
}

void Simple_XS::ResetSelector(ATOOLS::Selector_Data *_seldata)
{
  for (unsigned int i=0;i<m_xsecs.size();++i) m_xsecs[i]->ResetSelector(_seldata);
}

bool Simple_XS::PrepareXSecTables() 
{ 
  return 0; 
}

bool Simple_XS::LookUpXSec(double,bool,std::string) 
{ 
  return 0; 
}

void Simple_XS::SingleEvents() 
{
}

