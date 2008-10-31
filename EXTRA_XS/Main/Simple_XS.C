#include "Simple_XS.H"

#include "ISR_Handler.H"
#include "Run_Parameter.H"
#include "Random.H"
#include "MyStrStream.H"
#include "Data_Reader.H"
#include "Data_Reader.H"
#include "Model_Base.H"
#include "XS_Selector.H"
#include "Regulator_Base.H"
#include "Remnant_Base.H"
#include "Phase_Space_Handler.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Simple_XS::Simple_XS(const std::string &path,const std::string &file,
		     MODEL::Model_Base *const model):
  XS_Group(0,0,NULL,NULL),
  m_path(path), 
  m_file(file),
  m_minqcdjet(99), m_maxqcdjet(0), m_maxjet(2)
{
  m_nmax=0;
  m_atoms=1;
  p_dataread = new Data_Reader(" ",";","!","=");
  p_dataread->SetWordSeparator("\t");
  p_dataread->SetInputPath(m_path);
  p_dataread->SetInputFile(m_file);
  std::string modelfile=rpa.gen.Variable("MODELFILE");
  InitializeModel(model,m_path+modelfile);
}

Simple_XS::~Simple_XS() 
{
  if (p_selectordata) delete p_selectordata;
  delete p_dataread; 
}

bool Simple_XS::InitializeProcesses(BEAM::Beam_Spectra_Handler *const beamhandler,
				    PDF::ISR_Handler *const isrhandler,
				    const bool construct) 
{
  p_beamhandler=beamhandler; 
  p_isrhandler=isrhandler;
  p_remnants[0]=p_isrhandler->GetRemnant(0);
  p_remnants[1]=p_isrhandler->GetRemnant(1);
  if (p_remnants[0]==NULL || p_remnants[1]==NULL) {
    msg_Error()<<"Simple_XS::Simple_XS(..): "
		       <<"No beam remnant handler found."<<std::endl;
  }
  std::string processfile=rpa.gen.Variable("PROCESSFILE");
  std::string selectorfile=rpa.gen.Variable("SELECTORFILE");
  if (construct) p_selectordata = new Selector_Data(m_path,selectorfile);
  else p_selectordata = new Selector_Data();
  p_dataread->SetTags(Integrable_Base::ScaleTags());
  m_scalescheme=(PHASIC::scl::scheme)
    p_dataread->GetValue<int>("SCALE_SCHEME",1);
  p_dataread->SetTags(std::map<std::string,std::string>());
  m_muf2tag=p_dataread->GetValue<std::string>("FACTORIZATION_SCALE","");
  m_mur2tag=p_dataread->GetValue<std::string>("RENORMALIZATION_SCALE","");
  m_kfactorscheme=p_dataread->GetValue<int>("KFACTOR_SCHEME",1);
  double fix_scale=p_dataread->
    GetValue<double>("FIXED_SCALE",sqr(rpa.gen.Ecms()));
  int regulate=p_dataread->GetValue<int>("REGULATE_XS",0);
  if (regulate>0) {
    m_regulator=p_dataread->GetValue<std::string>
      ("XS_REGULATOR",std::string("QCD_Trivial"));
    double param=p_dataread->GetValue<double>("XS_REGULATION",0.71);
    m_regulation.push_back(param);
  }
  if (!construct) return true;
  My_In_File dummy(m_path, processfile);
  if (!dummy.Open()) {
    THROW(critical_error,"Can't open file '"+m_path+processfile+"'");
  }
  m_path=dummy.Path();
  processfile=dummy.File();
  dummy.Close();
  ifstream from((m_path+processfile).c_str());
  if (!from) {
    THROW(critical_error,"Cannot open file '"+m_path+processfile+"'");
  }
  size_t position, order_ew, order_strong, psmc;
  int         flag;
  string      buf,ini,fin;
  int         nIS,   nFS;
  std::vector<Flavour>  IS, FS, flavs;
  std::string efunc="1", printgraphs;
  double fixed_scale;
  Data_Reader reader(" ",";","!","=");
  reader.AddWordSeparator("\t");
  reader.AddIgnore(":");
  while(from) {
    getline(from,buf);
    if (buf[0] != '%' && buf.length()>0) {
      fixed_scale=fix_scale;
      order_ew=99;
      order_strong=99;
      psmc=0;
      efunc="1";
      printgraphs="";
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
	  if (!p_isrhandler->CheckConsistency(&IS.front())) {
	    msg_Error()<<"Error in initialising ISR_Handler."<<endl
		       <<" "<<p_isrhandler->Flav(0)<<" -> "<<IS[0]<<", "
		       <<" "<<p_isrhandler->Flav(1)<<" -> "<<IS[1]<<endl
		       <<"  Delete it and ignore the process."<<endl;
	  }
	  else {
	    do {
	      getline(from,buf);
	      if (buf[0] != '%' && buf.length()>0) {
		reader.SetString(buf);
		unsigned int order_ew_t, order_strong_t;
		double fixed_scale_t;
		std::string efunc_t="1", printgraphs_t;
		if (reader.ReadFromString(order_ew_t,"electroweak")) order_ew=order_ew_t;
		if (reader.ReadFromString(order_strong_t,"strong")) order_strong=order_strong_t;
		if (reader.ReadFromString(fixed_scale_t,"scale")) fixed_scale=fixed_scale_t;
		if (reader.ReadFromString(efunc_t,"Enhance_Function")) efunc=efunc_t;
		reader.SetTags(Integrable_Base::ColorSchemeTags());
		reader.SetTags(std::map<std::string,std::string>());
		if (reader.ReadFromString(printgraphs_t,"Print_Graphs"))
		  printgraphs=printgraphs_t;
		position = buf.find(string("process"));
		if (!from) {
		  msg_Error()<<METHOD<<"("<<m_path+processfile<<")."<<endl
			     <<"   EOF without 'End process' tag."<<endl;
		  position     = 0;
		}
	      }
	    }
	    while (position==std::string::npos);
	    if (nIS+nFS>(int)m_nmax) m_nmax = nIS+nFS;
	    flavs.resize(nIS+nFS);
	    for (int i=0;i<nIS;i++) flavs[i]     = IS[i]; 
	    for (int i=0;i<nFS;i++) flavs[i+nIS] = FS[i]; 
	    double inisum=0.0, finsum=0.0;
	    size_t qcdjets(0);
	    for (int i=0;i<nIS;i++) inisum+=flavs[i].Mass();
	    for (int i=0;i<nFS;i++) {
	      finsum+=flavs[i+nIS].Mass();
	      if (flavs[i+nIS].Strong()) ++qcdjets;
	    }
	    m_minqcdjet=Min(m_minqcdjet,qcdjets);
	    m_maxqcdjet=ATOOLS::Max(m_maxqcdjet,qcdjets);
	    if (inisum<rpa.gen.Ecms() && finsum<rpa.gen.Ecms()) {
	      InitializeProcess(&flavs.front(),efunc,printgraphs,psmc,
				inisum,finsum,order_ew,order_strong,
				nIS,nFS,fixed_scale);
	      if (m_xsecs.size()>0) p_selected=m_xsecs.back();
	    }
	  }
	}
      }
    }
  }
  if (msg_LevelIsTracking()) this->Print();
  m_maxjet=m_nmax-m_nin;
  SetCoreMaxJetNumber(m_maxjet);
  if (m_xsecs.size()>0) return Tests();
  msg_Error()<<METHOD<<"(): No valid process in '"
	     <<m_path+processfile<<"'."<<std::endl;
  return false;
}

void Simple_XS::InitializeProcess(ATOOLS::Flavour *flavs,std::string &efunc,
 				  std::string &printgraphs,bool psmc,
  				  double &inisum,double &finsum,
  				  size_t order_ew,size_t order_strong,
 				  size_t nin,size_t nout,
				  double &fixscale)
{
  size_t nt(0);
  for (size_t j=0;j<nin+nout;++j) nt+=flavs[j].Size();
  XS_Base *newxs(NULL);
  if (nt>nin+nout) {
    newxs = new XS_Group(nin,nout,flavs,m_scalescheme,m_kfactorscheme,
 			 p_beamhandler,p_isrhandler,p_selectordata,p_model);
    newxs->SetFactorizationScale(m_muf2tag);
    newxs->SetRenormalizationScale(m_mur2tag);
    if (!((XS_Group*)newxs)->
	ConstructProcesses(order_ew,order_strong,fixscale)) {
      delete newxs;
      return;
    }
  }
  else {
    newxs = XSSelector()->
      GetXS(nin,nout,flavs,false,order_ew,order_strong);
    if (newxs==NULL) return;
    newxs->SetFactorizationScale(m_muf2tag);
    newxs->SetRenormalizationScale(m_mur2tag);
    newxs->SetScales(fixscale);
    newxs->Initialize(m_scalescheme,m_kfactorscheme,
 		      p_beamhandler,p_isrhandler,p_selectordata);
  }
  newxs->SetISRThreshold(ATOOLS::Max(inisum,finsum));
  newxs->SetEnhanceFunction(efunc);
  Add(newxs);
}
  
int Simple_XS::ExtractFlavours(std::vector<ATOOLS::Flavour> &flavours,
			       std::string buffer)
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
  flavours.resize(count);
  for (i=0;i<count;i++) {
    flavours[i] = Flavour((kf_code)(iabs(ii[i])));
    if (ii[i]<0) flavours[i] = flavours[i].Bar();
  }
  return count;
}

bool Simple_XS::CalculateTotalXSec(const std::string &resultpath) 
{
  Reset();
  bool okay = 1;
  for (size_t i=0;i<m_xsecs.size();++i) {
    okay = okay && m_xsecs[i]->CalculateTotalXSec(resultpath,false);
    p_activepshandler=m_xsecs[i]->PSHandler(false);
    m_totalxs += m_xsecs[i]->TotalXS();
  }
  msg_Info()<<"In Simple_XS::CalculateTotalXSec() = "
	    <<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb."<<endl;
  return okay;
}

bool  Simple_XS::SelectOne() 
{
  DeSelect();
  if (m_totalxs==0) 
    p_selected = m_xsecs[ATOOLS::Min(size_t(ran.Get()*m_xsecs.size()),
				     m_xsecs.size())];
  else {
    double disc = m_totalxs * ran.Get(); 
    for (size_t i=0;i<m_xsecs.size();++i) {
      disc -= m_xsecs[i]->TotalXS();
      if (disc<0.) {
	p_selected = m_xsecs[i];
	p_selected->DeSelect();
	return true;
      }
    }
    msg_Error()<<"Process_Group::SelectOne(): "
	       <<"Total xsec = "<<m_totalxs<<std::endl;
    p_selected = m_xsecs[0];
    p_selected->DeSelect();
    return false;
  }
  return true;
}

ATOOLS::Blob_Data_Base *Simple_XS::OneEvent(const int mode)
{
  SelectOne();
  return p_selected->OneEvent();
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


void Simple_XS::ResetSelector(ATOOLS::Selector_Data *const selectordata)
{
  for (unsigned int i=0;i<m_xsecs.size();++i) m_xsecs[i]->ResetSelector(selectordata);
}
