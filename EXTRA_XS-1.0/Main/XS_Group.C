#include "XS_Group.H"

#include "Phase_Space_Handler.H"
#include "XS_Selector.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Random.H"
#include "MathTools.H"
#include "FSR_Channel.H"

using namespace EXTRAXS;

XS_Group::XS_Group(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours,
		   const int scalescheme,const int kfactorscheme,const double scalefactor,
		   BEAM::Beam_Spectra_Handler *const beamhandler,
		   PDF::ISR_Handler *const isrhandler,
		   ATOOLS::Selector_Data *const selectordata):
  XS_Base(nin,nout,flavours,scalescheme,kfactorscheme,scalefactor,
	  beamhandler,isrhandler,selectordata),
  m_atoms(false), m_channels(false), p_xsselector(new XS_Selector(this)) 
{
  p_selected=NULL;
}

XS_Group::XS_Group(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours):
  XS_Base(nin,nout,flavours), 
  m_atoms(false), m_channels(false), p_xsselector(new XS_Selector(this)) 
{
  p_selected=NULL;
}

XS_Group::XS_Group(const size_t nin,const size_t nout,const std::string &name):
  XS_Base(nin,nout,NULL), 
  m_atoms(false), m_channels(false), p_xsselector(new XS_Selector(this)) 
{
  p_selected=NULL;
  m_name=name;
}

XS_Group::XS_Group(): 
  m_atoms(false), m_channels(false), p_xsselector(new XS_Selector(this)) 
{
  p_selected=NULL;
}

XS_Group::~XS_Group()
{
  delete p_xsselector;
  Clear();
}

void XS_Group::Add(XS_Base *const xsec) 
{
  if (m_xsecs.size()==0) {
    m_nin=xsec->NIn();
    m_nout=xsec->NOut();
    p_flavours=new ATOOLS::Flavour[m_nin+m_nout];
    for (size_t i=0;i<m_nin+m_nout;++i) p_flavours[i]=xsec->Flavours()[i];
  }
  else {
    if (m_nin!=xsec->NIn() || m_nout!=xsec->NOut()) {
      ATOOLS::msg.Error()<<"XS_Group::Add("<<xsec<<"): ("<<this<<") Cannot add process '"
			 <<xsec->Name()<<"' to group '"<<m_name<<"' !"<<std::endl
			 <<"   Inconsistent number of external legs."<<std::endl; 
      return;
    }
  }  
  ATOOLS::msg.Tracking()<<"XS_Group::Add("<<xsec<<"): ("<<this<<") "
			<<"Adding '"<<xsec->Name()<<"' to group '"<<m_name<<"'."<<std::endl; 
  m_xsecs.push_back(xsec);
  m_nvector=ATOOLS::Max(m_nvector,xsec->NVector());
  p_selected=m_xsecs[0];
}

bool XS_Group::Delete(XS_Base *const xsec) 
{
  for (std::vector<XS_Base*>::iterator xsit=m_xsecs.begin();xsit!=m_xsecs.end();++xsit) {
    if (*xsit==xsec) {
      delete *xsit;
      m_xsecs.erase(xsit);
      return true;
    }
  }
  return false;
}

void XS_Group::Clear() 
{
  while (m_xsecs.size()>0) {
    delete *m_xsecs.begin();
    m_xsecs.erase(m_xsecs.begin());
  }
}

void XS_Group::SelectOne()
{
  DeSelect();
  if (m_totalxs==0) p_selected=m_xsecs[int(ATOOLS::ran.Get()*m_xsecs.size())];
  else {
    double disc;
    if (m_atoms) disc=m_max*ATOOLS::ran.Get();
    else disc=m_totalxs*ATOOLS::ran.Get();
    for (size_t i=0;i<m_xsecs.size();++i) {
      if (m_atoms) disc-=m_xsecs[i]->Max();
      else disc-=m_xsecs[i]->TotalXS();
      if (disc<=0.) {
	p_selected=m_xsecs[i];
	p_selected->SelectOne();
	return;
      }
    }
    if (disc>0.) { 
      ATOOLS::msg.Error()<<"Process_Group::SelectOne() : Cannot select process !"<<std::endl;
      if (m_atoms) ATOOLS::msg.Error()<<"   \\dsigma_{max} = "<<m_max<<std::endl;
      else ATOOLS::msg.Error()<<"   \\sigma_{tot} = "<<m_totalxs<<std::endl;
    }
  }
}

void XS_Group::WriteOutXSecs(std::ofstream &outfile)
{
  outfile<<m_name<<"  "<<m_totalxs<<"  "<<m_max<<"  "<<m_totalerr<<" "
	 <<m_totalsum<<" "<<m_totalsumsqr<<" "<<m_n<<std::endl;
  for (size_t i=0;i<m_xsecs.size();++i) m_xsecs[i]->WriteOutXSecs(outfile);
}

XS_Base *XS_Group::Matching(const std::string &name)
{
  if (name==m_name && !m_foundown) {
    m_foundown=true;
    return this;
  }
  else {
    for (size_t i=0;i<m_xsecs.size();++i) {
      if (name==m_xsecs[i]->Name()) return m_xsecs[i];
    }
  }
  return NULL;
}

bool XS_Group::CalculateTotalXSec(const std::string &resultpath)
{
  if (m_atoms) {
    bool okay=true;
    for (size_t i=0;i<m_xsecs.size();++i) {
      if (!m_xsecs[i]->CalculateTotalXSec(resultpath)) okay=false;
    }
    return okay;
  }
  else {
    Reset();
    if (p_isrhandler) {
      if (m_nin==2) {
	if (p_flavours[0].Mass()!=p_isrhandler->Flav(0).Mass() ||
	    p_flavours[1].Mass()!=p_isrhandler->Flav(1).Mass()) {
	  p_isrhandler->SetPartonMasses(p_flavours);
	}
      }
      SetISR(p_isrhandler);
    }
    CreateFSRChannels();
    if (!m_channels) {
      p_pshandler->CreateIntegrators();
      m_channels = true;
    }
    std::string filename=resultpath+std::string("/")+m_name+std::string(".xstotal"), singlename;
    double singlexs, singleerr, singlemax, singlesum, singlesumsqr;
    long unsigned int singlen;
    if (resultpath!=std::string("")) {
      m_foundown=false;
      std::ifstream infile;
      int hits=m_xsecs.size()+1;
      infile.open(filename.c_str());
      if (infile.good()) {
	infile>>singlename>>singlexs>>singlemax>>singleerr>>singlesum>>singlesumsqr>>singlen;
	do {
	  ATOOLS::msg.Events()<<"Found result: xs for "<<singlename<<" : "
			      <<singlexs*ATOOLS::rpa.Picobarn()<<" pb"
			      <<" +/- "<<singleerr/singlexs*100.<<"%,"<<std::endl
			      <<"         max : "<<singlemax<<std::endl;
	  XS_Base *xs=Matching(singlename);
	  if (xs!=NULL) {
	    xs->SetTotalXS(singlexs);
	    xs->SetTotalError(singleerr);
	    xs->SetMax(singlemax,false);
	    xs->SetSum(singlesum);
	    xs->SetSumSqr(singlesumsqr);
	    xs->SetPoints(singlen);
	    --hits;
	  }
	  infile>>singlename>>singlexs>>singlemax>>singleerr>>singlesum>>singlesumsqr>>singlen;
	} while (infile);
      }
      infile.close();
      if (hits==0) {
	p_pshandler->ReadIn(resultpath+std::string("/MC_")+m_name);
	if (p_pshandler->BeamIntegrator() != 0) p_pshandler->BeamIntegrator()->Print();
	if (p_pshandler->ISRIntegrator() != 0) p_pshandler->ISRIntegrator()->Print();
	if (p_pshandler->KMRZIntegrator() != 0) p_pshandler->KMRZIntegrator()->Print();
	if (p_pshandler->KMRKPIntegrator() != 0) p_pshandler->KMRKPIntegrator()->Print();
	if (p_pshandler->FSRIntegrator() != 0) p_pshandler->FSRIntegrator()->Print();
      }
    }
    if (m_foundown) SetTotal();
    m_resultpath=resultpath;
    m_resultfile=filename;
    ATOOLS::Exception_Handler::AddTerminatorObject(this);
    p_pshandler->InitIncoming();
    m_totalxs=p_pshandler->Integrate()/ATOOLS::rpa.Picobarn(); 
    if (!(ATOOLS::IsZero((m_n*m_totalxs-m_totalsum)/(m_n*m_totalxs+m_totalsum)))) {
      ATOOLS::msg.Error()<<"Result of PS-Integrator and internal summation do not coincide!"<<std::endl
			 <<"  "<<m_name<<" : "<<m_totalxs<<" vs. "<<m_totalsum/m_n<<std::endl;
    }
    if (m_totalxs>0.) {
      SetTotal();
      if (resultpath!=std::string("")) {
	std::ofstream to;
	to.open(filename.c_str(),std::ios::out);
	to.precision(12);
	ATOOLS::msg.Events()<<"Store result : xs for "<<m_name<<" : ";
	WriteOutXSecs(to);
	if (m_nin==2) ATOOLS::msg.Events()<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb";
	if (m_nin==1) ATOOLS::msg.Events()<<m_totalxs<<" GeV";
	ATOOLS::msg.Events()<<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<std::endl
			    <<"       max : "<<m_max<<std::endl;
	p_pshandler->WriteOut(resultpath+std::string("/MC_")+m_name);
	to.close();
      }
      ATOOLS::Exception_Handler::RemoveTerminatorObject(this);
      return 1;
    }
    ATOOLS::Exception_Handler::RemoveTerminatorObject(this);
    return 0;
  }
}

void XS_Group::PrepareTerminate()  
{
  if (m_resultpath.length()==0 && m_resultfile.length()==0) return;
  SetTotal();
  if (m_totalxs<=0.) return;
  std::ofstream to;
  to.open(m_resultfile.c_str(),std::ios::out);
  to.precision(12);
  ATOOLS::msg.Events()<<"Store result to "<<m_resultpath<<std::endl;
  ATOOLS::msg.Events()<<"Store result : xs for "<<m_name<<" : ";
  WriteOutXSecs(to);
  if (m_nin==2) ATOOLS::msg.Events()<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb";
  if (m_nin==1) ATOOLS::msg.Events()<<m_totalxs<<" GeV";
  ATOOLS::msg.Events()<<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<std::endl
		      <<"       max : "<<m_max<<std::endl;
  p_pshandler->WriteOut(m_resultpath+std::string("/MC_")+m_name);
  to.close();
}

void XS_Group::SetTotal()  
{ 
  m_totalxs=m_totalsum/m_n; 
  m_totalerr=sqrt((m_n*m_totalsumsqr-ATOOLS::sqr(m_totalsum))/(m_n-1))/m_n;
  if (p_selector) p_selector->Output();
  m_max=0.;
  ATOOLS::msg.Events()<<ATOOLS::om::bold<<"--------------------------------------------------"
		      <<ATOOLS::om::reset<<std::endl;
  for (size_t i=0;i<m_xsecs.size();++i) {
    m_xsecs[i]->SetTotal();
    m_max+=m_xsecs[i]->Max();
  }
  ATOOLS::msg.Events()<<ATOOLS::om::bold<<"--------------------------------------------------"
		      <<ATOOLS::om::reset<<std::endl;
  ATOOLS::msg.Events()<<"Total XS for "<<ATOOLS::om::bold<<m_name<<" : "
		      <<ATOOLS::om::blue<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb"
		      <<ATOOLS::om::reset<<" +/- ( "<<ATOOLS::om::red<<m_totalerr<<" pb = "
		      <<m_totalerr/m_totalxs*100.<<" %"<<ATOOLS::om::reset<<" )"<<std::endl
		      <<"      max = "<<m_max<<std::endl;
  ATOOLS::msg.Events()<<ATOOLS::om::bold<<"--------------------------------------------------"
		      <<ATOOLS::om::reset<<std::endl;
}

bool XS_Group::OneEvent() 
{
  if (m_atoms) {
    SelectOne();
    return p_selected->OneEvent();
  }
  return p_pshandler->OneEvent();
}

ATOOLS::Blob_Data_Base *XS_Group::WeightedEvent() 
{
  if (m_atoms) {
    SelectOne();
     return p_selected->WeightedEvent();
  }
  return p_pshandler->WeightedEvent();
}

void XS_Group::AddPoint(const double value) 
{
  m_n++;
  m_totalsum    += value;
  m_totalsumsqr += value*value;
  if (value>m_max) m_max = value;

  for (size_t i=0;i<m_xsecs.size();++i) {
    if (ATOOLS::dabs(m_last)>0.) {
      m_xsecs[i]->AddPoint(value*m_xsecs[i]->Last()/m_last);
    }
    else {
      m_xsecs[i]->AddPoint(0.);
    }  
  }
}

double XS_Group::Differential(double s,double t,double u)
{
  m_last = 0;
  for (size_t i=0;i<m_xsecs.size();++i) {
    m_xsecs[i]->SetMomenta(p_momenta);
    m_last+=m_xsecs[i]->Differential(s,t,u);
  }
  if (!(m_last<=0) && !(m_last>0)) {
    ATOOLS::msg.Error()<<"XS_Group::Differential("<<s<<","<<t<<","<<u<<"): "<<ATOOLS::om::red
		       <<"Cross section is 'nan'!"<<ATOOLS::om::reset<<std::endl;
  }
  return m_last;
}

double XS_Group::Differential2()
{
  if (p_isrhandler) {
    if (p_isrhandler->On()==0) return 0.;
    double tmp = 0.;
    for (size_t i=0;i<m_xsecs.size();++i) tmp += m_xsecs[i]->Differential2();

    if ((!(tmp<=0)) && (!(tmp>0))) {
      ATOOLS::msg.Error()<<"---- X_Group::Differential -------------------"<<std::endl;
    }
    m_last += tmp;
    return tmp;
  }
  return 0.;
}

void XS_Group::SetMax(const double max,const int flag) 
{
  if (flag==1) {
    for (size_t i=0;i<m_xsecs.size();++i) m_xsecs[i]->SetMax(max/(double)Size(),flag);
    m_max=max;
    return;
  }
  double sum=0.;
  m_max=0.;
  for (size_t i=0;i<m_xsecs.size();++i) {
    sum+=m_xsecs[i]->TotalXS();
    m_max+=m_xsecs[i]->Max();
  }
  if (m_totalxs!=0.) {
    if (!ATOOLS::IsEqual(sum,m_totalxs)) {
      ATOOLS::msg.Events().precision(12);
      ATOOLS::msg.Events()<<"XS_Group::SetMax(..): "
			  <<"'"<<m_name<<"' : Summation does not agree !"<<std::endl
			  <<"   sum = "<<sum<<" vs. total = "<<m_totalxs
			  <<" ("<<((sum-m_totalxs)/m_totalxs)<<")"<<std::endl;
    }
    m_totalxs=sum;
  }
}

void XS_Group::CreateFSRChannels() 
{
  p_pshandler->FSRIntegrator()->DropAllChannels();
  p_pshandler->FSRIntegrator()->Add(new PHASIC::S1Channel(2,2,p_flavours,
							  ATOOLS::Flavour(ATOOLS::kf::photon)));
  p_pshandler->FSRIntegrator()->Add(new PHASIC::T1Channel(2,2,p_flavours));
  p_pshandler->FSRIntegrator()->Add(new PHASIC::U1Channel(2,2,p_flavours));
}      

void XS_Group::DeSelect() 
{
  for (size_t i=0;i<m_xsecs.size();++i) m_xsecs[i]->DeSelect();
  p_selected=NULL;
}

void XS_Group::SetISR(PDF::ISR_Handler *const isrhandler) 
{
  for (size_t i=0;i<m_xsecs.size();++i) m_xsecs[i]->SetISR(isrhandler);
  p_isrhandler=isrhandler;
}

XS_Base *const XS_Group::operator[](const size_t i) const
{ 
  if (i<m_xsecs.size()) return m_xsecs[i];
  return NULL;
} 

void XS_Group::SetAtoms(bool _atoms) 
{ 
  m_atoms = _atoms; 
}

size_t XS_Group::Size() 
{ 
  return m_xsecs.size(); 
}

bool XS_Group::SetColours(const double s,const double t,const double u)
{
  return false;
}

double XS_Group::operator()(const double s,const double t,const double u)
{
  return 0.;
}

void XS_Group::Reset()
{
  for (size_t i=0;i<m_xsecs.size();++i) m_xsecs[i]->Reset();
  XS_Base::Reset();
}
