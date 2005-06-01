#include "XS_Group.H"

#include "ISR_Handler.H"
#include "Phase_Space_Handler.H"
#include "XS_Selector.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Random.H"
#include "MathTools.H"
#include "FSR_Channel.H"

#ifdef PROFILE__all
#define PROFILE__XS_Group
#endif
#ifdef PROFILE__XS_Group
#include "prof.hh"
#else
#define PROFILE_HERE
#endif

using namespace EXTRAXS;

XS_Group::XS_Group(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours,
		   const int scalescheme,const int kfactorscheme,
		   BEAM::Beam_Spectra_Handler *const beamhandler,
		   PDF::ISR_Handler *const isrhandler,
		   ATOOLS::Selector_Data *const selectordata):
  XS_Base(nin,nout,flavours,scalescheme,kfactorscheme,
	  beamhandler,isrhandler,selectordata),
  m_atoms(false), m_channels(false), p_xsselector(new XS_Selector(this)) 
{
  p_selected=NULL;
  p_selectordata=selectordata;
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
  if (xsec==NULL) return;
  if (m_xsecs.size()==0) {
    m_nin=xsec->NIn();
    m_nout=xsec->NOut();
    p_flavours=new ATOOLS::Flavour[m_nin+m_nout];
    for (size_t i=0;i<m_nin+m_nout;++i) p_flavours[i]=xsec->Flavours()[i];
    m_neweak=xsec->m_neweak;
    m_nstrong=xsec->m_nstrong;
    double inisum=0.0, finsum=0.0;
    for (size_t i=0;i<m_nin;i++) inisum+=p_flavours[i].Mass();
    for (size_t i=0;i<m_nout;i++) finsum+=p_flavours[i+m_nin].Mass();
    m_threshold=ATOOLS::Max(inisum,finsum);
  }
  if (m_neweak!=xsec->m_neweak) m_neweak=0;
  if (m_nstrong!=xsec->m_nstrong) m_nstrong=0;
  if (xsec->NVector()>m_nvector) CreateMomenta(xsec->NVector());
  if (xsec->NAddIn()>m_naddin || xsec->NAddOut()>m_naddout) {
    if (p_addmomenta!=NULL) delete [] p_addmomenta;
    m_naddin=xsec->NAddIn();
    m_naddout=xsec->NAddOut(); 
    p_addmomenta = new ATOOLS::Vec4D[m_naddin+m_naddout];
  }
  else {
    if (m_nin!=xsec->NIn() || m_nout!=xsec->NOut()) {
      ATOOLS::msg.Error()<<"XS_Group::Add("<<xsec<<"): ("<<this<<") Cannot add process '"
			 <<xsec->Name()<<"' to group '"<<m_name<<"' !"<<std::endl
			 <<"   Inconsistent number of external legs."<<std::endl; 
      return;
    }
  }  
  m_xsecs.push_back(xsec);
  if (!m_atoms) xsec->SetParent(this);
  for (size_t i=0;i<xsec->Resonances().size();++i) {
    bool present=false;
    for (size_t j=0;j<m_resonances.size();++j) {
      if (m_resonances[j]==xsec->Resonances()[i]) {
	present=true;
	break;
      }
    }
    if (!present) m_resonances.push_back(xsec->Resonances()[i]);
  }
  p_selected=m_xsecs[0];
  if (!m_atoms) xsec->SetPSHandler(p_activepshandler);
}

bool XS_Group::Remove(XS_Base *const xsec) 
{
  for (std::vector<XS_Base*>::iterator xsit=m_xsecs.begin();
       xsit!=m_xsecs.end();++xsit) {
    if (*xsit==xsec) {
      m_xsecs.erase(xsit);
      xsec->SetParent(xsec);
      return true;
    }
  }
  return false;
}

bool XS_Group::Delete(XS_Base *const xsec) 
{
  for (std::vector<XS_Base*>::iterator xsit=m_xsecs.begin();
       xsit!=m_xsecs.end();++xsit) {
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
    delete m_xsecs.back();
    m_xsecs.pop_back();
  }
}

bool XS_Group::SelectOne()
{
  DeSelect();
  if (m_totalxs==0.) p_selected=m_xsecs[int(ATOOLS::ran.Get()*m_xsecs.size())];
  else {
    double disc;
    if (m_atoms) disc=m_totalxs*ATOOLS::ran.Get();
    else disc=m_max*ATOOLS::ran.Get();
    for (size_t i=0;i<m_xsecs.size();++i) {
      if (m_atoms) disc-=m_xsecs[i]->TotalXS();
      else disc-=m_xsecs[i]->Max();
      if (disc<=0.) {
	p_selected=m_xsecs[i];
	p_selected->SetPSHandler(p_pshandler);
	p_selected->SelectOne();
	return true;
      }
    }
    if (disc>0.) { 
      ATOOLS::msg.Error()<<"XS_Group::SelectOne() : Cannot select process !"<<std::endl;
      if (m_atoms) ATOOLS::msg.Error()<<"   \\dsigma_{max} = "<<m_max<<std::endl;
      else ATOOLS::msg.Error()<<"   \\sigma_{tot} = "<<m_totalxs<<std::endl;
      return false;
    }
  }
  return true;
}

void XS_Group::WriteOutXSecs(std::ofstream &outfile)
{
  outfile.precision(12);
  outfile<<m_name<<"  "<<m_totalxs<<"  "<<m_max<<"  "<<m_totalerr<<" "
	 <<m_totalsum<<" "<<m_totalsumsqr<<" "<<m_n<<" "
	 <<m_ssum<<" "<<m_ssumsqr<<" "<<m_ssigma2<<" "<<m_sn<<" "<<m_wmin<<" "<<m_son<<std::endl; 
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

bool XS_Group::CalculateTotalXSec(const std::string &resultpath,
				  const bool create)
{
  if (m_atoms) {
    bool okay=true;
    for (size_t i=0;i<m_xsecs.size();++i) {
      if (!m_xsecs[i]->CalculateTotalXSec(resultpath,create)) okay=false;
    }
    return okay;
  }
  else {
    for (size_t i=0;i<m_xsecs.size();++i) {
      m_xsecs[i]->SetPSHandler(p_pshandler);
    }
    Reset();
    if (p_isrhandler) {
      if (m_nin==2) {
	if (p_flavours[0].Mass()!=p_isrhandler->Flav(0).Mass() ||
	    p_flavours[1].Mass()!=p_isrhandler->Flav(1).Mass()) {
	  p_isrhandler->SetPartonMasses(p_flavours);
	}
      }
      SetISR(p_isrhandler);
      p_pshandler->InitCuts();
      p_isrhandler->SetSprimeMin(p_pshandler->Cuts()->Smin());
    }
    CreateFSRChannels();
    if (!m_channels) {
      p_pshandler->CreateIntegrators();
      CreateISRChannels();
      m_channels = true;
    }
    std::string filename=resultpath+std::string("/")+m_name+std::string(".xs_tot"), singlename;
    double singlexs, singleerr, singlemax, singlesum, singlesumsqr,ssum,ssqrsum,ss2,wmin;
    long unsigned int singlen,sn,son;
    m_foundown=false;
    if (resultpath!=std::string("")) {
      std::ifstream infile;
      int hits=m_xsecs.size()+1;
      infile.open(filename.c_str());
      if (infile.good()) {
	infile>>singlename>>singlexs>>singlemax>>singleerr
	      >>singlesum>>singlesumsqr>>singlen
	      >>ssum>>ssqrsum>>ss2>>sn>>wmin>>son;
	do {
	  msg_Tracking()<<"Found result: xs for "<<singlename<<" : "
			<<singlexs*ATOOLS::rpa.Picobarn()<<" pb"
			<<" +/- "<<singleerr/singlexs*100.<<"%,"<<std::endl
			<<"         max : "<<singlemax<<std::endl;
	  XS_Base *xs=Matching(singlename);
	  if (xs!=NULL) {
	    xs->SetTotalXS(singlexs);
	    xs->SetTotalError(singleerr);
	    if (xs!=this) xs->SetMax(singlemax,true);
	    xs->SetSum(singlesum);
	    xs->SetSumSqr(singlesumsqr);
	    xs->SetPoints(singlen);
	    xs->SetSSum(ssum);
	    xs->SetSSumSqr(ssqrsum);
	    xs->SetSigmaSum(ss2);
	    xs->SetSPoints(sn);
	    xs->SetWMin(wmin);
	    xs->SetOptCounter(son);
	    --hits;
	  }
	  infile>>singlename>>singlexs>>singlemax>>singleerr
		>>singlesum>>singlesumsqr>>singlen
		>>ssum>>ssqrsum>>ss2>>sn>>wmin>>son;
	} while (infile);
      }
      infile.close();
      SetMax(0.,false);
      if (hits==0) {
	if (!p_pshandler->ReadIn(resultpath+std::string("/MC_")+m_name)) {
	  for (size_t i=0;i<m_xsecs.size();++i) m_xsecs[i]->Reset();
	  Reset();
	  m_foundown=false;
	}
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
    double var=TotalVar();
    m_totalxs=p_pshandler->Integrate()/ATOOLS::rpa.Picobarn(); 
    if (!(ATOOLS::IsZero((m_totalxs-TotalResult())/(m_totalxs+TotalResult())))) {
      ATOOLS::msg.Error()<<"Result of PS-Integrator and internal summation do not coincide!"<<std::endl
			 <<"  "<<m_name<<" : "<<m_totalxs<<" vs. "<<TotalResult()<<std::endl;
    }
    if (m_totalxs>0.) {
      SetTotal();
      if (var==TotalVar()) {
	ATOOLS::Exception_Handler::RemoveTerminatorObject(this);
	return 1;
      }
      if (resultpath!=std::string("")) {
	std::ofstream to;
	to.open(filename.c_str(),std::ios::out);
	to.precision(12);
	msg_Info()<<"Store result : xs for "<<m_name<<" : ";
	WriteOutXSecs(to);
	if (m_nin==2) msg_Info()<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb";
	if (m_nin==1) msg_Info()<<m_totalxs<<" GeV";
	msg_Info()<<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<std::endl
		  <<"       max : "<<m_max<<std::endl;
	p_pshandler->WriteOut(resultpath+std::string("/MC_")
			      +m_name,create);
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
  msg_Info()<<"Store result to "<<m_resultpath<<","<<std::endl
		    <<"   xs for "<<m_name<<" : ";
  WriteOutXSecs(to);
  if (m_nin==2) msg_Info()<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb";
  if (m_nin==1) msg_Info()<<m_totalxs<<" GeV";
  msg_Info()<<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<std::endl
		      <<"       max : "<<m_max<<std::endl;
  p_pshandler->WriteOut(m_resultpath+std::string("/MC_")+m_name);
  to.close();
}

void XS_Group::SetTotal()  
{ 
  m_totalxs=TotalResult();
  m_totalerr=TotalVar();
  if (p_selector) p_selector->Output();
  m_max=0.;
  for (size_t i=0;i<m_xsecs.size();++i) {
    m_xsecs[i]->SetTotal();
    m_max+=m_xsecs[i]->Max();
  }
  msg_Info()<<"Total XS for "<<ATOOLS::om::bold<<m_name<<" : "
	    <<ATOOLS::om::blue<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb"
	    <<ATOOLS::om::reset<<" +/- ( "<<ATOOLS::om::red
	    <<m_totalerr<<" pb = "<<m_totalerr/m_totalxs*100.
	    <<" %"<<ATOOLS::om::reset<<" )"<<std::endl
	    <<"      max = "<<m_max<<"\n"<<ATOOLS::om::bold<<m_name
	    <<ATOOLS::om::reset<<" : "<<ATOOLS::om::blue<<ATOOLS::om::bold
	    <<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb"<<ATOOLS::om::reset
	    <<" +/- "<<ATOOLS::om::reset<<ATOOLS::om::blue<<m_totalerr/m_totalxs*100.
	    <<" %,"<<ATOOLS::om::reset<<ATOOLS::om::bold<<" exp. eff: "
	    <<ATOOLS::om::red<<(100.*m_totalxs/m_max)<<" %."
	    <<ATOOLS::om::reset<<std::endl;
}

bool XS_Group::OneEvent() 
{
  if (m_atoms) {
    SelectOne();
    return p_selected->OneEvent();
  }
  return p_activepshandler->OneEvent();
}

ATOOLS::Blob_Data_Base *XS_Group::WeightedEvent(const int mode) 
{
  if (m_atoms) {
    SelectOne();
    return p_selected->WeightedEvent(mode);
  }
  return p_activepshandler->WeightedEvent(mode);
}

void XS_Group::AddPoint(const double value) 
{
  Integrable_Base::AddPoint(value);
  for (size_t i=0;i<m_xsecs.size();++i) {
    if (ATOOLS::dabs(m_last)>0.) {
      m_xsecs[i]->AddPoint(value*m_xsecs[i]->Last()/m_last);
    }
    else {
      m_xsecs[i]->AddPoint(0.);
    }  
  }
}

void XS_Group::ResetMax(int flag) 
{
  m_max = 0.;
  for (size_t i=0;i<m_xsecs.size();i++) {
    m_xsecs[i]->ResetMax(flag);
  }
}

void XS_Group::OptimizeResult() 
{
  double ssigma2 = ATOOLS::sqr(m_ssum/m_sn)/((m_ssumsqr/m_sn - ATOOLS::sqr(m_ssum/m_sn))/(m_sn-1));
  if (ssigma2>m_wmin) {
    m_ssigma2  += ssigma2; 
    m_totalsum += m_ssum*ssigma2/m_sn;
    m_totalsumsqr+= m_ssumsqr*ssigma2/m_sn;
    m_ssum     = 0.;
    m_ssumsqr  = 0.;
    m_sn       = 0;
    if (ssigma2/m_son>m_wmin) m_wmin = ssigma2/m_son;
    m_son      = 0;
  }
  m_son++;
  for (size_t i=0;i<m_xsecs.size();i++) m_xsecs[i]->OptimizeResult();
}

double XS_Group::Differential(double s,double t,double u)
{
  m_last = 0;
  for (size_t i=0;i<m_xsecs.size();++i) {
    m_xsecs[i]->SetMomenta(p_momenta);
    if (m_naddin>0 || m_naddout>0) 
      m_xsecs[i]->SetAddMomenta(p_addmomenta);
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
      ATOOLS::msg.Error().precision(12);
      ATOOLS::msg.Error()<<"XS_Group::SetMax(..): "
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
  if (m_nout>2) {
    p_pshandler->FSRIntegrator()->
      Add(new PHASIC::T1Channel(m_nin,m_nout,p_flavours));
  }
  else {
    p_pshandler->FSRIntegrator()->
      Add(new PHASIC::S1Channel(m_nin,m_nout,p_flavours));
    p_pshandler->FSRIntegrator()->
      Add(new PHASIC::T1Channel(m_nin,m_nout,p_flavours));
    p_pshandler->FSRIntegrator()->
      Add(new PHASIC::U1Channel(m_nin,m_nout,p_flavours));
  }
  for (size_t i=0;i<m_resonances.size();++i) {
    p_pshandler->FSRIntegrator()->
      Add(new PHASIC::S1Channel(m_nin,m_nout,p_flavours,m_resonances[i]));
  }
}      

void XS_Group::CreateISRChannels() 
{
}

void XS_Group::DeSelect() 
{
  for (size_t i=0;i<m_xsecs.size();++i) m_xsecs[i]->DeSelect();
  p_selected=NULL;
}

bool XS_Group::ReSelect(int) {
  return SelectOne();
}

void XS_Group::SetISR(PDF::ISR_Handler *const isrhandler) 
{
  for (size_t i=0;i<m_xsecs.size();++i) m_xsecs[i]->SetISR(isrhandler);
  p_isrhandler=isrhandler;
}

void XS_Group::SetAtoms(bool _atoms) 
{ 
  m_atoms = _atoms; 
}

XS_Base *const XS_Group::operator[](const size_t i) const
{ 
  if (i<m_xsecs.size()) return m_xsecs[i];
  return NULL;
} 

size_t XS_Group::Size() const
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

void XS_Group::SetPSHandler(PHASIC::Phase_Space_Handler *const pshandler) 
{
  for (size_t i=0;i<m_xsecs.size();++i) m_xsecs[i]->SetPSHandler(pshandler);
  p_activepshandler=pshandler;
} 

void XS_Group::ResetSelector(ATOOLS::Selector_Data *const selectordata)
{
  PROFILE_HERE;
  for (unsigned int i=0;i<m_xsecs.size();++i) {
    m_xsecs[i]->ResetSelector(selectordata);
  }
  XS_Base::ResetSelector(selectordata);
}

void XS_Group::Print()
{
  ATOOLS::msg.Out()<<(m_name!=""?m_name:"<no name>")<<" {\n";
  {
    msg_Indent();
    for (size_t i=0;i<m_xsecs.size();++i) m_xsecs[i]->Print();
  }
  ATOOLS::msg.Out()<<"} "<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb\n";
}

void XS_Group::SetEvents(const double number) 
{
  m_expevents=m_dicedevents=0;
  m_anasum=m_validanasum=0.0;
  for (size_t i=0;i<m_xsecs.size();++i) {
    m_xsecs[i]->SetEvents(number);
    m_expevents+=m_xsecs[i]->ExpectedEvents();
  }
}

bool XS_Group::SelectOneFromList()
{
  std::cout.precision(12);
  DeSelect();
  //  if (ATOOLS::IsEqual(m_xssum,m_dicedxssum)) return false;
  for (size_t i=0;i<m_xsecs.size();i++) {
    if (m_xsecs[i]->DicedEvents()<m_xsecs[i]->ExpectedEvents()) {
//       PRINT_INFO(m_xsecs[i]->DicedEvents()<<" "<<m_xsecs[i]->Events()<<" "<<m_xsecs[i]->Name());
      p_selected=m_xsecs[i];
      if (p_selected->SelectOneFromList()) break;
    }
  }
//    PRINT_INFO(Name()<<" selected "<<(p_selected==NULL?"NULL":p_selected->Name()));
  if (p_selected==NULL) return false;
  return true;
}

void XS_Group::AddEvent(const double xs,const double validxs,const int ncounts)
{ 
  m_dicedevents+=ncounts; 
  p_selected->AddEvent(xs,validxs,ncounts);
}

void XS_Group::ResetEvents()
{
}

void XS_Group::SetISRThreshold(const double threshold) 
{
  m_threshold=threshold;
  for (size_t i=0;i<m_xsecs.size();i++) 
    m_xsecs[i]->SetISRThreshold(threshold);
}
