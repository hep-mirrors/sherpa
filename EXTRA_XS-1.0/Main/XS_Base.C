#include "XS_Base.H"

#include "Phase_Space_Handler.H"
#include "ISR_Handler.H"
#include "Running_AlphaS.H"
#include "Combined_Selector.H"
#include "Standard_Selector.H"
#include "Run_Parameter.H"
#include "Regulator_Base.H"
#include "Message.H"
#include "MyStrStream.H"
#include "FSR_Channel.H"
#include "Color_Integrator.H"
#include "Sample_Multi_Channel.H"
#include "PreSample_Multi_Channel.H"
#include "VHAAG.H"
#include "Data_Reader.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;

XS_Base::XS_Base():
  p_colours(NULL), m_channels(false), m_psmc(false)
{
  m_name="Empty XS";
  m_scale[PHASIC::stp::sfs]=m_scale[PHASIC::stp::sis]=std::numeric_limits<double>::max();
} 

XS_Base::XS_Base(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours,
		 const PHASIC::scl::scheme scalescheme,const int kfactorscheme,
		 BEAM::Beam_Spectra_Handler *const beamhandler,
		 PDF::ISR_Handler *const isrhandler,
		 ATOOLS::Selector_Data *const selectordata):
  Integrable_Base(nin,nout,scalescheme,kfactorscheme,
		  beamhandler,isrhandler,selectordata),
  p_colours(NULL), m_channels(false), m_psmc(false)
{
  Init(flavours);
  Initialize(scalescheme,kfactorscheme,beamhandler,isrhandler,selectordata);
}

XS_Base::XS_Base(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours):
  Integrable_Base(nin,nout),
  p_colours(NULL), m_channels(false), m_psmc(false)
{
  Init(flavours);
  p_selector = new ATOOLS::No_Selector();
  m_scale[PHASIC::stp::sfs]=m_scale[PHASIC::stp::sis]=std::numeric_limits<double>::max();
}

XS_Base::~XS_Base() 
{
  if (p_colours!=NULL) { 
    for (size_t i=0;i<m_nin+m_nout;++i) delete [] p_colours[i];
    delete [] p_colours;
  }
}

class Order_KF {
public:
  bool operator()(const ATOOLS::Flavour &a,const ATOOLS::Flavour &b)
  { return a.Kfcode()<b.Kfcode(); }
};// end of class Order_KF

class Order_Anti {
public:
  bool operator()(const ATOOLS::Flavour &a,const ATOOLS::Flavour &b)
  { return a.IsFermion() && b.IsFermion()
      && (!a.IsAnti() && b.IsAnti()); }
};// end of class Order_Anti

class Order_SVFT {
public:
  bool operator()(const ATOOLS::Flavour &a,const ATOOLS::Flavour &b) 
  {
    if (a.IsScalar() && !b.IsScalar()) return true;
    if (a.IsVector() && !b.IsScalar() && !b.IsVector()) return true;
    if (a.IsFermion() && !b.IsFermion() && 
	!b.IsScalar() && !b.IsVector()) return true;
    return false;
  }
};// end of class Order_SVFT

class Order_Mass {
public:
  int operator()(const ATOOLS::Flavour &a,const ATOOLS::Flavour &b) 
  { return a.Mass()>b.Mass(); }
};// end of class Order_Mass

class Order_InvMass {
public:
  int operator()(const ATOOLS::Flavour &a,const ATOOLS::Flavour &b) 
  { return a.Mass()<b.Mass(); }
};// end of class Order_InvMass

class Order_Coupling {
public:
  int operator()(const ATOOLS::Flavour &a,const ATOOLS::Flavour &b) 
  { return !a.Strong() && b.Strong(); }
};// end of class Order_Coupling

void XS_Base::SortFlavours(ATOOLS::Flavour *flavs,const size_t &n)
{
  std::vector<ATOOLS::Flavour> fl(n);
  ATOOLS::Flavour heaviest(ATOOLS::kf::photon);
  for (size_t i(0);i<n;++i) {
    fl[i]=flavs[i];
    if (flavs[i].Mass()>heaviest.Mass()) heaviest=flavs[i];
    else if (flavs[i].Mass()==heaviest.Mass() &&
	     !flavs[i].IsAnti()) heaviest=flavs[i];
  }
  std::stable_sort(fl.begin(),fl.end(),Order_KF());
  std::stable_sort(fl.begin(),fl.end(),Order_Anti());
  std::stable_sort(fl.begin(),fl.end(),Order_SVFT());
  if (heaviest.IsAnti())  
    std::stable_sort(fl.begin(),fl.end(),Order_InvMass());
  else std::stable_sort(fl.begin(),fl.end(),Order_Mass());
  std::stable_sort(fl.begin(),fl.end(),Order_Coupling());
  for (size_t i(0);i<n;++i) flavs[i]=fl[i];
}

void XS_Base::Init(const ATOOLS::Flavour *flavours)
{
  if (m_nin+m_nout==0) return;
  p_flavours = new ATOOLS::Flavour[m_nin+m_nout];
  if (flavours!=NULL) {
    for (size_t i=0;i<m_nin+m_nout;++i) p_flavours[i]=flavours[i];
    SortFlavours(p_flavours,m_nin);
    SortFlavours(&p_flavours[m_nin],m_nout);
    m_name=GenerateName(m_nin,m_nout,p_flavours);
  }
  p_colours = new int*[m_nin+m_nout];
  for (size_t i=0;i<m_nin+m_nout;++i) { 
    p_colours[i] = new int[2]; 
    p_colours[i][0]=p_colours[i][1]=0; 
  }  
  double massin=0.0, massout=0.0;
  for (size_t i=0;i<m_nin;++i) massin+=p_flavours[i].Mass();
  for (size_t i=m_nin;i<m_nout;++i) massout+=p_flavours[i].Mass();
  if (massin>massout) m_threshold=ATOOLS::sqr(massin);
  else m_threshold=ATOOLS::sqr(massout);
}

std::string XS_Base::GenerateName(const size_t nin,const size_t nout,
				  const ATOOLS::Flavour *flavours) 
{
  std::string name(ATOOLS::ToString(nin)+"_"+ATOOLS::ToString(nout));
  for (size_t i(0);i<nin;++i) {
    name+="_"+std::string(flavours[i].IDName());
    if (flavours[i].Kfcode()==kf::quark && flavours[i].IsAnti())
      name+="b";
  }
  name+="__";
  for (size_t i(nin);i<nin+nout;++i) {
    name+="_"+std::string(flavours[i].IDName());
    if (flavours[i].Kfcode()==kf::quark && flavours[i].IsAnti())
      name+="b";
  }
  return name;
}

void XS_Base::Initialize(const PHASIC::scl::scheme scalescheme,
			 const int kfactorscheme,
			 BEAM::Beam_Spectra_Handler *const beamhandler,
			 PDF::ISR_Handler *const isrhandler,
			 ATOOLS::Selector_Data *const selectordata)
{
  p_beamhandler=beamhandler;
  p_isrhandler=isrhandler;
  ResetSelector(selectordata);
  SetScaleScheme(scalescheme);
  SetKFactorScheme(kfactorscheme);
  p_pshandler = new PHASIC::Phase_Space_Handler(this,isrhandler,beamhandler);
  p_activepshandler=p_pshandler;
  m_scale[PHASIC::stp::sfs]=m_scale[PHASIC::stp::sis]=
    std::numeric_limits<double>::max();
}

double XS_Base::Differential(const ATOOLS::Vec4D *momenta) 
{
  SetMomenta(momenta);
  SetSTU(momenta);
  return Differential(m_s,m_t,m_u);
}

bool XS_Base::SetColours(const ATOOLS::Vec4D *momenta) 
{
  SetMomenta(momenta);
  SetSTU(momenta);
  SetScale(momenta[2].PPerp2(),PHASIC::stp::ren);  
  return SetColours(m_s,m_t,m_u);
}

void XS_Base::SwapInOrder() 
{
  std::swap(p_flavours[0],p_flavours[1]);
  std::swap(p_momenta[0],p_momenta[1]);
  std::swap(p_colours[0],p_colours[1]);
  if (m_naddout>0) {
    std::swap(p_addmomenta[0],p_addmomenta[1]);
  }
  m_swaped=true;
}

void XS_Base::RestoreInOrder() 
{
  if (m_swaped) {
    std::swap(p_flavours[0],p_flavours[1]);
    std::swap(p_momenta[0],p_momenta[1]);
    std::swap(p_colours[0],p_colours[1]);
    if (m_naddout>0) {
      std::swap(p_addmomenta[0],p_addmomenta[1]);
    }
    m_swaped=false;
  }
}

void XS_Base::ResetSelector(ATOOLS::Selector_Data *const selectordata)
{
  if (p_selector!=NULL) delete p_selector;
  if (selectordata!=NULL) {
    p_selector = new ATOOLS::Combined_Selector(m_nin,m_nout,p_flavours,selectordata);
  }
  else {
    ATOOLS::msg.Error()<<"XS_Base::ResetSelector("<<selectordata<<"): "
		       <<"(\""<<m_name<<"\")"<<std::endl
		       <<"   No cuts specified. Initialize 'No_Selector'."<<std::endl;
    p_selector = new ATOOLS::No_Selector();
  }
  p_selector->SetProcessName(Name());
}

void XS_Base::SetSTU(const ATOOLS::Vec4D *momenta)
{
  m_s=(momenta[0]+momenta[1]).Abs2();
  m_t=(momenta[0]-momenta[2]).Abs2();
  m_u=(momenta[0]-momenta[3]).Abs2();
}

void XS_Base::SetMax(const double max,const int flag)             
{ 
  if (flag==1) m_max=max;
  SetMax();
}

void XS_Base::SetMax()             
{ 
}

XS_Base *const XS_Base::operator[](const size_t i) const 
{
  return NULL;
}

size_t XS_Base::Size() const
{ 
  return 0; 
}

bool XS_Base::SelectOne()
{ 
  return true;
}

bool XS_Base::ReSelect(int)
{ 
  return true;
}

bool XS_Base::SelectOneFromList()
{ 
  return true;
}

void XS_Base::DeSelect()
{ 
}

void XS_Base::Reset()
{
  m_n=0;
  m_last=m_lastlumi=m_lastdxs=0.0;
  m_totalxs=m_totalsum=m_totalsumsqr=m_totalerr=0.0;
  m_smax=m_max=m_wmin=m_ssigma2=m_ssumsqr=m_ssum=0.0;
  m_sn=0;
  m_son=1;
  m_vsmax.clear(); 
  m_vsn.clear();   
}

ATOOLS::Blob_Data_Base *XS_Base::SameWeightedEvent()
{
  return p_activepshandler->SameWeightedEvent();
}

void XS_Base::AssignRegulator(const std::string &regulator,
			      const std::vector<double> &parameters)
{
  PHASIC::Regulator_Base *function=NULL;
  if ((function=PHASIC::Regulator_Base::GetRegulator(this,regulator,parameters))!=NULL) {
    delete p_regulator;
    p_regulator=function;
  }
}

void XS_Base::Print()
{
  ATOOLS::msg.Out()<<m_name<<" {"<<m_colorscheme<<","<<m_helicityscheme
		   <<"} ("<<m_orderEW<<","<<m_orderQCD
		   <<")  ->  "<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb\n";
}

void XS_Base::UpdateIntegrator(Multi_Channel *&mc)
{
  if (m_nout<=2 || m_colorscheme!=cls::sample ||
      (!m_psmc && mc!=NULL)) return;
  int otfcs(0);
  Data_Reader read(" ",";","!","=");
  if (!read.ReadFromFile(otfcs,"CS_OTF")) otfcs=0;
  else msg_Info()<<METHOD<<"(): Set on the flight sampling "<<otfcs<<".\n";
  if (m_nout<7 && otfcs==0) 
    p_activepshandler->ColorIntegrator()->Initialize();
  p_activepshandler->ColorIntegrator()->SetOn(true);
  PreSample_Multi_Channel *psmc(dynamic_cast<PreSample_Multi_Channel*>(mc));
  Channel_Vector chs;
  if (psmc!=NULL) chs=psmc->ExtractChannels();
  Sample_Multi_Channel *smc
    (new Sample_Multi_Channel
     (m_nin,m_nout,p_flavours,p_activepshandler->ColorIntegrator()));
  smc->Initialize(chs);
  smc->Reset();
  if (mc!=NULL) {
    smc->SetValidN(mc->ValidN());
    smc->SetN(mc->N());
    mc->DropAllChannels();
    delete mc;
  }
  mc = smc;
  Reset();
}

void XS_Base::FillSIntegrator(Multi_Channel *&mc)
{
  if (m_nout>2 && m_colorscheme==cls::sample) {
    if (!m_psmc) {
      delete mc;
      mc=NULL;
      UpdateIntegrator(mc);
      return;
    }
    Color_Integrator *colint(p_activepshandler->ColorIntegrator());
    VHAAG *first(NULL);
    PreSample_Multi_Channel *psmc
      (new PreSample_Multi_Channel("Color Sample VHAAG",colint));
    for (size_t i(0);i<(m_nin+m_nout)/2;++i) {
      colint->GenerateType(i,true);
      Idx_Vector perm(colint->Orders().front());
      Multi_Channel *smc(new Multi_Channel(ToString(i),i));
      Idx_Vector rperm(perm.size());
      rperm.front()=perm.front();
      for (size_t j(1);j<rperm.size();++j) rperm[j]=perm[perm.size()-j];
      std::vector<size_t> hp(perm.size());
      for (size_t j(0);j<perm.size();++j) hp[j]=perm[j];
      VHAAG *sc(new VHAAG(m_nin,m_nout,hp,first));
      if (first==NULL) first=sc;
      smc->Add(sc);
      for (size_t j(0);j<perm.size();++j) hp[j]=rperm[j];
      sc = new VHAAG(m_nin,m_nout,hp,first);
      smc->Add(sc);
      psmc->AddMC(smc);
    }
    mc->DropAllChannels();
    delete mc;
    mc = psmc;
    colint->SetOn(false);
    return;
  }
  int otfcs(0);
  Data_Reader read(" ",";","!","=");
  if (!read.ReadFromFile(otfcs,"CS_OTF")) otfcs=0;
  else msg_Info()<<METHOD<<"(): Set on the flight sampling "<<otfcs<<".\n";
  if (m_colorscheme==cls::sample && otfcs==0) 
    p_activepshandler->ColorIntegrator()->Initialize();
  mc->DropAllChannels();
  if (p_isrhandler->KMROn()>0) {
    mc->Add(new PHASIC::T2Channel(m_nin,m_nout,p_flavours));
    mc->Add(new PHASIC::T3Channel(m_nin,m_nout,p_flavours));
  }
  else {
    mc->Add(new PHASIC::S1Channel(m_nin,m_nout,p_flavours));
    mc->Add(new PHASIC::T1Channel(m_nin,m_nout,p_flavours));
    mc->Add(new PHASIC::U1Channel(m_nin,m_nout,p_flavours));
  }
  for (size_t i=0;i<m_resonances.size();++i) {
    mc->Add(new PHASIC::S1Channel(m_nin,m_nout,p_flavours,m_resonances[i]));
  }
}      

void XS_Base::CreateISRChannels() 
{
}

void XS_Base::PrepareTerminate()  
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

