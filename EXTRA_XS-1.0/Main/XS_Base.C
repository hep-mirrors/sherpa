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
#include "XS_Model_Handler.H"
#include "Data_Reader.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;

bool XS_Base::s_sortflavours(true);

XS_Base::XS_Base():
  p_colours(NULL), 
  m_channels(false), m_ownmodel(false),
  p_model(NULL)
{
  m_name="Empty XS";
  m_scale[PHASIC::stp::sfs]=m_scale[PHASIC::stp::sis]=
    std::numeric_limits<double>::max();
  m_scale[stp::ren]=m_scale[stp::fac]=sqr(rpa.gen.Ecms());
} 

XS_Base::XS_Base(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours,
		 const PHASIC::scl::scheme scalescheme,const int kfactorscheme,
		 BEAM::Beam_Spectra_Handler *const beamhandler,
		 PDF::ISR_Handler *const isrhandler,
		 ATOOLS::Selector_Data *const selectordata,XS_Model_Base *const model):
  Integrable_Base(nin,nout,scalescheme,kfactorscheme,
		  beamhandler,isrhandler,selectordata),
  p_colours(NULL), 
  m_channels(false), m_ownmodel(false),
  p_model(model)
{
  Init(flavours);
  Initialize(scalescheme,kfactorscheme,beamhandler,isrhandler,selectordata);
  m_scale[stp::ren]=m_scale[stp::fac]=sqr(rpa.gen.Ecms());
}

XS_Base::XS_Base(const size_t nin,const size_t nout,
		 const ATOOLS::Flavour *flavours,XS_Model_Base *const model):
  Integrable_Base(nin,nout),
  p_colours(NULL), 
  m_channels(false), m_ownmodel(false),
  p_model(model)
{
  Init(flavours);
  p_selector = new ATOOLS::No_Selector();
  m_scale[PHASIC::stp::sfs]=m_scale[PHASIC::stp::sis]=
    std::numeric_limits<double>::max();
  m_scale[stp::ren]=m_scale[stp::fac]=sqr(rpa.gen.Ecms());
}

XS_Base::~XS_Base() 
{
  if (p_colours!=NULL) { 
    for (size_t i=0;i<m_nin+m_nout;++i) delete [] p_colours[i];
    delete [] p_colours;
  }
  if (m_ownmodel && p_model!=NULL) delete p_model;
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
    if (s_sortflavours) {
      SortFlavours(p_flavours,m_nin);
      SortFlavours(&p_flavours[m_nin],m_nout);
    }
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

void XS_Base::InitializeModel(MODEL::Model_Base *const model,
			      const std::string &file)
{
  if (m_ownmodel && p_model!=NULL) delete p_model;
  Data_Read read(file);
  std::string mt(read.GetValue<std::string>("SIGNAL_MODEL","SM"));
  p_model = XS_Model_Handler::GetModel(mt);
  if (p_model==NULL) THROW(fatal_error,"Interaction model not implemented");
  p_model->Initialize(model,file);
  m_ownmodel=true;
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
    msg_Error()<<"XS_Base::ResetSelector("<<selectordata<<"): "
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
  msg_Out()<<m_name<<" {"<<m_colorscheme<<","<<m_helicityscheme
		   <<"} ("<<m_orderEW<<","<<m_orderQCD
		   <<")  ->  "<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb\n";
}

void XS_Base::PrepareTerminate()  
{
  if (rpa.gen.BatchMode()) return;
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

void XS_Base::SetCoreMaxJetNumber(const int &n)
{
  m_coremaxjetnumber=n;
}

void XS_Base::CreateFSRChannels() 
{
  p_pshandler->FSRIntegrator()->DropAllChannels();
  p_pshandler->FSRIntegrator()->
    Add(new PHASIC::S1Channel(m_nin,m_nout,p_flavours));
  p_pshandler->FSRIntegrator()->
    Add(new PHASIC::T1Channel(m_nin,m_nout,p_flavours));
  p_pshandler->FSRIntegrator()->
    Add(new PHASIC::U1Channel(m_nin,m_nout,p_flavours));
  for (size_t i=0;i<m_resonances.size();++i) {
    p_pshandler->FSRIntegrator()->
      Add(new PHASIC::S1Channel(m_nin,m_nout,p_flavours,m_resonances[i]));
  }
}      

void XS_Base::CreateISRChannels() 
{
}

void XS_Base::SetScales(const double &scale)
{
  m_scale[stp::ren]=m_scale[stp::fac]=scale;
}
