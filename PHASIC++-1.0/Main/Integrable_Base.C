#include "Integrable_Base.H"
#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"

#include "Phase_Space_Handler.H"
#include "Message.H"

using namespace PHASIC;

Integrable_Base::Integrable_Base(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours,
				 const int scalescheme,const int kfactorscheme,const double scalefactor,
				 BEAM::Beam_Spectra_Handler *const beamhandler,
				 PDF::ISR_Handler *const isrhandler,
				 ATOOLS::Selector_Data *const selectordata):
  m_name(""), m_nin(nin), m_nout(nout), m_nvector(nin+nout), p_flavours(NULL), 
  p_momenta(new ATOOLS::Vec4D[nin+nout]), m_scalescheme(scalescheme), 
  m_kfactorscheme(kfactorscheme), m_scalefactor(scalefactor), m_threshold(0.), 
  m_overflow(0.), m_xinfo(std::vector<double>(4)),
  m_n(0), m_last(0.), m_lastlumi(0.), m_lastdxs(0.), m_max(0.),
  m_totalxs(0.),m_totalsum (0.), m_totalsumsqr(0.), m_totalerr(0.), 
  m_ssum(0.), m_ssumsqr(0.), m_smax(0.), m_ssigma2(0.), m_wmin(0.), m_sn(0), m_son(1), 
  m_swaped(false), p_selected(this), p_parent(this), 
  p_beamhandler(beamhandler), p_isrhandler(isrhandler), 
  p_pshandler(NULL), p_activepshandler(NULL), p_selector(NULL), p_cuts(NULL),
  m_ownselector(true) {}

Integrable_Base::~Integrable_Base()
{
  if (p_selector!=NULL && m_ownselector) delete p_selector;
  if (p_flavours!=NULL) delete [] p_flavours;
  if (p_momenta!=NULL) delete [] p_momenta;
}

Integrable_Base *const Integrable_Base::Selected()
{ 
  if (p_selected!=this && p_selected!=NULL) return p_selected->Selected();
  return this; 
}

Integrable_Base *const Integrable_Base::Parent()
{ 
  if (p_parent!=this && p_parent!=NULL) return p_parent->Parent();
  return this; 
}

double Integrable_Base::TotalResult()
{ 
//   if (m_ssigma2>0. && m_sn<1000) return m_totalsum/m_ssigma2; 
//   if (m_sn<1000) return m_ssum/m_sn;
//   double ssigma2 = (m_ssumsqr/m_sn - ATOOLS::sqr(m_ssum/m_sn))/(m_sn-1);
//   return (m_totalsum+m_ssum/ssigma2/m_sn)/(m_ssigma2+1./ssigma2);
  if (m_ssigma2>0. && m_sn<1000) return m_totalsum/m_ssigma2; 
  if (m_sn<1000) return m_ssum/m_sn;
  double ssigma2 = ATOOLS::sqr(m_ssum/m_sn)/((m_ssumsqr/m_sn - ATOOLS::sqr(m_ssum/m_sn))/(m_sn-1));
  return (m_totalsum+m_ssum*ssigma2/m_sn)/(m_ssigma2+ssigma2);
}

double Integrable_Base::TotalVar() 
{
  if (m_nin==1 && m_nout==2) return 0.;
  if (m_sn<1000) {
    if (m_ssigma2>0.) return m_totalsum/m_ssigma2/sqrt(m_ssigma2); 
    else return TotalResult(); 
  }

  double disc = ATOOLS::sqr(m_ssum/m_sn)/((m_ssumsqr/m_sn - ATOOLS::sqr(m_ssum/m_sn))/(m_sn-1));
  if (disc>0.) return TotalResult()/sqrt(m_ssigma2+disc);
  
  return m_totalsum/m_ssigma2/sqrt(m_ssigma2);
}

void Integrable_Base::SetMomenta(const ATOOLS::Vec4D *momenta) 
{ 
  if (!p_momenta) {
    ATOOLS::msg.Error()<<"Integrable_Base::SetMomenta("<<momenta<<"): "
		       <<"p_momenta = NULL. Abort."<<std::endl;
    abort();
  }
  for (size_t i=0;i<NVector();++i) p_momenta[i]=momenta[i];
}

void Integrable_Base::SetScale(const double scale) 
{ 
  std::cout<<"Integrable_Base::SetScale("<<scale<<"): Virtual function called."<<std::endl;
}

void Integrable_Base::SetMax(const double max, int depth) 
{
  if (max!=0.) m_max=max;
} 

void Integrable_Base::SetMax() 
{
  ATOOLS::msg.Error()<<"Integrable_Base::SetMax(): Virtual function called !"<<std::endl;
} 

void Integrable_Base::ResetMax(int) 
{
  ATOOLS::msg.Error()<<"Integrable_Base::ResetMax(): Virtual function called !"<<std::endl;
} 

bool Integrable_Base::OneEvent() 
{
  ATOOLS::msg.Error()<<"Integrable_Base::OneEvent(): Virtual function called !"<<std::endl;
  return false;
} 

bool Integrable_Base::Trigger(const ATOOLS::Vec4D *const momenta)
{
  return p_selector->Trigger(momenta);
} 

bool Integrable_Base::OneEvent(const double mass,const int mode) 
{
  return p_activepshandler->OneEvent(mass,mode);
} 

bool Integrable_Base::SameEvent() 
{
  return p_activepshandler->SameEvent();
  ATOOLS::msg.Error()<<"Integrable_Base::SameEvent(): Virtual function called !"<<std::endl;
  return false;
} 

ATOOLS::Blob_Data_Base *Integrable_Base::WeightedEvent(const int mode) 
{
  ATOOLS::msg.Error()<<"Integrable_Base::WeightedEvent(): Virtual function called !"<<std::endl;
  return NULL;
} 

ATOOLS::Blob_Data_Base *Integrable_Base::SameWeightedEvent() 
{
  ATOOLS::msg.Error()<<"Integrable_Base::SameWeightedEvent(): Virtual function called !"<<std::endl;
  return NULL;
} 

void Integrable_Base::SetPSHandler(Phase_Space_Handler *const pshandler) 
{
  p_activepshandler=pshandler;
} 

void Integrable_Base::OptimizeResult()
{
  ATOOLS::msg.Error()<<"Integrable_Base::OptimizeResult(): Virtual function called !"<<std::endl;
} 
