#include "Integrable_Base.H"

#include "Message.H"

using namespace PHASIC;

Integrable_Base::Integrable_Base(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours,
				 const int scalescheme,const int kfactorscheme,const double scalefactor,
				 BEAM::Beam_Spectra_Handler *const beamhandler,
				 PDF::ISR_Handler *const isrhandler,
				 ATOOLS::Selector_Data *const selectordata):
  m_name(""), m_nin(nin), m_nout(nout), m_nvector(nin+nout), p_flavours(NULL), 
  p_momenta(new ATOOLS::Vec4D[nin+nout]), m_scalescheme(scalescheme), m_kfactorscheme(kfactorscheme), 
  m_scalefactor(scalefactor), m_threshold(0.), m_overflow(0.), 
  m_xinfo(std::vector<double>(4)),
  m_n(0), m_last(0.), m_lastlumi(0.), m_lastdxs(0.), m_max(0.),
  m_totalxs(0.),m_totalsum (0.), m_totalsumsqr(0.), m_totalerr(0.), 
  m_swaped(false), p_selected(this), p_beamhandler(beamhandler), p_isrhandler(isrhandler), 
  p_pshandler(NULL), p_activepshandler(NULL), p_selector(NULL), p_cuts(NULL) {}

Integrable_Base::~Integrable_Base()
{
  if (p_selector!=NULL) delete p_selector;
  if (p_flavours!=NULL) delete [] p_flavours;
  if (p_momenta!=NULL) delete [] p_momenta;
}

Integrable_Base *const Integrable_Base::Selected()
{ 
  if (p_selected!=this) return p_selected->Selected(); 
  return this; 
}

void Integrable_Base::SetMomenta(const ATOOLS::Vec4D *momenta) 
{ 
  for (size_t i=0;i<NVector();++i) p_momenta[i]=momenta[i];
}

void Integrable_Base::SetMax(const double max) 
{
  m_max=max;
} 

void Integrable_Base::SetMax() 
{
  ATOOLS::msg.Error()<<"Integrable_Base::SetMax(): Virtual function called !"<<std::endl;
} 

bool Integrable_Base::OneEvent() 
{
  ATOOLS::msg.Error()<<"Integrable_Base::OneEvent(): Virtual function called !"<<std::endl;
  return false;
} 

bool Integrable_Base::SameEvent() 
{
  ATOOLS::msg.Error()<<"Integrable_Base::SameEvent(): Virtual function called !"<<std::endl;
  return false;
} 

ATOOLS::Blob_Data_Base *Integrable_Base::WeightedEvent() 
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



