#include "Spectrum_Generator_Base.H"
#include "Message.H"

using namespace MODEL;
using namespace ATOOLS;


DecayChannel::DecayChannel(const Flavour _flin) : 
  m_flin(_flin), m_width(0.) 
{ }

void DecayChannel::AddDecayProduct(const Flavour _flout) {
  m_flouts.push_back(_flout);
}

void DecayChannel::SetDecayWidth(const double _width)
{
  m_width = _width;
}

void DecayChannel::Output()
{
  msg.Out()<<" "<<m_flin<<" -> ";
  for (int i=0;i<m_flouts.size();++i) msg.Out()<<m_flouts[i]<<" ";
  msg.Out()<<" : "<<m_width<<" GeV."<<std::endl;
}

int DecayChannel::Nout() {
  return m_flouts.size();
}

Flavour DecayChannel::Flin() {
  return m_flin;
}

Flavour DecayChannel::Flout(const int i) {
  if (i>-1 && i<m_flouts.size()) return m_flouts[i];
  msg.Error()<<"Error in DecayChannel::Flout("<<i<<")"<<std::endl
	     <<"   Try to access out flavour for "<<m_flin<<". Return none. "<<std::endl;
  return Flavour(kf::none);
}

double DecayChannel::DecayWidth() {
  return m_width;
}




Spectrum_Generator_Base::Spectrum_Generator_Base(ATOOLS::Data_Read * _dataread,Model_Base * _model) :
  p_dataread(_dataread), p_model(_model)  {};

Spectrum_Generator_Base::~Spectrum_Generator_Base() 
{ 
  if (!decayproducts.empty()) 
    decayproducts.erase(decayproducts.begin(),decayproducts.end());
}




