#include "Hadron_Decay_Channel.H"
#include "HD_ME_Selector.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

Hadron_Decay_Channel::Hadron_Decay_Channel(Decay_Channel * _dc) :
  p_dc(_dc), p_me(NULL), m_metype(p_dc->METype()),
  Integrable_Base(1,_dc->NumberOfDecayProducts())
{
  m_resultpath      = string("./");
  m_resultfile      = m_histofile = string("");
  m_name            = p_dc->ProcessName();
  p_flavours        = new Flavour[m_nin+m_nout];
  p_flavours[0]     = p_dc->GetDecaying();
  for (int i=0;i<p_dc->NumberOfDecayProducts();i++) {
    p_flavours[i+1] = p_dc->GetDecayProduct(i);
  }
  HD_ME_Selector mesel;
  p_me = mesel.GetME(m_nin,m_nout,p_flavours,m_metype);
}

Hadron_Decay_Channel::~Hadron_Decay_Channel()
{
  if (p_dc) { delete p_dc; p_dc=NULL; }
}


void Hadron_Decay_Channel::InitialisePhaseSpace(vector<string> & PStype) 
{
  bool mustinit;
  p_ps = new HD_PS_Base(this,PStype,mustinit);
  if (mustinit) p_ps->Initialise();
}


double Hadron_Decay_Channel::Differential()
{
  p_momenta[0] = Vec4D(p_flavours[0].Mass(),0.,0.,0.);
  p_ps->GeneratePoint(p_momenta);
  p_ps->GenerateWeight(p_momenta);
  double value = (*p_me)(p_momenta), weight = p_ps->Weight();
  return value*weight;
}

