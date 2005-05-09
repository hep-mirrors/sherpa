#include "Hadron_Decay_Channel.H"
#include "HD_ME_Selector.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

Hadron_Decay_Channel::Hadron_Decay_Channel(Decay_Channel * _dc) :
  p_dc(_dc), p_me(NULL), m_metype(p_dc->METype()),
  Integrable_Base(1,_dc->NumberOfDecayProducts())
{
  m_resultpath      = string("./");								// where to write results
  m_resultfile      = m_histofile = string("");					// filename
  m_name            = p_dc->ProcessName();						// name of proces
  p_flavours        = new Flavour[m_nin+m_nout];					
  p_flavours[0]     = p_dc->GetDecaying();						// decaying particle
  for (int i=0;i<p_dc->NumberOfDecayProducts();i++) {			// decay products
    p_flavours[i+1] = p_dc->GetDecayProduct(i);			
  }
  HD_ME_Selector mesel;											// ME selector
  p_me = mesel.GetME(m_nin,m_nout,p_flavours,m_metype);			// get the appropr. ME
}

Hadron_Decay_Channel::~Hadron_Decay_Channel()
{
  if (p_dc) { delete p_dc; p_dc=NULL; }
  if (p_ps) { delete p_ps; p_ps=NULL; }
  if (p_me) { delete p_me; p_me=NULL; }
}


void Hadron_Decay_Channel::InitialisePhaseSpace(vector<string> & PStype) 
{
  bool mustinit;
  struct Model locmd;
  p_ps = new HD_PS_Base(this,PStype,mustinit,locmd);					// new PS
  p_me->SetModelParameters( locmd );
  if (mustinit) p_ps->Initialise();								// if in need => ini
}


double Hadron_Decay_Channel::Differential()
{
  p_momenta[0] = Vec4D(p_flavours[0].Mass(),0.,0.,0.);			// decay from rest
  p_ps->GeneratePoint(p_momenta);								// generate a PS point
  p_ps->GenerateWeight(p_momenta);								// calculate its weight factor
  double value = (*p_me)(p_momenta), weight = p_ps->Weight();	// get ME value & weight f
  return value*weight;
}

