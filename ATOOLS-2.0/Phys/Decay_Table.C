#include "Decay_Table.H"
#include "Random.H"
#include "Message.H"

using namespace ATOOLS;
using namespace std;

Decay_Channel::Decay_Channel(const Flavour & _flin) :
  m_flin(_flin), m_width(0.), m_minmass(0.), m_processname(string("")) { }

Decay_Channel::Decay_Channel(const Decay_Channel & _dec) :
  m_flin(_dec.m_flin), m_flouts(_dec.m_flouts), 
  m_width(_dec.m_width), m_minmass(_dec.m_minmass), 
  m_processname(_dec.m_processname) { }

void Decay_Channel::Output() const
{
  msg.Out()<<m_flin<<" -> ";
  for (FlSetIter fl=m_flouts.begin();fl!=m_flouts.end();++fl) msg.Out()<<(*fl)<<" ";
  msg.Out()<<" : "<<m_width<<" GeV."<<endl;
}






Decay_Table::Decay_Table(const Flavour _flin) :
  m_flin(_flin), m_width(0.), m_overwrite(0), m_smearing(0), m_fixdecay(0),
  m_generator(std::string("")) { }

void Decay_Table::AddDecayChannel(Decay_Channel * _dc)
{
  m_channels.push_back(_dc);
  m_width += _dc->Width();
}

void Decay_Table::SetSelectedChannel(const FlavourSet & _flouts)
{
  m_flouts           = _flouts;
  m_fixdecay         = true;
  Decay_Channel * dc = new Decay_Channel(m_flin);
  for (FlSetIter flit=m_flouts.begin();flit!=m_flouts.end();++flit) dc->AddDecayProduct((*flit));
  dc->SetWidth(0.);
  m_channels.push_back(dc);
  p_selected         = dc;
}

void Decay_Table::Select() {
  if (m_fixdecay) {
    if (p_selected==NULL) {
      for (int i=0;i<m_channels.size();i++) {
	if (m_channels[i]->GetDecayProducts()==m_flouts) {
	  p_selected = m_channels[i];
	  break;
	}
      }
    }
  }
  else {
    p_selected = NULL;
    if (m_channels.size()==1) {
      p_selected = m_channels[0];
      return;
    }
    double disc = m_width*ran.Get();
    for (int i=0;i<m_channels.size();i++) {
      disc -= m_channels[i]->Width();
      if (disc<0) {
	p_selected = m_channels[i];
	break;
      }
    }
  }
}

void Decay_Table::Output() {
  msg.Out()<<"Decay table for : "<<m_flin<<", total width is now "<<m_width<<" GeV,"<<endl
	   <<"   (instead of "<<m_flin.Width()<<" GeV), calculated by "<<m_generator<<endl
	   <<"----------------------------------------------------------------"<<endl;
  for (int i=0;i<m_channels.size();i++) m_channels[i]->Output();
  if (m_overwrite) msg.Out()<<" Value of Particle.dat has been overwritten by "
			    <<m_generator<<"."<<endl;
  msg.Out()<<"----------------------------------------------------------------"<<endl;
}

double Decay_Table::Width(const int i)
{
  if (i<0 || i>= m_channels.size()) {
    msg.Error()<<"Error in Decay_Table::Width("<<i<<")."<<endl
	       <<"   Out of bounds : 0 ... "<<m_channels.size()<<" ."<<endl
	       <<"   Return 0."<<endl;
    return 0.;
  }
  return m_channels[i]->Width();
}

Decay_Channel * Decay_Table::GetDecayChannel(const int i)
{
  if (i<0 || i>= m_channels.size()) {
    msg.Error()<<"Error in Decay_Table::Width("<<i<<")."<<endl
	       <<"   Out of bounds : 0 ... "<<m_channels.size()<<" ."<<endl
	       <<"   Return NULL"<<endl;
    return NULL;
  }
  return m_channels[i];
}

double Decay_Table::Width(const FlavourSet) { return 0.; }
Decay_Channel * Decay_Table::GetDecayChannel(const FlavourSet) { return NULL; }
