#include "Decay_Table.H"
#include "Message.H"

using namespace ATOOLS;

Decay_Channel::Decay_Channel(const Flavour _flin) :
  m_flin(_flin), m_width(0.) { }

void Decay_Channel::Output() 
{
  msg.Out()<<"  "<<m_flin<<" -> ";
  for (FlSetIter fl=m_flouts.begin();fl!=m_flouts.end();++fl) msg.Out()<<(*fl)<<" ";
  msg.Out()<<" : "<<m_width<<" GeV."<<endl;
}

Decay_Table::Decay_Table(const Flavour _flin) :
  m_flin(_flin), m_width(0.) { }

void Decay_Table::AddDecayChannel(Decay_Channel * _dc)
{
  m_channels.push_back(_dc);
  m_width += _dc->Width();
}

void Decay_Table::Output() {
  msg.Out()<<"Decay table for : "<<m_flin<<", total width : "<<m_width<<" GeV."<<endl
	   <<"----------------------------------------------------------------"<<endl;
  for (int i=0;i<m_channels.size();i++) m_channels[i]->Output();
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
