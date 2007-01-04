#include "Decay_Table.H"
#include "Random.H"
#include "Message.H"

using namespace ATOOLS;
using namespace std;

Decay_Channel::Decay_Channel(const Flavour & _flin) :
  m_metype(string("")), m_psfile(string("")), 
  m_width(0.), m_minmass(0.), m_flin(_flin)
 { }

Decay_Channel::Decay_Channel(const Decay_Channel & _dec) :
  m_processname(_dec.m_processname),
  m_metype(_dec.m_metype), m_psfile(_dec.m_psfile),
  m_width(_dec.m_width), m_minmass(_dec.m_minmass), 
  m_flin(_dec.m_flin), m_flouts(_dec.m_flouts) { }

void Decay_Channel::SetProcessName(const std::string _name) 
{ 
  if (_name!=string("")) {
    m_processname = _name;
    return;
  }
  m_processname = m_flin.Name()+string("->");
  for (FlSetConstIter fl=m_flouts.begin();fl!=m_flouts.end();++fl) 
    m_processname += fl->Name()+string("_");
  m_processname = m_processname.substr(0,m_processname.size()-1);
}

void Decay_Channel::Output() const
{
  msg.Out()<<m_flin<<" -> ";
  for (FlSetConstIter fl=m_flouts.begin();fl!=m_flouts.end();++fl) msg.Out()<<(*fl)<<" ";
  msg.Out()<<" : "<<m_width<<" GeV";
  if (m_metype!=string("")) msg.Out()<<", ME : "<<m_metype;
  if (m_psfile!=string("")) msg.Out()<<", PS : "<<m_psfile;
  msg.Out()<<"."<<endl;
}






Decay_Table::Decay_Table(const Flavour _flin) :
  m_overwrite(0), m_smearing(0), m_fixdecay(0), m_width(0.), m_flin(_flin)
{ }

void Decay_Table::AddDecayChannel(Decay_Channel * _dc)
{
  m_channels.push_back(_dc);
  m_width += _dc->Width();
}

void Decay_Table::SetSelectedChannel(const FlavourSet & _flouts,const bool bar)
{
  m_fixdecay         = true;
  Decay_Channel * dc = new Decay_Channel(m_flin);
  dc->SetWidth(0.);
  for (FlSetConstIter flit=_flouts.begin();flit!=_flouts.end();++flit) {
    if (bar) dc->AddDecayProduct((*flit).Bar());
        else dc->AddDecayProduct((*flit));
  }
  if (bar) {
    m_selectedchannelsbar.push_back(dc);
  }
  else {
    m_selectedchannels.push_back(dc);
  }
  m_channels.push_back(dc);
}

void Decay_Table::Select(const int flag) {
  if ((flag==1 && !m_selectedchannels.empty()) ||
      (flag==-1 && !m_selectedchannelsbar.empty())) {
    if (m_fixdecay) {
      if (flag==1) {
	p_selected  = (*m_seliter);
	m_seliter++;
	if (m_seliter==m_selectedchannels.end()) m_seliter = m_selectedchannels.begin();
      }
      if (flag==-1) {
	p_selected  = (*m_seliterbar);
	m_seliterbar++;
	if (m_seliterbar==m_selectedchannelsbar.end()) m_seliterbar = m_selectedchannelsbar.begin();
      }
      return;
    }
  }
  p_selected = NULL;
  if (m_channels.size()==1) {
    p_selected = m_channels[0];
    return;
  }
  double disc = m_width*ran.Get();
  for (size_t i=0;i<m_channels.size();++i) {
    disc -= m_channels[i]->Width();
    if (disc<0) {
      p_selected = m_channels[i];
      break;
    }
  }
}

void Decay_Table::Reset() {
  if (!m_selectedchannels.empty())    m_seliter    = m_selectedchannels.begin();
  if (!m_selectedchannelsbar.empty()) m_seliterbar = m_selectedchannelsbar.begin(); 
}

void Decay_Table::Output() {
  msg.Out()<<"Decay table for : "<<m_flin<<", total width is now "<<m_width<<" GeV,"<<endl
	   <<"   (instead of "<<m_flin.Width()<<" GeV), calculated by "<<m_generator<<endl
	   <<"----------------------------------------------------------------"<<endl;
  for (size_t i=0;i<m_channels.size();i++) m_channels[i]->Output();
  if (m_overwrite) msg.Out()<<" Value of Particle.dat has been overwritten by "
			    <<m_generator<<"."<<endl;
  msg.Out()<<"----------------------------------------------------------------"<<endl;
}

double Decay_Table::Width(const int i)
{
  if (i<0 || i>= (int)m_channels.size()) {
    msg.Error()<<"Error in Decay_Table::Width("<<i<<")."<<endl
	       <<"   Out of bounds : 0 ... "<<m_channels.size()<<" ."<<endl
	       <<"   Return 0."<<endl;
    return 0.;
  }
  return m_channels[i]->Width();
}

Decay_Channel * Decay_Table::GetDecayChannel(const int i)
{
  if (i<0 || i>= (int)m_channels.size()) {
    msg.Error()<<"Error in Decay_Table::Width("<<i<<")."<<endl
	       <<"   Out of bounds : 0 ... "<<m_channels.size()<<" ."<<endl
	       <<"   Return NULL"<<endl;
    return NULL;
  }
  return m_channels[i];
}

Decay_Channel * Decay_Table::GetOneDecayChannel()     
{ 
  if (p_selected) return p_selected; 
  return NULL;
}


double Decay_Table::Width(const FlavourSet) { return 0.; }
Decay_Channel * Decay_Table::GetDecayChannel(const FlavourSet) { return NULL; }
