#include "Decay_Table.H"
#include "Random.H"
#include "Message.H"

using namespace ATOOLS;
using namespace std;

Decay_Channel::Decay_Channel(const Flavour & _flin) :
  m_metype(string("")), m_psfile(string("")), 
  m_width(0.), m_deltawidth(-1.), m_minmass(0.), m_origin(""), m_flin(_flin)
 { }

Decay_Channel::Decay_Channel(const Decay_Channel & _dec) :
  m_processname(_dec.m_processname),
  m_metype(_dec.m_metype), m_psfile(_dec.m_psfile),
  m_width(_dec.m_width), m_deltawidth(_dec.m_deltawidth), 
  m_minmass(_dec.m_minmass), m_origin(_dec.m_origin),
  m_flin(_dec.m_flin), m_flouts(_dec.m_flouts) { }

void Decay_Channel::SetProcessName(const std::string _name) 
{ 
  if (_name!=string("")) {
    m_processname = _name;
    return;
  }
  m_processname = m_flin.IDName()+string("->");
  for (FlSetConstIter fl=m_flouts.begin();fl!=m_flouts.end();++fl) 
    m_processname += fl->IDName()+string("_");
  m_processname = m_processname.substr(0,m_processname.size()-1);
}

void Decay_Channel::Output() const
{
  msg_Out()<<m_flin<<" -> ";
  for (FlSetConstIter fl=m_flouts.begin();fl!=m_flouts.end();++fl) msg_Out()<<(*fl)<<" ";
  msg_Out()<<" : "<<m_width;
  if (m_deltawidth>=0.) msg_Out()<<" (+/-"<<m_deltawidth<<")";
  msg_Out()<<" GeV";
  if (m_metype!=string("")) msg_Out()<<", ME : "<<m_metype;
  if (m_psfile!=string("")) msg_Out()<<", PS : "<<m_psfile;
  msg_Out()<<"."<<endl;
}






Decay_Table::Decay_Table(const Flavour _flin) :
  m_overwrite(0), m_smearing(0), m_fixdecay(0), m_width(0.), m_flin(_flin)
{ }

void Decay_Table::AddDecayChannel(Decay_Channel * _dc)
{
  for(size_t i=0;i<m_channels.size();i++) {
    if(m_channels[i]->GetDecayProducts()==_dc->GetDecayProducts() &&
       _dc->Width()!=0.0) {
      msg_Error()<<METHOD<<" Warning: Duplicate decaychannel: ";
      _dc->Output();
    }
  }
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
	if (m_seliter==m_selectedchannels.end()) 
	  m_seliter = m_selectedchannels.begin();
      }
      if (flag==-1) {
	p_selected  = (*m_seliterbar);
	m_seliterbar++;
	if (m_seliterbar==m_selectedchannelsbar.end()) 
	  m_seliterbar = m_selectedchannelsbar.begin();
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
  msg_Out()<<"Decay table for : "<<m_flin<<", total width is now "<<m_width<<" GeV,"<<endl
	   <<"   (instead of "<<m_flin.Width()<<" GeV), calculated by "<<m_generator<<endl
	   <<"----------------------------------------------------------------"<<endl;
  for (size_t i=0;i<m_channels.size();i++) {
    m_channels[i]->Output();
    if( m_channels[i]->DeltaWidth() >= 0. ) {
      double wanted_br = m_channels[i]->Width()/m_flin.Width()*100.;
      double upper_br  = (m_channels[i]->Width()+m_channels[i]->DeltaWidth())/m_flin.Width()*100.;
      double lower_br  = (m_channels[i]->Width()-m_channels[i]->DeltaWidth())/m_flin.Width()*100.;
      double exp_br    = m_channels[i]->Width()/m_width*100.;
      if( exp_br > upper_br+Accu() || exp_br < lower_br-Accu() ) {  
        msg_Out()<<om::red<<"     WARNING: branching ratio "
          <<exp_br<<"% is out of bounds ("<<wanted_br <<" +/- "<<(upper_br-lower_br)/2.<<" %)."<<om::reset<<endl;
      }
    }
  }
  if (m_overwrite) msg_Out()<<" Value of Particle.dat has been overwritten by "
			    <<m_generator<<"."<<endl;
  msg_Out()<<"----------------------------------------------------------------"<<endl;
}

void Decay_Table::UpdateWidth() {
  m_width = 0.;
  for (size_t i=0;i<m_channels.size();i++) m_width+= m_channels[i]->Width();
}

void Decay_Table::ScaleToWidth() {
  if(m_flin.Width()/m_width!=1.0) {
    double delta_tot(0.0);
    for (size_t i=0;i<m_channels.size();i++)
      delta_tot+=m_channels[i]->DeltaWidth();
    if (delta_tot>0.0) {
      for (size_t i=0;i<m_channels.size();i++) {
        double scale_fac=m_channels[i]->DeltaWidth()/delta_tot;
        m_channels[i]->SetWidth(m_channels[i]->Width()+
                                scale_fac*(m_flin.Width()-m_width));
      }
      UpdateWidth();
    }
  }
}

double Decay_Table::Width(const int i)
{
  if (i<0 || i>= (int)m_channels.size()) {
    msg_Error()<<"Error in Decay_Table::Width("<<i<<")."<<endl
	       <<"   Out of bounds : 0 ... "<<m_channels.size()<<" ."<<endl
	       <<"   Return 0."<<endl;
    return 0.;
  }
  return m_channels[i]->Width();
}

Decay_Channel * Decay_Table::GetDecayChannel(const int i)
{
  if (i<0 || i>= (int)m_channels.size()) {
    msg_Error()<<"Error in Decay_Table::Width("<<i<<")."<<endl
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

Decay_Channel * Decay_Table::GetDecayChannel(const FlavourSet decayproducts)
{
  for(size_t i=0;i<m_channels.size();i++) {
    if(m_channels[i]->GetDecayProducts() == decayproducts) return m_channels[i];
  }
  return NULL;
}

void Decay_Table::EraseDecayChannel(const int i) {
  delete m_channels[i];
  for (size_t j=i;j<m_channels.size()-1;j++) m_channels[j] = m_channels[j+1];
  m_channels.pop_back();
}
