#include "PHASIC++/Decays/Decay_Table.H"
#include "PHASIC++/Decays/Decay_Channel.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Blob.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Decay_Table::Decay_Table(const Flavour _flin, const ATOOLS::Mass_Selector* ms) :
  vector<Decay_Channel*>(0),
  m_activewidth(0.), m_totalwidth(0.), m_flin(_flin), p_ms(ms)
{ }

Decay_Table::~Decay_Table()
{
  for (size_t i=0; i<size(); i++) {
    delete at(i); at(i)=NULL;
  }
}

void Decay_Table::AddDecayChannel(Decay_Channel * _dc)
{
  for(size_t i=0;i<size();i++) {
    if(at(i)->Flavs()==_dc->Flavs() &&
       _dc->Width()!=0.0 && at(i)->Width()!=0.0 &&
       _dc->Active()>-1 && at(i)->Active()>-1) {
      msg_Error()<<METHOD<<" Warning: Duplicate decaychannel: ";
      _dc->Output();
      msg_Error()<<endl;
    }
  }
  push_back(_dc);
  if (_dc->Active()>=0) m_totalwidth += _dc->Width();
  if (_dc->Active()>0) m_activewidth += _dc->Width();
}

void Decay_Table::RemoveDecayChannel(size_t i)
{
  if (at(i)->Active()>=0) m_totalwidth -= at(i)->Width();
  if (at(i)->Active()>0) m_activewidth -= at(i)->Width();
  erase(begin()+i);
}

void Decay_Table::Output() {
  msg_Out()<<(*this);
}

namespace PHASIC {
  std::ostream &operator<<(std::ostream &os,const Decay_Table &dt)
  {
    os<<"Decay table for : "<<dt.m_flin<<"."<<endl;
    os<<setw(30)<<"Total width: "<<dt.TotalWidth()<<" GeV"<<endl;
    if (dt.ActiveWidth()!=dt.TotalWidth())
      os<<setw(30)<<"Active width: "<<dt.ActiveWidth()<<" GeV"<<endl;
    if (dt.Flav().Width()!=dt.TotalWidth())
      os<<setw(30)<<"Flavour width: "<<dt.Flav().Width()<<" GeV"<<endl;
    os<<"----------------------------------------"<<endl;
    for (size_t i=0;i<dt.size();i++) {
      if (dt.at(i)->Active()!=-1) {
	os<<*dt.at(i);
	if (dt.TotalWidth()>0. && dt.at(i)->Width()>0.) 
	  os<<", BR= "<<setw(5)<<(dt.at(i)->Width()/dt.TotalWidth()*100.)
	    <<" %";
	os<<endl;
      }
    }
    os<<"----------------------------------------"<<endl;
    return os;
  }
}

void Decay_Table::UpdateWidth() {
  m_activewidth = m_totalwidth = 0.;
  for (size_t i=0;i<size();i++) {
    if (at(i)->Active()>=0) m_totalwidth += at(i)->Width();
    if (at(i)->Active()>0) m_activewidth += at(i)->Width();
  }
}

void Decay_Table::UpdateWidth(Decay_Channel * hdc,const double &width)
{
  if (hdc->Active()>=0) m_totalwidth -= hdc->Width();
  if (hdc->Active()>0) m_activewidth -= hdc->Width();
  hdc->SetWidth(width);
  if (hdc->Active()>=0) m_totalwidth += hdc->Width();
  if (hdc->Active()>0) m_activewidth += hdc->Width();
}

Decay_Channel * Decay_Table::GetDecayChannel
    (const Flavour_Vector& flavs) const
{
  for(size_t i=0;i<size();i++) {
    if(at(i)->Flavs() == flavs && at(i)->Active()>0) {
      return at(i);
    }
  }
  return NULL;
}

void Decay_Table::EraseDecayChannel(const int i) {
  delete at(i);
  erase(begin()+i);
}

Decay_Channel* Decay_Table::Select() const
{
  DEBUG_FUNC(m_flin);
  Decay_Channel* selected(NULL);
  if (size()==1) return at(0);
  double disc = m_activewidth*ran->Get();
  for (size_t i=0;i<size();++i) {
    if (at(i)->Active()<1) continue;
    disc -= at(i)->Width();
    if (disc<0) {
      selected=at(i);
      break;
    }
  }
  DEBUG_VAR(*selected);
  return selected;
}
