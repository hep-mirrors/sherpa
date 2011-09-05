#include "ATOOLS/Phys/Decay_Table.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;
using namespace std;

Decay_Table::Decay_Table(const Flavour _flin) :
  vector<Decay_Channel*>(0),
  m_width(0.), m_flin(_flin)
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
    if(at(i)->GetDecayProducts()==_dc->GetDecayProducts() &&
       _dc->Width()!=0.0) {
      msg_Error()<<METHOD<<" Warning: Duplicate decaychannel: ";
      _dc->Output();
    }
  }
  push_back(_dc);
  m_width += _dc->Width();
}

void Decay_Table::Output() {
  msg_Out()<<(*this);
}

std::ostream &ATOOLS::operator<<(std::ostream &os,const Decay_Table &dt)
{
  os<<"Decay table for : "<<dt.m_flin<<", total width is now "
    <<dt.m_width<<" GeV,"<<endl
    <<"   (instead of "<<dt.m_flin.Width()<<" GeV)"<<endl
    <<"----------------------------------------"<<endl;
  for (size_t i=0;i<dt.size();i++) {
    os<<*(dt.at(i));
    if( dt.at(i)->DeltaWidth() >= 0. ) {
      double wanted_br=dt.at(i)->Width()/dt.m_flin.Width()*100.;
      double upper_br=(dt.at(i)->Width()+dt.at(i)->DeltaWidth())/
        dt.m_flin.Width()*100.;
      double lower_br=(dt.at(i)->Width()-dt.at(i)->DeltaWidth())/
        dt.m_flin.Width()*100.;
      double exp_br=dt.at(i)->Width()/dt.m_width*100.;
      if( exp_br > upper_br+Accu() || exp_br < lower_br-Accu() ) {  
        msg_Out()<<om::red<<"     WARNING: branching ratio "
                 <<exp_br<<"% is out of bounds ("<<wanted_br <<" +/- "
                 <<(upper_br-lower_br)/2.<<" %)."<<om::reset<<endl;
      }
    }
  }
  os<<"----------------------------------------"<<endl;
  return os;
}

void Decay_Table::UpdateWidth() {
  m_width = 0.;
  for (size_t i=0;i<size();i++) m_width+= at(i)->Width();
}

void Decay_Table::ScaleToWidth() {
  if(m_flin.Width()/m_width!=1.0) {
    double delta_tot(0.0);
    for (size_t i=0;i<size();i++)
      delta_tot+=at(i)->DeltaWidth();
    if (delta_tot>0.0) {
      for (size_t i=0;i<size();i++) {
        double scale_fac=at(i)->DeltaWidth()/delta_tot;
        at(i)->SetWidth(at(i)->Width()+
                                scale_fac*(m_flin.Width()-m_width));
      }
      UpdateWidth();
    }
  }
}

Decay_Channel * Decay_Table::GetDecayChannel(const FlavourSet decayproducts)
{
  for(size_t i=0;i<size();i++) {
    if(at(i)->GetDecayProducts() == decayproducts) return at(i);
  }
  return NULL;
}

void Decay_Table::EraseDecayChannel(const int i) {
  delete at(i);
  for (size_t j=i;j<size()-1;j++) at(j) = at(j+1);
  pop_back();
}

Decay_Channel * Decay_Table::Select()
{
  Decay_Channel* selected(NULL);
  if (size()==1) selected=at(0);
  double disc = m_width*ran->Get();
  for (size_t i=0;i<size();++i) {
    disc -= at(i)->Width();
    if (disc<0) {
      selected=at(i);
      break;
    }
  }
  return selected;
}
