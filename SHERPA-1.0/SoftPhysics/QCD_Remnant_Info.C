#include "QCD_Remnant_Info.H"

#include "Color_Tester.H"
#include "Message.H"
#include <iomanip>

#ifdef PROFILE__all
#define PROFILE__QCD_Remnant_Info
#endif
#ifdef PROFILE__QCD_Remnant_Info
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

using namespace SHERPA;

std::ostream &SHERPA::operator<<(std::ostream &str,const qri::type type)
{
  switch (type) {
  case qri::real: return str<<"real";
  case qri::anti: return str<<"anti";
  default: break;
  }
  return str;
}

std::ostream &SHERPA::operator<<(std::ostream &str,
				 const QCD_Remnant_Info &info)
{
  str<<"QCD_Remnant_Info("<<&info<<"): ";
  for (short unsigned int i=0;i<2;++i) {
    qri::type type((qri::type)i);
    str<<"\n\n   p_begin["<<type<<"]    = "<<*info.p_begin[type]
       <<"\n   p_end["<<type<<"]      = "<<*info.p_end[type];
    size_t i=0;
    for (QCD_Remnant_Info::Particle_Flow_Map::const_iterator 
	   fit=info.m_flows[type].begin();
	 fit!=info.m_flows[type].end();++fit) {
      str<<"\n   m_flows["<<type<<"]["<<i++<<"] = "<<*fit->first<<" -> ("
	 <<std::setw(3)<<fit->second->Code(COLOR(qri::real))<<","
	 <<std::setw(3)<<fit->second->Code(COLOR(qri::anti))<<")";
    }
  }
  str<<"\n";
  size_t i=0;
  for (QCD_Remnant_Info::Particle_Flow_Map::const_iterator 
	 fit=info.m_treated.begin();fit!=info.m_treated.end();++fit) {
    str<<"\n   m_treated["<<i++<<"]     = "<<*fit->first<<" -> ("
       <<std::setw(3)<<fit->second->Code(COLOR(qri::real))<<","
       <<std::setw(3)<<fit->second->Code(COLOR(qri::anti))<<")";
  }
  return str<<"\n\n}"<<std::endl;
}

QCD_Remnant_Info::QCD_Remnant_Info(ATOOLS::Particle *const begin,
				   ATOOLS::Particle_List *const companions)
{
  SelectCompanion(begin,companions);
  CollectString(qri::real);
  CollectString(qri::anti);
}

QCD_Remnant_Info::~QCD_Remnant_Info()
{
  for (short unsigned int i=0;i<2;++i) {
    qri::type type((qri::type)i);
    while (m_flows[type].size()>0) {
      delete m_flows[type].begin()->second;
      m_flows[type].erase(m_flows[type].begin());
    }
    if (m_delete[type]) delete p_begin[type];
  }
}

void QCD_Remnant_Info::SelectCompanion(ATOOLS::Particle *const begin,
				       ATOOLS::Particle_List *const companions)
{
  qri::type anti((qri::type)(begin->Flav().IsAnti()^
			     begin->Flav().IsDiQuark()));
  m_delete[ANTI(anti)]=m_delete[anti]=false;
  p_begin[anti]=begin;
  if (p_begin[anti]->Flav().IsGluon()) {
    p_begin[ANTI(anti)]=begin;
    begin->GetFlow()->SetCode(COMPANIONTAG,0);
  }
  else {
    p_begin[ANTI(anti)] = 
      new ATOOLS::Particle(1,p_begin[anti]->Flav().Bar());
    ATOOLS::Flow *flow=p_begin[ANTI(anti)]->GetFlow();
    flow->SetCode(COLOR(ANTI(anti)),
		  p_begin[anti]->GetFlow()->Code(COLOR(anti)));
    flow->SetCode(COMPANIONTAG,-1);
    if (companions!=NULL) companions->push_back(p_begin[ANTI(anti)]);
    else m_delete[ANTI(anti)]=true;
  }
}

void QCD_Remnant_Info::CollectString(const qri::type type) 
{
  ATOOLS::Color_Tester tester(COLOR(type),p_begin[type]->GetFlow(COLOR(type)));
  ATOOLS::Parton_Finder finder(tester);
  p_end[type]=finder.FindConnected(p_begin[type],true);
  if (p_end[type]==NULL) p_end[type]=finder.FindConnected(p_begin[type],false);
  for (std::set<const ATOOLS::Particle *>::const_iterator 
	 pit=finder.Track().begin(); pit!=finder.Track().end();++pit) {
    m_flows[type][(ATOOLS::Particle*)*pit] = 
      new ATOOLS::Flow(*(*pit)->GetFlow());
  }
}

bool QCD_Remnant_Info::AssignColor(Particle_Flow_Map::iterator fit,
				   const unsigned int oldc,
				   const unsigned int newc)
{
  if (fit==m_flows[qri::real].end() ||
      fit==m_flows[qri::anti].end()) return true;
  int index=fit->second->Index(oldc);
  if (index<0) {
    ATOOLS::msg.Error()<<"QCD_Remnant_Info::AssignColor(..): "
		       <<"Incorrect color flow."<<std::endl//;
    //    msg_Tracking()
		       <<"   "<<*fit->first<<" => "
		       <<oldc<<" -> "<<newc<<std::endl;
    return false;
  }
  if (fit->second->Code(3-index)!=newc) {
    Particle_Flow_Map::iterator next=fit;
    if (AssignColor(++next,oldc,newc)) {
      fit->second->SetCode(index,newc);
      return true;
    }
  }
  return false;
}

bool QCD_Remnant_Info::AssignColors(const qri::type type,const int color)
{
  unsigned int oldc=p_begin[type]->GetFlow()->Code(COLOR(type));
  return m_status[type]=AssignColor(m_flows[type].begin(),oldc,color);
}

void QCD_Remnant_Info::SetColors()
{
  for (short unsigned int i=0;i<2;++i) {
    qri::type type((qri::type)i);
    for (Particle_Flow_Map::iterator fit=m_flows[type].begin();
	 fit!=m_flows[type].end();++fit) {
      fit->first->GetFlow()->SetCode(*fit->second);
    }
  }
  for (Particle_Flow_Map::iterator fit=m_treated.begin();
       fit!=m_treated.end();++fit) {
    fit->first->GetFlow()->SetCode(*fit->second);
  }
}

void QCD_Remnant_Info::Cat(Particle_Flow_Map::iterator fit)
{
  Particle_Flow_Map::iterator pos=m_treated.find(fit->first);
  if (pos==m_treated.end()) {
    m_treated[fit->first] = new ATOOLS::Flow(*fit->second);
  }
  else {
    ATOOLS::Flow *flow=pos->second;
    if (flow->Code(COLOR(qri::real))==
	fit->first->GetFlow()->Code(COLOR(qri::real))) 
      flow->SetCode(COLOR(qri::real),fit->second->Code(COLOR(qri::real)));
    else flow->SetCode(COLOR(qri::anti),fit->second->Code(COLOR(qri::anti)));
    flow->SetCode(COMPANIONTAG,fit->second->Code(0));
  }
  PRINT_INFO(*fit->first<<" "<<*m_treated[fit->first]);
}

bool QCD_Remnant_Info::Cat(QCD_Remnant_Info *const info,const qri::type type)
{
  qri::type anti=ANTI(type);

//   std::cout<<*info<<*this<<type<<ANTI(type)<<COLOR(anti)<<std::endl;

  AssignColors(type,info->Begin(anti)->GetFlow()->Code(COLOR(anti)));
  for (Particle_Flow_Map::iterator fit=info->m_flows[anti].begin();
       fit!=info->m_flows[anti].end();++fit) Cat(fit);
  for (Particle_Flow_Map::iterator fit=m_flows[type].begin();
       fit!=m_flows[type].end();++fit) Cat(fit);
  m_flows[type].clear();
  for (Particle_Flow_Map::iterator fit=info->m_flows[type].begin();
       fit!=info->m_flows[type].end();++fit) {
    m_flows[type][fit->first] = new ATOOLS::Flow(*fit->second);
  }
  p_begin[type]=info->p_begin[type];
  p_end[type]=info->p_end[type];
  delete info;

//   std::cout<<*this<<std::endl;
//   abort();

  return true;
}

