#include "Color_Dipole.H"

#include "Color_Tester.H"
#include "Message.H"
#include <iomanip>

#ifdef PROFILE__all
#define PROFILE__Color_Dipole
#endif
#ifdef PROFILE__Color_Dipole
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

#define COMPANIONTAG 0

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
				 const Color_Dipole &info)
{
  str<<"Color_Dipole("<<&info<<"): {";
  for (short unsigned int i=0;i<2;++i) {
    qri::type type((qri::type)i);
    str<<"\n\n   p_begin["<<type<<"]    = ";
    if (info.p_begin[type]!=NULL) str<<*info.p_begin[type]; 
    else str<<"NULL";
    str<<"\n   p_end["<<type<<"]      = ";
    if (info.p_end[type]!=NULL) str<<*info.p_end[type];
    else str<<"NULL";
    size_t i=0;
    for (Color_Dipole::Particle_Flow_Map::const_iterator 
	   fit=info.m_flows[type].begin();
	 fit!=info.m_flows[type].end();++fit) {
      str<<"\n   m_flows["<<type<<"]["<<i++<<"] = "<<*fit->first<<" -> ("
	 <<std::setw(3)<<fit->second->Code(COLOR(qri::real))<<","
	 <<std::setw(3)<<fit->second->Code(COLOR(qri::anti))<<")";
    }
  }
  return str<<"\n\n}"<<std::endl;
}

Color_Dipole::Color_Dipole():
  p_companions(NULL)
{
  p_begin[qri::real]=p_end[qri::real]=NULL;
  p_begin[qri::anti]=p_end[qri::anti]=NULL;
}

Color_Dipole::Color_Dipole(ATOOLS::Particle *const begin,
			   ATOOLS::Particle_List *const companions):
  p_companions(companions)
{
  SelectCompanion(begin);
  CollectString(qri::real);
  CollectString(qri::anti);
//   msg_Tracking()<<*this<<std::endl;
}

Color_Dipole::~Color_Dipole()
{
  Clear(qri::real);
  Clear(qri::anti);
}

void Color_Dipole::Clear(const qri::type type)
{
  while (m_flows[type].size()>0) {
    delete m_flows[type].begin()->second;
    m_flows[type].erase(m_flows[type].begin());
  }
  p_end[type]=p_begin[type]=NULL;
}

void Color_Dipole::SelectCompanion(ATOOLS::Particle *const begin)
{
  qri::type anti((qri::type)(begin->Flav().IsAnti()^
			     begin->Flav().IsDiQuark()));
  p_begin[anti]=begin;
  if (p_begin[anti]->Flav().IsGluon()) {
    p_begin[ANTI(anti)]=begin;
    begin->GetFlow()->SetCode(COMPANIONTAG,0);
  }
  else {
    if (p_companions==NULL) {
      p_end[ANTI(anti)]=p_begin[ANTI(anti)]=NULL;
      return;
    }
    p_begin[ANTI(anti)] = 
      new ATOOLS::Particle(-1,p_begin[anti]->Flav().Bar());
    p_begin[ANTI(anti)]->SetNumber(0);
    ATOOLS::Flow *flow=p_begin[ANTI(anti)]->GetFlow();
    flow->SetCode(COLOR(ANTI(anti)),
		  p_begin[anti]->GetFlow()->Code(COLOR(anti)));
    flow->SetCode(COMPANIONTAG,1);
    p_companions->push_back(p_begin[ANTI(anti)]);
  }
}

void Color_Dipole::CollectString(const qri::type type) 
{
  if (p_begin[type]==NULL) return;
  ATOOLS::Color_Tester tester(COLOR(type),p_begin[type]->GetFlow(COLOR(type)));
  ATOOLS::Parton_Finder finder(tester);
  p_end[type]=finder.FindConnected(p_begin[type],true);
  if (p_end[type]==NULL) p_end[type]=finder.FindConnected(p_begin[type],false);
  for (std::vector<const ATOOLS::Particle *>::const_iterator 
	 pit=finder.Track().begin(); pit!=finder.Track().end();++pit) {
    m_flows[type][(ATOOLS::Particle*)*pit] = 
      new ATOOLS::Flow(*(*pit)->GetFlow());
  }
}

void Color_Dipole::Prepend(ATOOLS::Particle *const part,
			   const qri::type type,const bool same)
{
  unsigned int color=Begin(type)->GetFlow()->Code(COLOR(type));
  unsigned int index=COLOR(type);
  if (!same) index=3-index;
  part->GetFlow()->SetCode(index,color);
  m_flows[type][part] = new ATOOLS::Flow(*part->GetFlow());
  p_begin[type]=part;
}

void Color_Dipole::Append(ATOOLS::Particle *const part,
			  const qri::type type,const bool same)
{
  unsigned int color=Begin(type)->GetFlow()->Code(COLOR(type));
  unsigned int index=End(type)->GetFlow()->Index(color);
  if (!same) index=3-index;
  part->GetFlow()->SetCode(index,color);
  m_flows[type][part] = new ATOOLS::Flow(*part->GetFlow());
  p_end[type]=part;
}

bool Color_Dipole::Cat(Color_Dipole *const dipole,const qri::type type)
{
  unsigned int color=p_begin[type]->GetFlow()->Code(COLOR(type));
//   msg_Debugging()<<"Color_Dipole::Cat(..): ["<<type<<"] ("<<color<<") {"
// 		 <<"\n   t->b = "<<*p_begin[type]
// 		 <<"\n   d->e = "<<*dipole->p_end[ANTI(type)]
// 		 <<"\n   t->e = "<<*p_end[type]
// 		 <<"\n   d->b = "<<*dipole->p_begin[ANTI(type)]
// 		 <<"\n}"<<std::endl;
  dipole->AssignColors(ANTI(type),ATOOLS::Flow::Counter());
  dipole->SetColors(ANTI(type));
  Clear(type);
  p_begin[type]=dipole->p_begin[type];
  p_end[type]=dipole->p_end[type];
  for (Particle_Flow_Map::iterator fit=dipole->m_flows[type].begin();
       fit!=dipole->m_flows[type].end();++fit) 
    m_flows[type][fit->first] = new ATOOLS::Flow(*fit->first->GetFlow());
  AssignColors(type,color);
  SetColors(type);
//   msg_Debugging()<<*this<<std::endl;
  return true;
}

bool Color_Dipole::Cat(Color_Dipole *const dipole)
{
  bool success=false;
  for (short unsigned int i=0;i<2;++i) {
    if (Begin((qri::type)i)==dipole->End(ANTI(i)))
      if (Cat(dipole,(qri::type)i)) success=true;
  }
  return success;
}

bool Color_Dipole::Includes(ATOOLS::Particle *const part)
{
  return (m_flows[qri::real].find(part)!=m_flows[qri::real].end() ||
	  m_flows[qri::anti].find(part)!=m_flows[qri::anti].end());
}

bool Color_Dipole::AssignColor(Particle_Flow_Map::iterator fit,
				   const unsigned int oldc,
				   const unsigned int newc)
{
  if (fit==m_flows[qri::real].end() ||
      fit==m_flows[qri::anti].end()) return true;
  int index=fit->second->Index(oldc);
  if (index<0) {
    ATOOLS::msg.Error()<<"Color_Dipole::AssignColor(..): "
		       <<"Invalid color {\n   "<<*fit->second
		       <<" => ("<<oldc<<" -> "<<newc<<")\n   "
		       <<*fit->first<<"\n}"<<std::endl;
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

bool Color_Dipole::AssignColors(const qri::type type,const int color)
{
  unsigned int oldc=m_flows[type][p_begin[type]]->Code(COLOR(type));
  return m_status[type]=AssignColor(m_flows[type].begin(),oldc,color);
}

void Color_Dipole::UnDo(const qri::type type)
{
  for (Particle_Flow_Map::iterator fit=m_flows[type].begin();
       fit!=m_flows[type].end();++fit) {
    fit->second->SetCode(*fit->first->GetFlow());
  }
}

void Color_Dipole::UnDo()
{
  for (short unsigned int i=0;i<2;++i) UnDo((qri::type)i);
}

void Color_Dipole::SetColors(const qri::type type)
{
  unsigned int color=p_begin[type]->GetFlow()->Code(COLOR(type));
  for (Particle_Flow_Map::iterator fit=m_flows[type].begin();
       fit!=m_flows[type].end();++fit) {
    unsigned int index=fit->first->GetFlow()->Index(color);
    fit->first->GetFlow()->SetCode(index,fit->second->Code(index));
  }
}

void Color_Dipole::SetColors()
{
  for (short unsigned int i=0;i<2;++i) SetColors((qri::type)i);
}

void Color_Dipole::Split(const qri::type type)
{
  ATOOLS::Particle *gluon = new ATOOLS::Particle(-1,ATOOLS::kf::gluon);
  gluon->SetNumber(0);
  unsigned int color=ATOOLS::Flow::Counter();
  gluon->GetFlow()->SetCode(COLOR(type),p_begin[ANTI(type)]->
			    GetFlow()->Code(COLOR(ANTI(type))));
  gluon->GetFlow()->SetCode(COLOR(ANTI(type)),color);
  AssignColors(type,color);
  SetColors(type);
  Clear(type);
  p_companions->push_back(gluon);
  p_end[type]=p_begin[type]=gluon;
  m_flows[type][gluon] = new ATOOLS::Flow(*gluon->GetFlow());
//   msg_Debugging()<<"Color_Dipole::Split("<<type<<"): "
//  		 <<*this<<std::endl;
}

bool Color_Dipole::Insert(Color_Dipole *const info,const qri::type type)
{
  qri::type anti=ANTI(type);
  int oldc=m_flows[type][Begin(type)]->Code(COLOR(type));
  int newc=info->m_flows[anti][info->Begin(anti)]->Code(COLOR(anti));
  if (!AssignColors(type,newc)) {
    Split(type);
    AssignColors(type,newc);
  }
  if (!info->AssignColors(type,oldc)) {
    info->Split(type);
    info->AssignColors(type,oldc);
  }
//   msg_Debugging()<<"Color_Dipole::Insert(..): ["<<type<<"] => ("
// 		 <<oldc<<","<<newc<<") -> ("<<") {\n   /"
// 		 <<*info->Begin(type)
// 		 <<" -> ("<<std::setw(3)
// 		 <<info->m_flows[type][info->Begin(type)]->Code(1)
// 		 <<","<<std::setw(3)
// 		 <<info->m_flows[type][info->Begin(type)]->Code(2)<<")"
// 		 <<"\n   \\"<<*info->Begin(anti)
// 		 <<" -> ("<<std::setw(3)
// 		 <<info->m_flows[anti][info->Begin(anti)]->Code(1)
// 		 <<","<<std::setw(3)
// 		 <<info->m_flows[anti][info->Begin(anti)]->Code(2)<<")"
// 		 <<"\n   /"<<*Begin(type)
// 		 <<" -> ("<<std::setw(3)<<m_flows[type][Begin(type)]->Code(1)
// 		 <<","<<std::setw(3)<<m_flows[type][Begin(type)]->Code(2)<<")"
// 		 <<"\n   \\"<<*Begin(anti)
// 		 <<" -> ("<<std::setw(3)<<m_flows[anti][Begin(anti)]->Code(1)
// 		 <<","<<std::setw(3)<<m_flows[anti][Begin(anti)]->Code(2)<<")"
// 		 <<"\n}"<<std::endl;
  return true;
}

