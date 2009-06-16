#include "PDF/Remnant/Color_Dipole.H"

#include "ATOOLS/Phys/Color_Tester.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
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

using namespace PDF;

Color_Dipole::Dipole_Vector Color_Dipole::s_dipoles;
Color_Dipole::Particle_Flow_Map Color_Dipole::s_flows;
Color_Dipole::Particle_Flow_Map Color_Dipole::s_oldflows;

std::set<ATOOLS::Particle*> PDF::Color_Dipole::s_partons;

std::ostream &PDF::operator<<(std::ostream &str,const qri::type type)
{
  switch (type) {
  case qri::real: return str<<"real";
  case qri::anti: return str<<"anti";
  default: break;
  }
  return str;
}

std::ostream &PDF::operator<<(std::ostream &str,
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
    size_t j=0;
    for (Color_Dipole::Particle_Flow_Map::const_iterator 
	   fit=info.m_flows[type].begin();
	 fit!=info.m_flows[type].end();++fit) {
      str<<"\n   m_flows["<<type<<"]["<<j++<<"] = "<<*fit->first<<" -> ("
	 <<std::setw(3)<<fit->second->Code(COLOR(qri::real))<<","
	 <<std::setw(3)<<fit->second->Code(COLOR(qri::anti))<<")";
    }
  }
  return str<<"\n\n}"<<std::endl;
}

Color_Dipole::Color_Dipole():
  p_companions(NULL)
{
  s_dipoles.push_back(this);
  p_begin[qri::anti]=p_begin[qri::real]=NULL;
  p_end[qri::anti]=p_end[qri::real]=NULL;
  p_next[qri::anti]=p_next[qri::real]=NULL;
}

Color_Dipole::Color_Dipole(ATOOLS::Particle *const begin,
			   ATOOLS::Particle_List *const companions):
  p_companions(companions)
{
  s_dipoles.push_back(this);
  SelectCompanion(begin);
  p_next[qri::anti]=p_next[qri::real]=NULL;
}

Color_Dipole::~Color_Dipole()
{
  for (Dipole_Vector::iterator dit=s_dipoles.begin();
       dit!=s_dipoles.end();++dit) 
    if (*dit==this) {
      s_dipoles.erase(dit);
      break;
    }
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
    p_begin[ANTI(anti)]->SetStatus(ATOOLS::part_status::active);
    p_begin[ANTI(anti)]->SetFinalMass(p_begin[ANTI(anti)]->Flav().Mass());
    p_begin[ANTI(anti)]->SetNumber(0);
    p_begin[ANTI(anti)]->SetInfo('F');
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
  ATOOLS::Color_Tester tester(COLOR(type),p_begin[type]->
			      GetFlow(COLOR(type)));
  ATOOLS::Parton_Finder finder(tester);
  p_end[type]=finder.FindConnected(p_begin[type],true);
  if (p_end[type]==NULL) p_end[type]=
    finder.FindConnected(p_begin[type],false);
  m_flows[type].clear();
  p_next[type]=NULL;
  for (ATOOLS::Particle_List::const_iterator 
	 pit=finder.Track().begin(); pit!=finder.Track().end();++pit) {
    if (s_flows.find(*pit)==s_flows.end()) {
      s_flows[*pit] = new ATOOLS::Flow(*(*pit)->GetFlow());
      s_oldflows[*pit] = new ATOOLS::Flow(*(*pit)->GetFlow());
    }
    m_flows[type][*pit]=s_flows[*pit];
    if ((*pit)->DecayBlob()==NULL &&
	(*pit)->ProductionBlob()==NULL) s_partons.insert(*pit);
  }
}

void Color_Dipole::CollectStrings() 
{
  CollectString(qri::real);
  CollectString(qri::anti);
}

void Color_Dipole::DetectLoop(const qri::type type)
{
  if (p_begin[type]!=p_end[type] &&
      p_begin[type]->ProductionBlob()!=NULL &&
      p_begin[type]->ProductionBlob()==
      p_end[type]->ProductionBlob()) {
    p_next[type]=this;
  }
}

void Color_Dipole::DetectLoops()
{
  for (short unsigned int i=0; i<2;++i) DetectLoop((qri::type)i);
}

bool Color_Dipole::AssignColor(Particle_Flow_Map::iterator fit,
			       const unsigned int oldc,
			       const unsigned int newc)
{
  if (fit==m_flows[qri::real].end() ||
      fit==m_flows[qri::anti].end()) return true;
  int index=fit->second->Index(oldc);
  if (index<0) {
    msg_Error()<<"Color_Dipole::AssignColor(..): "
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
  return AssignColor(m_flows[type].begin(),oldc,color);
}

bool Color_Dipole::Singlet(const qri::type type) const 
{
  for (Particle_Flow_Map::const_iterator fit=m_flows[type].begin();
       fit!=m_flows[type].end();++fit) {
    if (fit->first->GetFlow(COLOR(qri::real))==
	fit->first->GetFlow(COLOR(qri::anti))) return true;
  }
  return false;
}

bool Color_Dipole::Cat(Color_Dipole *const dipole,const qri::type type)
{
  unsigned int color=ATOOLS::Flow::Counter();
  dipole->AssignColors(ANTI(type),color);
  if (!AssignColors(type,color)) return false;
  dipole->p_next[ANTI(type)]=this;
  p_next[type]=dipole;
  return true;
}

bool Color_Dipole::
Connected(const Color_Dipole *dipole,
	  const qri::type type,const size_t catcher) const
{
  if (catcher>1000) 
    THROW(fatal_error,"Dipole nesting deeper than 1000 levels.");
  if (this==dipole) return true;
  if (p_next[ANTI(type)]!=NULL && p_next[ANTI(type)]!=this) 
    return p_next[ANTI(type)]->Connected(dipole,type,catcher+1);
  return false;
}

void Color_Dipole::SetAllColors()
{
  for (Particle_Flow_Map::iterator fit=s_flows.begin();
       fit!=s_flows.end();++fit) {
    ATOOLS::Flow *flow=fit->first->GetFlow();
    for (short unsigned int i=0;i<2;++i) 
      flow->SetCode(COLOR((qri::type)i),
		    fit->second->Code(COLOR((qri::type)i)));
  }
}

void Color_Dipole::ResetAllColors()
{
  for (Particle_Flow_Map::iterator fit=s_oldflows.begin();
       fit!=s_oldflows.end();++fit) {
    ATOOLS::Flow *oflow=fit->first->GetFlow(), *nflow=s_flows[fit->first];
    for (short unsigned int i=0;i<2;++i) {
      oflow->SetCode(COLOR((qri::type)i),
		     fit->second->Code(COLOR((qri::type)i)));
      nflow->SetCode(COLOR((qri::type)i),
		     fit->second->Code(COLOR((qri::type)i)));
    }
  }
}

void Color_Dipole::ClearAll()
{
  while (s_flows.size()>0) {
    delete s_flows.begin()->second;
    s_flows.erase(s_flows.begin());
  }
  while (s_oldflows.size()>0) {
    delete s_oldflows.begin()->second;
    s_oldflows.erase(s_oldflows.begin());
  }
  s_partons.clear();
}

