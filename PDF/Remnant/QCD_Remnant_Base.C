#include "PDF/Remnant/QCD_Remnant_Base.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"
#include <algorithm>
#include <iomanip>

#ifdef PROFILE__all
#define PROFILE__QCD_Remnant_Base
#endif
#ifdef PROFILE__QCD_Remnant_Base
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

using namespace PDF;

QCD_Remnant_Base::
QCD_Remnant_Base(PDF::ISR_Handler *isrhandler,
		 const unsigned int beam,const rtp::code type):
  Remnant_Base(type,beam), p_start(NULL), m_deltax(0.0125),
  m_xscheme(1), m_maxtrials(100), p_string(new double[2])
{
  m_scale=4.0;
  if (isrhandler==NULL) {
    THROW(fatal_error,"QCD remnant needs ISR Handler.");
  }
  p_pdfbase=isrhandler->PDF(m_beam)->GetBasicPDF();
}

QCD_Remnant_Base::~QCD_Remnant_Base()
{
  delete [] p_string;
}

void QCD_Remnant_Base::Clear()
{
  PROFILE_HERE;
  if (p_start!=NULL) delete p_start;
  while (m_connected.size()>0) {
    if (m_connected.front()!=p_start) delete m_connected.front();
    m_connected.erase(m_connected.begin());
  }
  Color_Dipole::ClearAll();
  Remnant_Base::Clear();
  p_start=NULL;
}

void QCD_Remnant_Base::AssignRemnants() 
{
  PROFILE_HERE;
  ATOOLS::Particle *startreal=p_start->Begin(qri::real);
  ATOOLS::Particle *startanti=p_start->Begin(qri::anti);
  for (ATOOLS::Particle_List::iterator pit=m_extracted.begin();
       pit!=m_extracted.end();++pit) {
    if (*pit==startreal || *pit==startanti) continue;
    if ((*pit)->Flav().Strong())
      m_connected.push_back(new Color_Dipole(*pit,&m_companions));
  }
}

Color_Dipole *QCD_Remnant_Base::FindClosest(const Color_Dipole *dipole,
					    const qri::type type)
{
  PROFILE_HERE;
  Color_Dipole *closest=p_start;
  const ATOOLS::Vec4D &ref=dipole->End(type)->Momentum();
  double min=std::numeric_limits<double>::max();
  int poss=0;
  for (Dipole_Vector::iterator dit=m_connected.begin();
       dit!=m_connected.end();++dit) {
    if ((*dit)->Next(ANTI(type))!=NULL ||
	(*dit)->Connected(dipole,ANTI(type))) continue;
    const ATOOLS::Vec4D &p=(*dit)->End(ANTI(type))->Momentum();
    double cur=p.PPerp(ref);
    if (p==ATOOLS::Vec4D()) cur=ref.PPerp();
    if (cur<=min) {
      min=cur;
      closest=*dit;
    }
    ++poss;
  }
  return closest;
}

Color_Dipole *QCD_Remnant_Base::FindRandom(const Color_Dipole *dipole,
					   const qri::type type)
{
  PROFILE_HERE;
  Dipole_Vector cand;
  for (Dipole_Vector::iterator dit=m_connected.begin();
       dit!=m_connected.end();++dit) {
    if ((*dit)->Next(ANTI(type))==NULL &&
	!(*dit)->Connected(dipole,ANTI(type))) cand.push_back(*dit);
  }
  if (cand.empty()) return p_start;
  double ran=ATOOLS::ran->Get()*cand.size();
  return cand[ATOOLS::Max((size_t)ran,cand.size()-1)];
}

Color_Dipole *QCD_Remnant_Base::Find(const Color_Dipole *dipole,
				     const qri::type type)
{
  if (p_string[1]==1.0) return FindRandom(dipole,type);
  return FindClosest(dipole,type);
}

class Compare_PT {
public:
  bool operator()(const std::pair<qri::type,Color_Dipole *> i1,
		  const std::pair<qri::type,Color_Dipole *> i2);
};

bool Compare_PT::operator()(const std::pair<qri::type,Color_Dipole *> i1,
			    const std::pair<qri::type,Color_Dipole *> i2)
{
  double pt21=i1.second->End(i1.first)->Momentum().PPerp2();
  double pt22=i2.second->End(i2.first)->Momentum().PPerp2();
  if (pt21==pt22) return i1.first<i2.first;
  return pt21>pt22;
}

bool QCD_Remnant_Base::Connect(const bool sort) 
{
  PROFILE_HERE;
  std::vector<std::pair<qri::type,Color_Dipole*> > sorted;
  for (Dipole_Vector::iterator dit=m_connected.begin();
       dit!=m_connected.end();++dit) {
    if ((*dit)->Begin(qri::real)==(*dit)->End(qri::real) ||
	(*dit)->End(qri::real)->ProductionBlob()!=p_beamblob) 
      sorted.push_back(std::pair<qri::type,Color_Dipole*>(qri::real,*dit));
    if ((*dit)->Begin(qri::anti)==(*dit)->End(qri::anti) ||
	(*dit)->End(qri::anti)->ProductionBlob()!=p_beamblob) 
      sorted.push_back(std::pair<qri::type,Color_Dipole*>(qri::anti,*dit));
  }
  std::stable_sort(sorted.begin(),sorted.end(),Compare_PT());
  sorted.push_back(std::pair<qri::type,Color_Dipole*>(qri::real,p_start));
  sorted.push_back(std::pair<qri::type,Color_Dipole*>(qri::anti,p_start));
  for (std::vector<std::pair<qri::type,Color_Dipole*> >::iterator
	 sit=sorted.begin();sit!=sorted.end();++sit) {
    if (sit->second->Next(sit->first)!=NULL) continue;
    Color_Dipole *next=sort?Find(sit->second,sit->first):
      FindRandom(sit->second,sit->first);
    if (!sit->second->Cat(next,sit->first)) {
      Color_Dipole::ResetAllColors();
      return false;
    }
  }
  return true;
}

bool QCD_Remnant_Base::AdjustColors() 
{
  PROFILE_HERE;
  if (!m_active) return true;
  QCD_Remnant_Base *partner=dynamic_cast<QCD_Remnant_Base*>(p_partner);
  if (partner==NULL) {
    for (size_t i(1);i<m_maxtrials;++i) {
      for (Dipole_Vector::iterator dit=m_connected.begin();
	   dit!=m_connected.end();++dit) {
	(*dit)->CollectStrings();
	(*dit)->DetectLoops();
      }
      p_start->CollectStrings();
      p_start->DetectLoops();
      if (Connect(i==1)) {
	Color_Dipole::SetAllColors();
	for (Dipole_Vector::iterator rit=m_connected.begin();
	     rit!=m_connected.end();++rit) {
	  if ((*rit)->Singlet(qri::real) || (*rit)->Singlet(qri::anti)) 
	    msg_Error()<<"QCD_Remnant_Base::AdjustColors(): "
		       <<"Colour singlet."<<std::endl;
	}
	return true;
      }
    }
    return false;
  }
  QCD_Remnant_Base *self=this;
  size_t i=0;
  while (++i<m_maxtrials) {
    bool order=i==1;
    for (short unsigned int j=0;j<2;++j) {
      if (j>0) std::swap<QCD_Remnant_Base*>(self,partner);
      for (Dipole_Vector::iterator dit=partner->m_connected.begin(); 
	   dit!=partner->m_connected.end();++dit) {
	(*dit)->CollectStrings();
	(*dit)->DetectLoops();
      }
      partner->p_start->CollectStrings();
      partner->p_start->DetectLoops();
      if (partner->Connect(order)) {
	Color_Dipole::SetAllColors();
	for (Dipole_Vector::iterator dit=self->m_connected.begin();
	     dit!=self->m_connected.end();++dit) {
	  (*dit)->CollectStrings();
	  (*dit)->DetectLoops();
	}
	self->p_start->CollectStrings();
	self->p_start->DetectLoops();
	if (self->Connect(order)) {
	  Color_Dipole::SetAllColors();
	  for (Dipole_Vector::iterator rit=partner->m_connected.begin();
	       rit!=partner->m_connected.end();++rit) {
	    if ((*rit)->Singlet(qri::real) || (*rit)->Singlet(qri::anti)) 
	      msg_Error()<<"QCD_Remnant_Base::AdjustColors(): "
				 <<"Colour singlet."<<std::endl;
	  }
	  for (Dipole_Vector::iterator rit=self->m_connected.begin();
	       rit!=self->m_connected.end();++rit) {
	    if ((*rit)->Singlet(qri::real) || (*rit)->Singlet(qri::anti)) 
	      msg_Error()<<"QCD_Remnant_Base::AdjustColors(): "
				 <<"Colour singlet."<<std::endl;
	  }
	  return true;
	}
      }
    }
  }
  msg_Tracking()<<"QCD_Remnant_Base::AdjustColors(): "
		<<"No solution in event ["
		<<ATOOLS::rpa->gen.NumberOfGeneratedEvents()<<"]."
		<<std::endl;
  return false;
}

void QCD_Remnant_Base::FillRemnants()
{
  PROFILE_HERE;
  for (ATOOLS::Particle_List::iterator pit=m_extracted.begin();
       pit!=m_extracted.end();++pit) {
    p_beamblob->AddToOutParticles(*pit);
  }
  for (ATOOLS::Particle_List::iterator pit=m_companions.begin();
       pit!=m_companions.end();++pit) {
    p_beamblob->AddToOutParticles(*pit);
  }
}
