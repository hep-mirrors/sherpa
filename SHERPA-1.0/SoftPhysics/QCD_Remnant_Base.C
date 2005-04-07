#include "QCD_Remnant_Base.H"

#include "Exception.H"
#include "Random.H"
#include "Run_Parameter.H"
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

using namespace SHERPA;

QCD_Remnant_Base::QCD_Remnant_Base(PDF::ISR_Handler *isrhandler,
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
  if (p_start!=NULL) delete p_start;
  while (m_connected.size()>0) {
    if (m_connected.front()!=p_start) delete m_connected.front();
    m_connected.erase(m_connected.begin());
  }
  Color_Dipole::s_partons.clear();
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
    Color_Dipole *dipole = new Color_Dipole(*pit,&m_companions);
    if (p_start->Cat(dipole)) delete dipole;
    else {
      Dipole_Vector::iterator dit=m_connected.begin();
      for (;dit!=m_connected.end();++dit) 
	if ((*dit)->Cat(dipole)) {
	  delete dipole;
	  break;
	}
      if (dit==m_connected.end()) m_connected.push_back(dipole);
    }
  }
}

Color_Dipole *QCD_Remnant_Base::FindClosest(const Color_Dipole *dipole,
					    const qri::type type)
{
  Color_Dipole *closest=p_start;
  const ATOOLS::Vec4D &ref=dipole->End(ANTI(type))->Momentum();
  double min=std::numeric_limits<double>::max();
  std::multimap<double,Color_Dipole*> sorted;
  for (Dipole_Vector::iterator dit=m_attached.begin();
       dit!=m_attached.end();++dit) {
    if (*dit==dipole) continue;
    const ATOOLS::Vec4D &p=(*dit)->End(type)->Momentum();
    double cur=p.PPerp(ref);
    if (p==ATOOLS::Vec4D()) cur=ref.PPerp();
    if (p_string[0]!=1.0) 
      sorted.insert(std::pair<double,Color_Dipole*>(cur,*dit));
    if (cur<=min) {
      min=cur;
      closest=*dit;
    }
  }
  if (p_string[0]!=1.0) {
    double pos=(1.-p_string[0])*(sorted.size()-1), i=0.0;
    for (std::multimap<double,Color_Dipole*>::const_iterator dit=sorted.begin();
	 dit!=sorted.end();++dit) {
      //     std::cout<<pos<<" vs "<<i<<" "<<sorted.size()<<std::endl;
      if (i++>pos) return dit->second;
    }
  }
  return closest;
}

Color_Dipole *QCD_Remnant_Base::FindRandom(const Color_Dipole *dipole,
					   const qri::type type)
{
  double ran=ATOOLS::ran.Get()*(m_attached.size()-1), i=0.0;
  for (Dipole_Vector::iterator dit=m_attached.begin();
       dit!=m_attached.end();++dit) {
    if (*dit==dipole) continue;
    if (i++>ran) return *dit;
  }
  return p_start;
}

Color_Dipole *QCD_Remnant_Base::Find(const Color_Dipole *dipole,
				     const qri::type type)
{
  if (p_string[1]==1.0) return FindRandom(dipole,type);
  return FindClosest(dipole,type);
}

class Compare_PT {
public:
  bool operator()(const Color_Dipole *i1,const Color_Dipole *i2);
};

bool Compare_PT::operator()(const Color_Dipole *i1,const Color_Dipole *i2) 
{
  double pp21=ATOOLS::Max(i1->End(qri::real)->Momentum().PPerp2(),
			  i1->End(qri::anti)->Momentum().PPerp2());
  double pp22=ATOOLS::Max(i2->End(qri::real)->Momentum().PPerp2(),
			  i2->End(qri::anti)->Momentum().PPerp2());
  return (pp21<pp22);
}

bool QCD_Remnant_Base::Connect(const bool sorted) 
{
  m_attached.clear();
  m_attached.push_back(p_start);
  std::stable_sort(m_connected.begin(),m_connected.end(),Compare_PT());
  for (Dipole_Vector::iterator dit=m_connected.begin();
       dit!=m_connected.end();++dit) {
    qri::type type=(qri::type)((*dit)->End(qri::real)->Momentum().PPerp2()>=
			       (*dit)->End(qri::anti)->Momentum().PPerp2());
    if (!Find(*dit,type)->Insert(*dit,type)) {
      for (Dipole_Vector::iterator uit=m_connected.begin();
	   uit!=dit;++uit) (*uit)->UnDo();
      return false;
    }
    m_attached.push_back(*dit);
  }
  return true;
}

bool QCD_Remnant_Base::ConnectRemnants() 
{
  PROFILE_HERE;
  if (Connect(true)) return true;
  ATOOLS::msg.Error()<<"QCD_Remnant_Base::ConnectRemnants(): "
		     <<"No solution in event ["
		     <<ATOOLS::rpa.gen.NumberOfDicedEvents()<<"]. Abort."
		     <<std::endl;
  return false;
}

void QCD_Remnant_Base::SplitSinglet(Color_Dipole *const singlet)
{
  for (int i=0;i<2;++i) {
    qri::type type=(qri::type)i;
    if (!singlet->Singlet(type)) continue;
    msg_Debugging()<<"QCD_Remnant_Base::SplitSinglet(..): {\n"
		   <<*singlet<<"\n";
    m_attached.clear();
    Color_Dipole *dipole=Find(singlet,type);
    if (dipole!=NULL) {
      singlet->Cross(dipole,type);
    }
    else {
      for (size_t i=0;i<m_connected.size();++i) {
	dipole=m_connected[i];
	if (singlet->Cross(dipole,type)) break;
      }
    }
    if (dipole==NULL) return;
    dipole->SetColors();
    singlet->SetColors();
    msg_Debugging()<<*singlet<<*dipole<<"}"<<std::endl;
  }
}

void QCD_Remnant_Base::FillRemnants()
{
  PROFILE_HERE;
  for (Dipole_Vector::iterator rit=m_connected.begin();
       rit!=m_connected.end();++rit) {
    (*rit)->SetColors();
  }
  p_start->SetColors();
  for (size_t i=0;i<m_connected.size();++i) SplitSinglet(m_connected[i]);
  for (ATOOLS::Particle_List::iterator pit=m_extracted.begin();
       pit!=m_extracted.end();++pit) {
    p_beamblob->AddToOutParticles(*pit);
  }
  for (ATOOLS::Particle_List::iterator pit=m_companions.begin();
       pit!=m_companions.end();++pit) {
    p_beamblob->AddToOutParticles(*pit);
  }
}

