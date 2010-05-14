#include "PHASIC++/Process/NLO_Helpers.H"
#include "ATOOLS/Phys/Blob.H"

using namespace PHASIC;

NLO_subevt::NLO_subevt() : n(0), p_fl(0), p_mom(0),
			   m_result(0.), m_me(0.), m_mewgt(0.),
			   m_facscale(0.), m_renscale(0.), m_alpha(0.),
			   m_ID("")
{}

NLO_subevt::NLO_subevt(size_t nn,const ATOOLS::Flavour* fl,const ATOOLS::Vec4D* mom,std::string ID) : 
  n(nn), p_fl(fl), p_mom(mom),
  m_result(0.), m_me(0.), m_mewgt(0.),
  m_facscale(0.), m_renscale(0.), m_alpha(0.),
  m_ID(ID)
{ }

NLO_subevt::NLO_subevt(NLO_subevt* nlos)
{ 
  *this = *nlos; 
  m_result=m_me; 
}

void NLO_subevt::Print() 
{
  std::cout<<m_ID<<": "<<m_result<<std::endl;
}
  
ATOOLS::Particle_List* NLO_subevt::CreateParticleList() {
  ATOOLS::Particle_List * pl = new ATOOLS::Particle_List;
  for (size_t i=2;i<n;i++) {
    pl->push_back(new ATOOLS::Particle(i,p_fl[i],p_mom[i]));
  }
  if (m_flip) pl->Flip();
  return pl;
}

void NLO_subevt::Mult(const double scal) 
{
  m_result*=scal;
  m_me*=scal;
  m_mewgt*=scal;
}
    
void NLO_subevt::MultMEwgt(const double scal) {
  m_mewgt*=scal;
}

NLO_subevt& NLO_subevt::operator*= (const double scal) {
  m_result*=scal;
  return *this;
}

NLO_subevt& NLO_subevt::operator+= (const double sm) {
  m_result+=sm;
  return *this;
}

NLO_subevt& NLO_subevt::operator=(const NLO_subevt& cp) {
  if (this!=&cp) {
    n        = cp.n;
    p_fl     = cp.p_fl;
    p_mom    = cp.p_mom;
    m_result = cp.m_result;
    m_mewgt  = cp.m_mewgt;
    m_me     = cp.m_me;
    m_facscale  = cp.m_facscale;
    m_renscale  = cp.m_renscale;
    m_alpha  = cp.m_alpha;
    m_ID     = cp.m_ID;
  }
  return *this;
}

void NLO_subevtlist::Print()
{
  for (std::vector<NLO_subevt*>::iterator it=begin();it!=end();it++){
    (*it)->Print();
  }
}

void NLO_subevtlist::Mult(const double scal)
{
  for (std::vector<NLO_subevt*>::iterator it=begin();it!=end();it++){
    (*it)->Mult(scal);
  }
}

NLO_subevtlist& NLO_subevtlist::operator*=(const double scal) {
  for (std::vector<NLO_subevt*>::iterator it=begin();it!=end();it++) (*it)->m_result*=scal;
  return *this;
}

void NLO_subevtlist::MultMEwgt(const double scal)
{
  for (std::vector<NLO_subevt*>::iterator it=begin();it!=end();it++){
    (*it)->MultMEwgt(scal);
  }
}


namespace ATOOLS {
  template<PHASIC::NLO_subevt*>
  std::ostream& operator<<(std::ostream& s, const PHASIC::NLO_subevtlist &sel) {
    s<<sel.size()<<"\n";
    return s;
  }
}
namespace ATOOLS {
  template <> Blob_Data<PHASIC::NLO_subevtlist*>::~Blob_Data() {}
  template class Blob_Data<PHASIC::NLO_subevtlist*>;
}

std::ostream &PHASIC::operator<<(std::ostream &str,const nlo_type::code &c) 
{
  std::string out="";
  if (c&nlo_type::born) out+="B";
  if (c&nlo_type::loop) out+="V";
  if (c&nlo_type::vsub) out+="I";
  if (c&nlo_type::real) out+="R";
  if (c&nlo_type::rsub) out+="S";
  if (c&nlo_type::polecheck) out+="P";
  return str<<out;
}

std::istream &PHASIC::operator>>(std::istream &str,nlo_type::code &c) 
{
  std::string tag;
  str>>tag;
  c=nlo_type::lo;
  if (tag.find('B')!=std::string::npos) c|=nlo_type::born;
  if (tag.find('V')!=std::string::npos) c|=nlo_type::loop;
  if (tag.find('I')!=std::string::npos) c|=nlo_type::vsub;
  if (tag.find('R')!=std::string::npos) c|=nlo_type::real;
  if (tag.find('S')!=std::string::npos) c|=nlo_type::rsub;
  if (tag.find('P')!=std::string::npos) c|=nlo_type::polecheck;
  return str;
}
