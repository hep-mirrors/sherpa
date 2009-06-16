#include "PHASIC++/Process/NLO_Helpers.H"
#include "ATOOLS/Phys/Blob.H"

using namespace PHASIC;

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

ATOOLS::Particle_List* NLO_subevt::CreateParticleList() {
  ATOOLS::Particle_List * pl = new ATOOLS::Particle_List;
  for (size_t i=2;i<n;i++) {
    pl->push_back(new ATOOLS::Particle(i,p_fl[i],p_mom[i]));
  }
  if (m_flip) pl->Flip();
  return pl;
}

NLO_subevtlist& NLO_subevtlist::operator*=(const double scal) {
  for (std::vector<NLO_subevt*>::iterator it=begin();it!=end();it++) (*it)->m_result*=scal;
  return *this;
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
