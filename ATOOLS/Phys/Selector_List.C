#include "ATOOLS/Phys/Selector_List.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;
using namespace std;

Selector_Particle::~Selector_Particle() {}

Selector_List::Selector_List() : m_nin(0) {}

Selector_List::Selector_List(const Flavour_Vector &fl,
                             const Vec4D_Vector &p, const size_t &nin)
 : m_nin(nin)
{
  if (fl.size()!=p.size())
    THROW(fatal_error,"Number of flavours and momenta does not match.");
  for (size_t i(0);i<fl.size();++i) push_back(Selector_Particle(fl[i],p[i]));
}

Selector_List::Selector_List(const Flavour *fl, size_t n,
                             const Vec4D_Vector &p, const size_t &nin)
 : m_nin(nin)
{
  if (n!=p.size())
    THROW(fatal_error,"Number of flavours and momenta does not match.");
  for (size_t i(0);i<n;++i) push_back(Selector_Particle(fl[i],p[i]));
}

void Selector_List::AddFlavour(Flavour fl)
{
  push_back(Selector_Particle(fl,Vec4D(0.,0.,0.,0.)));
}

void Selector_List::SetMomenta(Vec4D_Vector p)
{
  if (size()!=p.size()) THROW(fatal_error,"Wrong number of momenta.");
  for (size_t i(0);i<size();++i) at(i).SetMomentum(p[i]);
}

Vec4D_Vector Selector_List::ExtractMomenta()
{
  Vec4D_Vector p;
  for (size_t i(0);i<size();++i) p.push_back(at(i).Momentum());
}

namespace ATOOLS {
  std::ostream &operator<<(std::ostream &ostr,const Selector_Particle &p)
  {
    ostr<<p.Flavour()<<": "<<p.Momentum();
    return ostr;
  }

  std::ostream &operator<<(std::ostream &ostr,const Selector_List &l)
  {
    ostr<<"Selector list:\n";
    for (size_t i(0);i<l.size();++i) {
      ostr<<"  "<<l[i]<<std::endl;
    }
    return ostr;
  }
}

