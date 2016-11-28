#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Permutation.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

#include <sstream>

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Single_Vertex::Single_Vertex(): order(2,0), dec(0)
{ 
}

Single_Vertex::~Single_Vertex()
{
}

bool Single_Vertex::PureQCD() const
{
  for (size_t i(1); i<order.size(); i++)
    if(order[i]) return false;
  return true;
}

int Single_Vertex::Compare(const Single_Vertex *v) const
{
  if (NLegs()!=v->NLegs()) return 1;
  if (cpl.size()!=v->cpl.size()) return 2;
  for (size_t i(0);i<cpl.size();++i)
    if (cpl[i].Value()!=v->cpl[i].Value()) return 2;
  for (size_t i(0);i<Lorentz.size();++i) {
    if (!(Color[i]==v->Color[i])) return 3;
    if (!(Lorentz[i]==v->Lorentz[i])) return 4;
  }
  return 0;
}

std::string Single_Vertex::PID() const
{
  std::string name(in[0].IDName());
  for (int i(1);i<NLegs();++i) name+='|'+in[i].IDName();
  return name;
}

bool Single_Vertex::operator==(const Single_Vertex& probe) 
{
  return (in==probe.in &&
	  cpl==probe.cpl &&
	  order==probe.order &&
	  Color==probe.Color &&
	  Lorentz==probe.Lorentz);
}

namespace MODEL{ 

  std::ostream &operator<<(std::ostream& s, const Single_Vertex& sv)
  {
    if (sv.in.size()) {
      s<<'('<<sv.in[0];
      for (size_t i(1);i<sv.in.size();++i) s<<','<<sv.in[i];
      s<<')';
    }
    if (sv.id.size()) {
      s<<'['<<sv.id[0];
      for (size_t i(1);i<sv.id.size();++i) s<<','<<sv.id[i];
      s<<']';
    }
    if (sv.cpl.size() && sv.Color.size() && sv.Lorentz.size()) {
      s<<"{{"<<sv.Coupling(0)<<"*"<<sv.Color[0].FullString()<<"*"<<sv.Lorentz[0];
      for (size_t i(1);i<sv.cpl.size();++i)
       	s<<"}{"<<sv.Coupling(i)<<"*"<<sv.Color[i].FullString()<<"*"<<sv.Lorentz[i];
      s<<"}}";
    }
    return s;
  }

}
