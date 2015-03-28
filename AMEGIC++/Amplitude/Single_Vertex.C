#include "AMEGIC++/Amplitude/Single_Vertex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

using namespace ATOOLS;
using namespace AMEGIC;

// Constructors and Destructors
Single_Vertex::Single_Vertex(): order(2,0)
{ t = 0; nleg=3; cpl.resize(4); order[1]=1; dec=0; }

Single_Vertex::Single_Vertex(const Single_Vertex& v): 
  t(0)
{ 
  *this=v; 
}

Single_Vertex::~Single_Vertex()
{
  for (size_t i(0);i<Lorentz.size();++i) delete Lorentz[i];
}

int Single_Vertex::Compare(const Single_Vertex *v) const
{
  if (nleg!=v->nleg) return 1;
  if (cpl.size()!=v->cpl.size()) return 2;
  for (size_t i(0);i<cpl.size();++i)
    if (cpl[i].Value()!=v->cpl[i].Value()) return 2;
  for (size_t i(0);i<Lorentz.size();++i) {
    if (!(Color[i]==v->Color[i])) return 3;
    if (!(*Lorentz[i]==*v->Lorentz[i])) return 4;
  }
  return 0;
}

Complex Single_Vertex::Coupling(size_t i) const
{
  return cpl[i].Value();
}
 
std::string Single_Vertex::PID() const
{
  std::string name;
  for (int i(1);i<nleg;++i) name+='{'+in[i].IDName()+'}';
  return name+'{'+in[0].IDName()+'}';
}

int Single_Vertex::CheckCoupling() const
{
  for (size_t i(0);i<cpl.size();++i) 
    if (cpl[i].Value()!=Complex(0.,0.)) return 1;
  return 0;
}


// Operators
Single_Vertex& Single_Vertex::operator=(const Single_Vertex& v) 
{
  for (size_t i(0);i<Lorentz.size();++i) delete Lorentz[i];
  Lorentz=std::vector<MODEL::Lorentz_Function*>();
    
      if (this!=&v) {
	for (short int i=0;i<4;i++) in[i]  = v.in[i];
	cpl.clear();
	for (size_t j=0;j<v.cpl.size();j++) cpl.push_back(v.cpl[j]);
	
	nleg = v.nleg;
	on   = v.on;
	t=v.t;
	order=v.order;
	dec=v.dec;
	Color.resize(v.Color.size());
	for (size_t i(0);i<Color.size();++i)
	  Color[i] = v.Color[i];
	Lorentz.resize(v.Lorentz.size());
	for (size_t i(0);i<Lorentz.size();++i)
	  Lorentz[i] = v.Lorentz[i]->GetCopy();
      }
      return *this;
    }


 
bool Single_Vertex::operator==(const Single_Vertex& probe) 
{
  switch (nleg) // different checks for 3-leg and 4-leg vertices
  {
    case 4: return (probe.nleg == 4) &&
       	           (in[0] == probe.in[0]) &&
                   (in[1] == probe.in[1]) &&
                   (in[2] == probe.in[2]) &&
                   (in[3] == probe.in[3]);
    case 3:  return (probe.nleg == 3) &&
	            (in[0] == probe.in[0]) &&
                    (in[1] == probe.in[1]) &&
                    (in[2] == probe.in[2]);
    default: return 0; 
  }
}

namespace AMEGIC{ 
std::ostream &operator<<(std::ostream& s, const Single_Vertex& sv)
{
  return s<<'('<<sv.in[0]<<','<<sv.in[1]<<','<<sv.in[2]<<','<<sv.in[3]
          <<") with cpl["<<sv.Coupling(0)<<','<<sv.Coupling(1)<<','<<sv.Coupling(2)<<','<<sv.Coupling(3)<<']'
	  <<" is "<<((sv.on) ? "on" : "off");
}
}
