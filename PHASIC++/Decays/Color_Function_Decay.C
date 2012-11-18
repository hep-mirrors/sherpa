#include "PHASIC++/Decays/Color_Function_Decay.H"

#include "MODEL/Interaction_Models/Color_Function.H"
#include "ATOOLS/Phys/Color.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <algorithm>

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

size_t Color_Function_Decay::m_nid(100);

Color_Function_Decay::Color_Function_Decay() : m_max(-1)
{
  push_back(make_pair("1",vector<int>()));
}

Color_Function_Decay::Color_Function_Decay(const MODEL::Color_Function& c1) :
  m_max(-1)
{
  push_back(make_pair(string(""),vector<int>()));
  std::string& name = back().first;
  std::vector<int>& inds = back().second;
  switch(c1.Type()) {
  case MODEL::cf::T:
    name="T";
    for (size_t i(0); i<3; ++i)
      inds.push_back(c1.ParticleArg(i));
    break;
  case MODEL::cf::F:
    name="F";
    for (size_t i(0); i<3; ++i)
      inds.push_back(c1.ParticleArg(i));
    break;
  case MODEL::cf::D:
    name="D";
    for (size_t i(0); i<2; ++i)
      inds.push_back(c1.ParticleArg(i));
    break;
  case MODEL::cf::G:
    name="G";
    for (size_t i(0); i<2; ++i)
      inds.push_back(c1.ParticleArg(i));
    break;
  case MODEL::cf::None:
    name="1";
    break;
  default:
    THROW(not_implemented,"Encountered unknown type of color function: "+
          ToString(c1.Type()));
  }
  for (size_t i(0); i<back().second.size(); ++i) {
    if (back().second[i]>m_max) m_max=back().second[i];
  }
}

void Color_Function_Decay::Conjugate()
{
  for (size_t i=0; i<size(); ++i) {
    if (at(i).first=="T") {
      int help=at(i).second[1];
      at(i).second[1]=at(i).second[2];
      at(i).second[2]=help;
    }
    else if (at(i).first=="D") {
      int help=at(i).second[0];
      at(i).second[0]=at(i).second[1];
      at(i).second[1]=help;
    }
    else if (at(i).first=="G") {
      int help=at(i).second[0];
      at(i).second[0]=at(i).second[1];
      at(i).second[1]=help;
    }
  }
}

std::string Color_Function_Decay::String() const
{
  std::string ret;
  for (size_t i(0); i<size(); ++i) {
    if (i>0) ret+="*";
    if (at(i).first=="G") {
      ret+="2*T["+ToString(at(i).second[0])+","
	+ToString(m_nid)+","+ToString(m_nid+1)+"]"
	+"*T["+ToString(at(i).second[1])+","
	+ToString(m_nid+1)+","+ToString(m_nid)+"]";
      m_nid+=2;
    }
    else {
    ret+=at(i).first;
    for (size_t j(0); j<at(i).second.size(); ++j) {
      if (j==0) ret+="["+ToString(at(i).second[j]);
      else ret+=","+ToString(at(i).second[j]);

      if (j==at(i).second.size()-1) ret+="]";
    }
    }
  }
  return ret;
}

void Color_Function_Decay::BumpIndices(size_t i, int bump)
{
  for (size_t j(0); j<at(i).second.size(); ++j) at(i).second[j]+=bump;
  for (size_t j(0); j<at(i).second.size(); ++j) {
    if (at(i).second[j]>m_max) m_max=at(i).second[j];
  }
}

std::vector<int> Color_Function_Decay::Multiply(const Color_Function_Decay& c)
{
  DEBUG_FUNC(String()<<" * "<<c.String());
  // returns a vector of what has been added to the indices of each color
  // function in c to make the indices unique.
  std::vector<int> ret(c.size());
  for (size_t i(0); i<c.size(); ++i) {
    push_back(make_pair(c[i].first, c[i].second));
    // make indices unique
    ret[i]=m_max+1;
    BumpIndices(size()-1, ret[i]);
  }
  return ret;
}

void Color_Function_Decay::Contract(int i1, int i2)
{
  DEBUG_FUNC(i1<<" "<<i2);
  int has_i1(HasIndex(i1)), has_i2(HasIndex(i2));
  if (has_i1==0 || has_i2==0) return;
  else if (has_i1==1 && has_i2==1) {
    for (size_t i(0); i<size(); ++i) {
      for (size_t j(0); j<at(i).second.size(); ++j) {
        if (at(i).second[j]==i2) at(i).second[j]=i1;
      }
    }
    m_internal.push_back(i1);
  }
  else
    THROW(fatal_error, "Internal error.");
}

void Color_Function_Decay::ReplaceIndex(int i1, int i2)
{
  for (size_t i(0); i<size(); ++i) {
    for (size_t j(0); j<at(i).second.size(); ++j) {
      if (at(i).second[j]==i1) at(i).second[j]=i2;
    }
  }
}

int Color_Function_Decay::HasIndex(int i1) const
{
  int ret(0);
  for (size_t i(0); i<size(); ++i) {
    for (size_t j(0); j<at(i).second.size(); ++j) {
      if (at(i).second[j]==i1) ++ret;
    }
  }
  return ret;
}

Complex Color_Function_Decay::Contract(const Color_Function_Decay& c)
{
  DEBUG_FUNC(String()+" * "+c.String());
  DEBUG_VAR(this->m_max);
  DEBUG_VAR(c.m_max);
  double max=std::max(m_max, c.m_max);
  Color_Function_Decay c1(c);
  for (size_t i(0); i<c1.Internal().size(); ++i) {
    DEBUG_VAR(i);
    DEBUG_VAR(c1.m_internal[i]);
    if (c1.m_internal[i]<=max) continue;
    c1.ReplaceIndex(c1.m_internal[i], max+1+i);
  }
  c1.Conjugate();
  DEBUG_VAR(String());
  DEBUG_VAR(c1.String());

  Expression expr(String()+"*"+c1.String());
  if (msg_LevelIsDebugging()) expr.Print();
  if (expr.Evaluate()) {
    DEBUG_VAR(expr.Result());
    return expr.Result();
  }
  else {
    THROW(fatal_error, "Color factor evaluation failed: "+String()+
          " * "+c1.String());
  }
  return Complex(1.0,0.0);
}

namespace PHASIC {
  std::ostream &operator<<(std::ostream &os,const Color_Function_Decay &c)
  {
    os<<c.String();
    return os;
  }
}
