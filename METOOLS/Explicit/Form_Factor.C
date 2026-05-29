#include "METOOLS/Explicit/Form_Factor.H"
#include "METOOLS/Explicit/Current.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE METOOLS::Vertex_Key
#define OBJECT_TYPE METOOLS::Form_Factor
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

using namespace METOOLS;
using namespace ATOOLS;


Form_Factor::Form_Factor(const std::string &id,const Vertex_Key &key):
  m_id(id), p_v(key.p_v) {}

Form_Factor::~Form_Factor()
{}

std::ostream &METOOLS::operator<<(std::ostream &str,const Form_Factor &c)
{
  return str<<c.ID();
}

namespace METOOLS {
  class FFNone: public Form_Factor {
  public:
    FFNone(const Vertex_Key &key): Form_Factor("0",key) {}
    Complex FF() { return Complex(0.,0.); }
  };// end of class FFNone


  class FFPoint: public Form_Factor {
  public:
    FFPoint(const Vertex_Key &key): Form_Factor("Point",key) {}
    Complex FF() { return Complex(1.,0.); }
  };// end of class FFPoint
}// end of namespace METOOLS

DECLARE_GETTER(FFNone,"None",Form_Factor,Vertex_Key);
Form_Factor *Getter<Form_Factor,Vertex_Key,FFNone>::
operator()(const Vertex_Key &args) const {
  msg_Out()<<METHOD<<"\n";
  return new FFNone(args);
}

void ATOOLS::Getter<Form_Factor,Vertex_Key,FFNone>::
PrintInfo(std::ostream &str,const size_t width) const { str<<"0"; }

DECLARE_GETTER(FFPoint,"Point",Form_Factor,Vertex_Key);
Form_Factor *Getter<Form_Factor,Vertex_Key,FFPoint>::
operator()(const Vertex_Key &args) const {
  msg_Out()<<METHOD<<"\n";
  return new FFPoint(args);
}

void ATOOLS::Getter<Form_Factor,Vertex_Key,FFPoint>::
PrintInfo(std::ostream &str,const size_t width) const { str<<"Point"; }


namespace METOOLS {
  class FFDipole: public Form_Factor {
  private:
    int    m_pos;
    double m_norm, m_Lambda2;
  public:
    FFDipole(const Vertex_Key &key);
    Complex FF();
  };// end of class FFPoint
}// end of namespace METOOLS


#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include <set>

/////////////////////////////////////////////////////////////////////
//
// Simple dipole form factor
//
/////////////////////////////////////////////////////////////////////
FFDipole::FFDipole(const Vertex_Key &key):
  Form_Factor("Dipole",key), m_pos(-1), m_norm(0.), m_Lambda2(0.) {
  // Find the position of the (elementary) particle:
  // Normally this would be the photon, but I'll assume it is anything
  // that is not a hadron.
  // Also: extract the flavours in the vertex to look up parameters
  // in the model look-up table.
  std::set<kf_code> kfs;
  for (size_t i(0);i<key.m_j.size();++i) {
    if (!key.m_j[i]->Flav().IsHadron()) m_pos=i;
    kfs.insert(key.m_j[i]->Flav().Kfcode());
  }
  std::string tag = "";
  for (size_t i(0);i<key.p_mv->in.size();++i) 
    tag += std::string("_")+ToString(key.p_mv->in[i].Kfcode());
  msg_Out()<<METHOD<<": tag = {"<<tag<<"}\n";
  m_Lambda2 = key.p_model->ScalarConstant("Lambda2"+tag);
  m_norm    = key.p_model->ScalarConstant("Q"+tag);
}

Complex FFDipole::FF() {
  if (m_norm==0. || m_Lambda2==0.) return 0.;
  Current *j(m_pos<0?p_v->JC():p_v->J(m_pos));
  double Q2(-j->P().Abs2());
  return Complex(m_norm/sqr(1.+Q2/m_Lambda2),0.);
}

DECLARE_GETTER(FFDipole,"Dipole",Form_Factor,Vertex_Key);
Form_Factor *Getter<Form_Factor,Vertex_Key,FFDipole>::
operator()(const Vertex_Key &args) const {
  msg_Out()<<METHOD<<"\n";
  return new FFDipole(args);
}

void ATOOLS::Getter<Form_Factor,Vertex_Key,FFDipole>::
PrintInfo(std::ostream &str,const size_t width) const { str<<"Dipole"; }


