#include "METOOLS/Explicit/Form_Factor.H"
#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Vertex.H"

#include "MODEL/Main/Single_Vertex.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include <set>

namespace METOOLS {
  class FFDirac_F1: public Form_Factor {
  private:
    int    m_pos;
    double m_q, m_mu, m_Lambda2, m_mass2;
  public:
    FFDirac_F1(const Vertex_Key &key);
    Complex FF();
  };// end of class FFPoint

  class FFDirac_F2: public Form_Factor {
  private:
    int    m_pos;
    double m_q, m_mu, m_Lambda2, m_mass2;
  public:
    FFDirac_F2(const Vertex_Key &key);
    Complex FF();
  };// end of class FFPoint
}// end of namespace METOOLS


using namespace METOOLS;
using namespace ATOOLS;
using namespace std;
/////////////////////////////////////////////////////////////////////
//
// Dirac F_1 and F_2 form factors, assuming a dipole form
//
// Normalisations:
// F_1(0) = q      for proton, neutron,
//                 where q is the charge (=1, 0)
// F_2(0) = mu-q   for proton, neutron,
//                 where mu is anomalous magnetic moment
//
/////////////////////////////////////////////////////////////////////
FFDirac_F1::FFDirac_F1(const Vertex_Key &key):
  Form_Factor("Dirac_F1",key),
  m_pos(-1), m_q(0.), m_mu(0.), m_Lambda2(0.), m_mass2(0.) {
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
  m_Lambda2 = key.p_model->ScalarConstant("Lambda2"+tag);
  m_q       = key.p_model->ScalarConstant("Q"+tag);
  m_mu      = key.p_model->ScalarConstant("Mu"+tag);
  m_mass2   = sqr(key.m_j[0]->Flav().HadMass());
  msg_Out()<<METHOD<<"("<<key.m_j.size()<<"): tag = {"<<tag<<"}: "
	   <<"q = "<<m_q<<", mu = "<<m_mu<<", "
	   <<"Lambda2 = "<<m_Lambda2<<", mass = "<<sqrt(m_mass2)<<".\n";
}

Complex FFDirac_F1::FF() {
  if (m_q==0. && m_mu==0.) return 0.;
  Current *j = m_pos<0?p_v->JC():p_v->J(m_pos);
  // there was a minus sign before in Q2.
  double Q2  = j->P().Abs2();
  double dip = 1./sqr(1.+Q2/m_Lambda2);
  double tau = Q2/(4.*m_mass2);
  return Complex((m_q+tau*m_mu)/(1.+tau)*dip,0.);
}

FFDirac_F2::FFDirac_F2(const Vertex_Key &key):
  Form_Factor("Dirac_F2",key),
  m_pos(-1), m_q(0.), m_mu(0.), m_Lambda2(0.), m_mass2(0.) {
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
  m_Lambda2 = key.p_model->ScalarConstant("Lambda2"+tag);
  m_q       = key.p_model->ScalarConstant("Q"+tag);
  m_mu      = key.p_model->ScalarConstant("Mu"+tag);
  m_mass2   = sqr(key.m_j[0]->Flav().HadMass());
  msg_Out()<<METHOD<<"("<<key.m_j.size()<<"): tag = {"<<tag<<"}: "
	   <<"q = "<<m_q<<", mu = "<<m_mu<<", "
	   <<"Lambda2 = "<<m_Lambda2<<", mass = "<<sqrt(m_mass2)<<".\n";
}

Complex FFDirac_F2::FF() {
  if (m_q==0. && m_mu==0.) return 0.;
  Current *j = m_pos<0?p_v->JC():p_v->J(m_pos);
  // there was a minus sign before in Q2.
  double Q2  = j->P().Abs2();
  double dip = 1./sqr(1.+Q2/m_Lambda2);
  double tau = Q2/(4.*m_mass2);
  return Complex((m_mu-m_q)/(1.+tau)*dip,0.);
}

DECLARE_GETTER(FFDirac_F1,"Dirac_F1",Form_Factor,Vertex_Key);
Form_Factor *Getter<Form_Factor,Vertex_Key,FFDirac_F1>::
operator()(const Vertex_Key &args) const {
  msg_Out()<<METHOD<<"\n";
  return new FFDirac_F1(args);
}

void ATOOLS::Getter<Form_Factor,Vertex_Key,FFDirac_F1>::
PrintInfo(std::ostream &str,const size_t width) const { str<<"Dirac_F1"; }

DECLARE_GETTER(FFDirac_F2,"Dirac_F2",Form_Factor,Vertex_Key);
Form_Factor *Getter<Form_Factor,Vertex_Key,FFDirac_F2>::
operator()(const Vertex_Key &args) const {
  msg_Out()<<METHOD<<"\n";
  return new FFDirac_F2(args);
}

void ATOOLS::Getter<Form_Factor,Vertex_Key,FFDirac_F2>::
PrintInfo(std::ostream &str,const size_t width) const { str<<"Dirac_F2"; }
