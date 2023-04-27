#include "METOOLS/Explicit/Form_Factor.H"

#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;

namespace METOOLS {

  class F1p: public Form_Factor {
  private:
    double m_lambda, m_mp, m_muP;
    int m_mode;
  public:
    F1p(const Vertex_Key &key): Form_Factor("F1p",key), m_mode(-1)
    {
      for (size_t i(0);i<key.m_j.size();++i)
	if (key.m_j[i]->Flav().IsPhoton()) m_mode=i;
      m_mp=Flavour(kf_p_plus).Mass();
      m_lambda=key.p_model->ScalarConstant("lambda");
      m_muP=key.p_model->ScalarConstant("Mu Proton");
    }
    double FF()
    {
      Current *j(m_mode<0?p_v->JC():p_v->J(m_mode));
      double Q2(-j->P().Abs2());
#ifdef DEBUG__BG
      msg_Debugging()<<METHOD<<"("<<j->Id()<<","<<j->Flav()<<"): {\n"
		     <<"  p = "<<j->P()<<" -> Q^2 = "<<Q2<<"\n"
		     <<"  \\lambda = "<<m_lambda<<", \\mu_P = "<<m_muP<<"\n";
#endif
      double tau=Q2/4/sqr(m_mp);
      double Gep=1.0/sqr(1.0+Q2/sqr(m_lambda));
      double Gmp=m_muP*Gep;
#ifdef DEBUG__BG
      msg_Debugging()<<"  Gep = "<<Gep<<", Gmp = "<<Gmp<<"\n"
		     <<"  F_1 = "<<(Gep+tau*Gmp)/(1+tau)<<"\n}\n";
#endif
      return (Gep+tau*Gmp)/(1+tau);
    }
  };// end of class F1p
    
  class F2p: public Form_Factor {
  private:
    double m_lambda, m_mp, m_muP;
    int m_mode;
  public:
    F2p(const Vertex_Key &key): Form_Factor("F2p",key), m_mode(-1)
    {
      for (size_t i(0);i<key.m_j.size();++i)
	if (key.m_j[i]->Flav().IsPhoton()) m_mode=i;
      m_mp=Flavour(kf_p_plus).Mass();
      m_lambda=key.p_model->ScalarConstant("lambda");
      m_muP=key.p_model->ScalarConstant("Mu Proton");
    }
    double FF()
    {
      Current *j(m_mode<0?p_v->JC():p_v->J(m_mode));
      double Q2(-j->P().Abs2());
#ifdef DEBUG__BG
      msg_Debugging()<<METHOD<<"("<<j->Id()<<","<<j->Flav()<<"): {\n"
		     <<"  p = "<<j->P()<<" -> Q^2 = "<<Q2<<"\n"
		     <<"  \\lambda = "<<m_lambda<<", \\mu_P = "<<m_muP<<"\n";
#endif
      double tau=Q2/4/sqr(m_mp);
      double Gep=1.0/sqr(1.0+Q2/sqr(m_lambda));
      double Gmp=m_muP*Gep;
#ifdef DEBUG__BG
      msg_Debugging()<<"  Gep = "<<Gep<<", Gmp = "<<Gmp<<"\n"
		     <<"  F_2 = "<<(Gmp-Gep)/(1+tau)<<"\n}\n";
#endif
      return (Gmp-Gep)/(1+tau);
    }
  };// end of class F2p

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(F1p,"F1p",Form_Factor,Vertex_Key);
Form_Factor *Getter<Form_Factor,Vertex_Key,F1p>::
operator()(const Vertex_Key &args) const
{ return new F1p(args); }
void ATOOLS::Getter<Form_Factor,Vertex_Key,F1p>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"F1p"; }

DECLARE_GETTER(F2p,"F2p",Form_Factor,Vertex_Key);
Form_Factor *Getter<Form_Factor,Vertex_Key,F2p>::
operator()(const Vertex_Key &args) const
{ return new F2p(args); }
void ATOOLS::Getter<Form_Factor,Vertex_Key,F2p>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"F2p"; }
