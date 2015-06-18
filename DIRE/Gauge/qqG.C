#include "DIRE/Shower/Alpha_QCD.H"

#include "MODEL/Main/Single_Vertex.H"
#include "DIRE/Shower/Shower.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

namespace DIRE {
  
  class qqG: public Alpha_QCD {
  private:

    double m_C[3], m_N;

    int m_mode;

  public:

    inline qqG(const Kernel_Key &key):
      Alpha_QCD(key), m_N(1.0),
      m_mode(key.p_v->in[1+key.m_mode].IsGluon())
    {
      m_C[0]=(m_Nc-1.0)/(2.0*m_CF);
      m_C[1]=m_C[0]/(m_Nc*m_Nc);
      m_C[2]=m_C[1]*(m_Nc-1.0);
      if ((key.m_type&1) && m_mode) m_N=m_Nc/(m_Nc*m_Nc-1.0);
    }

    double Scale(const Splitting &s) const
    {
      if (s.m_cpl) return s.m_t;
      if ((m_type&1) && m_mode) return s.m_Q2*s.m_y/s.m_x;
      return s.m_t;
    }

    bool Allowed(const Splitting &s) const
    {
      if (s.p_n)
	return s.p_s->Flav().StrongCharge()==8 ||
	  s.p_s->Flav().StrongCharge()<0;
      Color cij(s.p_c->Col()), ck(s.p_s->Col());
      if (s.m_cm==0 && cij.m_i==ck.m_j) return true;
      return false;
    }

    double Charge(const Splitting &s) const
    {
      return m_CF*m_N;
    }

    double Weight(const Splitting &s) const
    {
      return m_CF*m_N;
    }

    bool GeneratePoint(Splitting &s) const
    {
      s.ResetCol();
      Color ci(-1,0), cj;
      cj.m_i=s.p_c->Col().m_i;
      if (m_mode) std::swap<Color>(ci,cj);
      s.AddCol(ci,cj);
      return true;
    }

    int Construct(Splitting &s) const
    {
      if (m_mode) s.m_ci[0].m_j=s.m_cj[0].m_i=Flow::Counter();
      else s.m_ci[0].m_i=s.m_cj[0].m_j=Flow::Counter();
      s.p_c->SetColor(s.m_ci[0]);
      s.p_n->SetColor(s.m_cj[0]);
      return 1;
    }

  };// end of class qqG

  class qGq {};

}// end of namespace DIRE

using namespace DIRE;

DECLARE_GETTER(qqG,"QCD{-3}{3}{8}",Gauge,Kernel_Key);

Gauge *ATOOLS::Getter<Gauge,Kernel_Key,qqG>::
operator()(const Parameter_Type &key) const
{
  return new qqG(key);
}

void ATOOLS::Getter<Gauge,Kernel_Key,qqG>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"q->qG";
}

DECLARE_GETTER(qGq,"QCD{-3}{8}{3}",Gauge,Kernel_Key);

Gauge *ATOOLS::Getter<Gauge,Kernel_Key,qGq>::
operator()(const Parameter_Type &key) const
{
  return new qqG(key);
}

void ATOOLS::Getter<Gauge,Kernel_Key,qGq>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"q->Gq";
}
