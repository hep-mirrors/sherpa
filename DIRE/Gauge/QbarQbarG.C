#include "DIRE/Shower/Alpha_QCD.H"

#include "MODEL/Main/Single_Vertex.H"
#include "DIRE/Shower/Shower.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

namespace DIRE {
  
  class QQG: public Alpha_QCD {
  private:

    double m_C[3], m_N;

    int m_mode;

  public:

    inline QQG(const Kernel_Key &key):
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
	  s.p_s->Flav().StrongCharge()>0;
      Color cij(s.p_c->Col()), ck(s.p_s->Col());
      if (s.m_cm==1 && cij.m_j==ck.m_i) return true;
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
      Color ci(0,-1), cj;
      cj.m_j=s.p_c->Col().m_j;
      if (m_mode) std::swap<Color>(ci,cj);
      s.AddCol(ci,cj);
      return true;
    }

    int Construct(Splitting &s) const
    {
      if (m_mode) s.m_ci[0].m_i=s.m_cj[0].m_j=Flow::Counter();
      else s.m_ci[0].m_j=s.m_cj[0].m_i=Flow::Counter();
      s.p_c->SetColor(s.m_ci[0]);
      s.p_n->SetColor(s.m_cj[0]);
      return 1;
    }

  };// end of class QQG

  class QGQ {};

}// end of namespace DIRE

using namespace DIRE;

DECLARE_GETTER(QQG,"QCD{3}{-3}{8}",Gauge,Kernel_Key);

Gauge *ATOOLS::Getter<Gauge,Kernel_Key,QQG>::
operator()(const Parameter_Type &key) const
{
  return new QQG(key);
}

void ATOOLS::Getter<Gauge,Kernel_Key,QQG>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Q->QG";
}

DECLARE_GETTER(QGQ,"QCD{3}{8}{-3}",Gauge,Kernel_Key);

Gauge *ATOOLS::Getter<Gauge,Kernel_Key,QGQ>::
operator()(const Parameter_Type &key) const
{
  return new QQG(key);
}

void ATOOLS::Getter<Gauge,Kernel_Key,QGQ>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Q->GQ";
}
