#include "DIRE/Shower/Alpha_QCD.H"

#include "MODEL/Main/Single_Vertex.H"
#include "DIRE/Shower/Shower.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

namespace DIRE {
  
  class GqQ: public Alpha_QCD {
  private:

    double m_N;

    int m_mode;

  public:

    inline GqQ(const Kernel_Key &key):
      Alpha_QCD(key), m_N(1.0),
      m_mode(key.p_v->in[1+key.m_mode].IsAnti())
    {
      if (key.m_type&1) m_N=(m_Nc*m_Nc-1.0)/m_Nc;
    }

    double Scale(const Splitting &s) const
    {
      return s.m_t;
    }

    bool Allowed(const Splitting &s) const
    {
      if (s.p_n) return s.p_s->Flav().Strong();
      Color cij(s.p_c->Col()), ck(s.p_s->Col());
      if (s.m_cm==0 && cij.m_i==ck.m_j) return true;
      if (s.m_cm==1 && cij.m_j==ck.m_i) return true;
      return false;
    }

    double Charge(const Splitting &s) const
    {
      return m_TR/2.0*m_N;
    }

    double Weight(const Splitting &s) const
    {
      return m_TR/2.0*m_N;
    }

    bool GeneratePoint(Splitting &s) const
    {
      s.ResetCol();
      Color ci(-1,0), cj(0,-1);
      ci.m_i=s.p_c->Col().m_i;
      cj.m_j=s.p_c->Col().m_j;
      if (m_mode) std::swap<Color>(ci,cj);
      s.AddCol(ci,cj);
      return true;
    }

    int Construct(Splitting &s) const
    {
      s.p_c->SetColor(s.m_ci[0]);
      s.p_n->SetColor(s.m_cj[0]);
      return 1;
    }

  };// end of class GqQ

  class GQq;

}// end of namespace DIRE

using namespace DIRE;

DECLARE_GETTER(GqQ,"QCD{8}{3}{-3}",Gauge,Kernel_Key);

Gauge *ATOOLS::Getter<Gauge,Kernel_Key,GqQ>::
operator()(const Parameter_Type &key) const
{
  return new GqQ(key);
}

void ATOOLS::Getter<Gauge,Kernel_Key,GqQ>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"G->qQ";
}

DECLARE_GETTER(GQq,"QCD{8}{-3}{3}",Gauge,Kernel_Key);

Gauge *ATOOLS::Getter<Gauge,Kernel_Key,GQq>::
operator()(const Parameter_Type &key) const
{
  return new GqQ(key);
}

void ATOOLS::Getter<Gauge,Kernel_Key,GQq>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"G->Qq";
}
