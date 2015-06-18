#include "DIRE/Shower/Alpha_QCD.H"

#include "DIRE/Shower/Shower.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

namespace DIRE {
  
  class GGG: public Alpha_QCD {
  private:

    int m_mode;

  public:

    inline GGG(const Kernel_Key &key):
      Alpha_QCD(key), m_mode(key.m_mode) {}

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
      return m_CA/2.0;
    }

    double Weight(const Splitting &s) const
    {
      return m_CA/2.0;
    }

    bool GeneratePoint(Splitting &s) const
    {
      s.ResetCol();
      int mcs(s.p_c->Col().m_j==s.p_s->Col().m_i);
      int msc(s.p_c->Col().m_i==s.p_s->Col().m_j);
      if (mcs&&msc) {
	if (ran->Get()>0.5) mcs=0; else msc=0;
      }
      if (mcs) {
	Color ci(s.p_c->Col().m_i,-1);
	Color cj(-1,s.p_c->Col().m_j);
	s.AddCol(ci,cj);
      }
      if (msc) {
	Color ci(-1,s.p_c->Col().m_j);
	Color cj(s.p_c->Col().m_i,-1);
	s.AddCol(ci,cj);
      }
      return true;
    }

    int Construct(Splitting &s) const
    {
      int nc=Flow::Counter();
      if (s.m_ci[0].m_i<0) s.m_ci[0].m_i=nc;
      if (s.m_ci[0].m_j<0) s.m_ci[0].m_j=nc;
      if (s.m_cj[0].m_i<0) s.m_cj[0].m_i=nc;
      if (s.m_cj[0].m_j<0) s.m_cj[0].m_j=nc;
      if (m_mode && (m_type&1))
	std::swap<Color>(s.m_ci[0],s.m_cj[0]);
      for (size_t i(0);i<s.m_ci.size();++i) {
	s.p_c->SetColor(s.m_ci[i]);
	s.p_n->SetColor(s.m_cj[i]);
      }
      return 1;
    }

  };// end of class GGG

}// end of namespace DIRE

using namespace DIRE;

DECLARE_GETTER(GGG,"QCD{8}{8}{8}",Gauge,Kernel_Key);

Gauge *ATOOLS::Getter<Gauge,Kernel_Key,GGG>::
operator()(const Parameter_Type &key) const
{
  return new GGG(key);
}

void ATOOLS::Getter<Gauge,Kernel_Key,GGG>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"G->GG";
}
