#include "DIRE/Shower/Alpha_QCD.H"

#include "MODEL/Main/Single_Vertex.H"
#include "DIRE/Shower/Shower.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

namespace DIRE {
  
  class qqqQ: public Alpha_QCD {
  private:

    double m_N;

    int m_mode;

  public:

    inline qqqQ(const Kernel_Key &key):
      Alpha_QCD(key), m_N(1.0),
      m_mode(0)
    {
    }

    double Scale(const Splitting &s) const
    {
      return s.m_t;
    }

    bool Allowed(const Splitting &s) const
    {
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
      s.p_c->SetColor(Color(s.m_cj[0].m_i,0));
      s.p_l->SetColor(Color(0,s.m_cj[0].m_j));
      s.p_n->SetColor(s.m_ci[0]);
      return 1;
    }

  };// end of class qqqQ

  class qQqq: public Alpha_QCD {
  private:

    double m_N;

    int m_mode;

  public:

    inline qQqq(const Kernel_Key &key):
      Alpha_QCD(key), m_N(1.0),
      m_mode(0)
    {
    }

    double Scale(const Splitting &s) const
    {
      return s.m_t;
    }

    bool Allowed(const Splitting &s) const
    {
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
      s.p_l->SetColor(Color(s.m_cj[0].m_i,0));
      s.p_c->SetColor(Color(0,s.m_cj[0].m_j));
      s.p_n->SetColor(s.m_ci[0]);
      return 1;
    }

  };// end of class qQqq

  class QqQQ: public Alpha_QCD {
  private:

    double m_N;

    int m_mode;

  public:

    inline QqQQ(const Kernel_Key &key):
      Alpha_QCD(key), m_N(1.0),
      m_mode(0)
    {
    }

    double Scale(const Splitting &s) const
    {
      return s.m_t;
    }

    bool Allowed(const Splitting &s) const
    {
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
      s.p_c->SetColor(Color(s.m_cj[0].m_i,0));
      s.p_l->SetColor(Color(0,s.m_cj[0].m_j));
      s.p_n->SetColor(s.m_ci[0]);
      return 1;
    }

  };// end of class QqQQ

  class QQQq: public Alpha_QCD {
  private:

    double m_N;

    int m_mode;

  public:

    inline QQQq(const Kernel_Key &key):
      Alpha_QCD(key), m_N(1.0),
      m_mode(0)
    {
    }

    double Scale(const Splitting &s) const
    {
      return s.m_t;
    }

    bool Allowed(const Splitting &s) const
    {
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
      s.p_l->SetColor(Color(s.m_cj[0].m_i,0));
      s.p_c->SetColor(Color(0,s.m_cj[0].m_j));
      s.p_n->SetColor(s.m_ci[0]);
      return 1;
    }

  };// end of class QqQQ

}// end of namespace DIRE

using namespace DIRE;

DECLARE_GETTER(qqqQ,"QCD{-3}{3}{3}{-3}",Gauge,Kernel_Key);

Gauge *ATOOLS::Getter<Gauge,Kernel_Key,qqqQ>::
operator()(const Parameter_Type &key) const
{
  return new qqqQ(key);
}

void ATOOLS::Getter<Gauge,Kernel_Key,qqqQ>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"q->qqQ";
}

DECLARE_GETTER(qQqq,"QCD{-3}{-3}{3}{3}",Gauge,Kernel_Key);

Gauge *ATOOLS::Getter<Gauge,Kernel_Key,qQqq>::
operator()(const Parameter_Type &key) const
{
  return new qQqq(key);
}

void ATOOLS::Getter<Gauge,Kernel_Key,qQqq>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"q->Qqq";
}

DECLARE_GETTER(QqQQ,"QCD{3}{3}{-3}{-3}",Gauge,Kernel_Key);

Gauge *ATOOLS::Getter<Gauge,Kernel_Key,QqQQ>::
operator()(const Parameter_Type &key) const
{
  return new QqQQ(key);
}

void ATOOLS::Getter<Gauge,Kernel_Key,QqQQ>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Q->qQQ";
}

DECLARE_GETTER(QQQq,"QCD{3}{-3}{-3}{3}",Gauge,Kernel_Key);

Gauge *ATOOLS::Getter<Gauge,Kernel_Key,QQQq>::
operator()(const Parameter_Type &key) const
{
  return new QQQq(key);
}

void ATOOLS::Getter<Gauge,Kernel_Key,QQQq>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Q->QQq";
}
