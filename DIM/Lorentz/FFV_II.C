#include "DIM/Shower/Lorentz_II.H"

#include "MODEL/Main/Single_Vertex.H"
#include "DIM/Shower/Shower.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace DIM {
  
  class FFV_II: public Lorentz_II {
  private:

    double m_jmax;

  public:

    inline FFV_II(const Kernel_Key &key):
      Lorentz_II(key), m_jmax(m_fl[0].Kfcode()<3?5.0:2.0) {}

    double Value(const Splitting &s) const
    {
      double A=2.0*(1.0-s.m_z)/(sqr(1.0-s.m_z)+s.m_t/s.m_Q2);
      double B=-(1.0+s.m_z);
      return A*(1.0+p_sk->GF()->K(s))+B;
    }

    double Integral(const Splitting &s) const
    {
      double I=log(1.0+sqr(1.0-s.m_eta)*s.m_Q2/s.m_t0);
      return I*(1.0+p_sk->GF()->KMax(s))*m_jmax;
    }

    double Estimate(const Splitting &s) const
    {
      double E=2.0*(1.0-s.m_z)/(sqr(1.0-s.m_z)+s.m_t0/s.m_Q2);
      return E*(1.0+p_sk->GF()->KMax(s))*m_jmax;
    }

    bool GeneratePoint(Splitting &s) const
    {
      double k2(s.m_t0/s.m_Q2);
      s.m_z=1.0-sqrt(k2*(pow(1.0+sqr(1.0-s.m_eta)/k2,ran->Get())-1.0));
      s.m_phi=2.0*M_PI*ran->Get();
      return true;
    }

  };// end of class FFV_II

  class FVF_II: public Lorentz_II {
  private:

    double m_jmax;

  public:

    inline FVF_II(const Kernel_Key &key):
      Lorentz_II(key), m_jmax(5.0) {}

    double Value(const Splitting &s) const
    {
      double V=2.0*s.m_z/(sqr(s.m_z)+s.m_t/s.m_Q2)-(2.0-s.m_z);
      return V;
    }

    double Integral(const Splitting &s) const
    {
      double I=log((s.m_Q2+s.m_t0)/(s.m_Q2*sqr(s.m_eta)+s.m_t0));
      return I*m_jmax*PDFEstimate(s);
    }

    double Estimate(const Splitting &s) const
    {
      double E=2.0*s.m_z/(sqr(s.m_z)+s.m_t0/s.m_Q2);
      return E*m_jmax*PDFEstimate(s);
    }

    bool GeneratePoint(Splitting &s) const
    {
      double k2(s.m_t0/s.m_Q2);
      s.m_z=sqrt(pow((1.0+k2)/(sqr(s.m_eta)+k2),-ran->Get())*(1.0+k2)-k2);
      s.m_phi=2.0*M_PI*ran->Get();
      return true;
    }

  };// end of class FVF_II

  class VFF_II: public Lorentz_II {
  private:

    double m_jmax;

  public:

    inline VFF_II(const Kernel_Key &key):
      Lorentz_II(key), m_jmax(5.0) {}

    double Value(const Splitting &s) const
    {
      double B=1.0-2.0*s.m_z*(1.0-s.m_z);
      return B;
    }

    double Integral(const Splitting &s) const
    {
      return (1.0-s.m_eta)*m_jmax*PDFEstimate(s);
    }

    double Estimate(const Splitting &s) const
    {
      return m_jmax*PDFEstimate(s);
    }

    bool GeneratePoint(Splitting &s) const
    {
      s.m_z=s.m_eta+(1.0-s.m_eta)*ran->Get();
      s.m_phi=2.0*M_PI*ran->Get();
      return true;
    }

  };// end of class VFF_II

}// end of namespace DIM

using namespace DIM;

DECLARE_GETTER(FFV_II,"II_FFV",Lorentz,Kernel_Key);

Lorentz *ATOOLS::Getter<Lorentz,Kernel_Key,FFV_II>::
operator()(const Parameter_Type &args) const
{
  if (args.m_type!=3) return NULL;
  if ((args.m_mode==0 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[1].IntSpin()==1 &&
       args.p_v->in[2].IntSpin()==2) ||
      (args.m_mode==1 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[2].IntSpin()==1 &&
       args.p_v->in[1].IntSpin()==2)) {
    return new FFV_II(args);
  }
  if ((args.m_mode==0 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[1].IntSpin()==2 &&
       args.p_v->in[2].IntSpin()==1) ||
      (args.m_mode==1 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[2].IntSpin()==2 &&
       args.p_v->in[1].IntSpin()==1)) {
    return new VFF_II(args);
  }
  if (args.p_v->in[0].IntSpin()==2 &&
      args.p_v->in[1].IntSpin()==1 &&
      args.p_v->in[2].IntSpin()==1) {
    return new FVF_II(args);
  }
  return NULL;
}

void ATOOLS::Getter<Lorentz,Kernel_Key,FFV_II>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFV Lorentz Function";
}
