#include "DIRE/Shower/Lorentz_IF.H"

#include "MODEL/Main/Single_Vertex.H"
#include "DIRE/Shower/Shower.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace DIRE {
  
  class FFV_IF: public Lorentz_IF {
  private:

    double m_jmax;

  public:

    inline FFV_IF(const Kernel_Key &key):
      Lorentz_IF(key), m_jmax(m_fl[0].Kfcode()<3?5.0:2.0) {}

    double Value(const Splitting &s) const
    {
      double A1=2.0*(1.0-s.m_z)/(sqr(1.0-s.m_z)+s.m_t/s.m_Q2);
      double B1=-(1.0+s.m_z);
      return A1*(1.0+p_sk->GF()->K(s))+B1;
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

  };// end of class FFV_IF

  class FVF_IF: public Lorentz_IF {
  private:

    double m_jmax;

  public:

    inline FVF_IF(const Kernel_Key &key):
      Lorentz_IF(key), m_jmax(5.0) {}

    double Value(const Splitting &s) const
    {
      double V=2.0*s.m_z/(sqr(s.m_z)+s.m_t/s.m_Q2)-(2.0-s.m_z);
      if (s.m_mk2==0.0) return V;
      V-=2.0*s.m_mk2/(s.m_Q2+s.m_mk2)*s.m_y/(1.0-s.m_y);
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

  };// end of class FVF_IF

  class VFF_IF: public Lorentz_IF {
  private:

    double m_jmax;

  public:

    inline VFF_IF(const Kernel_Key &key):
      Lorentz_IF(key), m_jmax(5.0) {}

    double Value(const Splitting &s) const
    {
      double V=1.0-2.0*s.m_z*(1.0-s.m_z);
      return V;
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

  };// end of class VFF_IF

}// end of namespace DIRE

using namespace DIRE;

DECLARE_GETTER(FFV_IF,"IF_FFV",Lorentz,Kernel_Key);

Lorentz *ATOOLS::Getter<Lorentz,Kernel_Key,FFV_IF>::
operator()(const Parameter_Type &args) const
{
  if (args.m_type!=1) return NULL;
  if ((args.m_mode==0 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[1].IntSpin()==1 &&
       args.p_v->in[2].IntSpin()==2) ||
      (args.m_mode==1 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[2].IntSpin()==1 &&
       args.p_v->in[1].IntSpin()==2)) {
    return new FFV_IF(args);
  }
  if ((args.m_mode==0 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[1].IntSpin()==2 &&
       args.p_v->in[2].IntSpin()==1) ||
      (args.m_mode==1 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[2].IntSpin()==2 &&
       args.p_v->in[1].IntSpin()==1)) {
    return new VFF_IF(args);
  }
  if (args.p_v->in[0].IntSpin()==2 &&
      args.p_v->in[1].IntSpin()==1 &&
      args.p_v->in[2].IntSpin()==1) {
    return new FVF_IF(args);
  }
  return NULL;
}

void ATOOLS::Getter<Lorentz,Kernel_Key,FFV_IF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFV Lorentz Function";
}
