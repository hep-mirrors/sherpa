#include "ALARIC++/Shower/Lorentz_IS.H"

#include "MODEL/Main/Single_Vertex.H"
#include "ALARIC++/Shower/Shower.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace ALARIC {
  
  class FFV_IS: public Lorentz_IS {
  private:

    double m_jmax;

  public:

    inline FFV_IS(const Kernel_Key &key):
      Lorentz_IS(key), m_jmax(m_fl[0].Kfcode()<3?5.0:2.0) {}

    double Value(const Splitting &s) const
    {
      double B=1.0-s.m_zi;
      if (s.m_mec&1) B=1.-s.m_x-2.*s.m_y*(1.-s.m_y/(1.-s.m_x));
      B*=1.0+p_sk->GF()->K(s);
      return B;
    }

    double Integral(const Splitting &s) const
    {
      double I=0.5*sqr(1.0-s.m_eta);
      return I*m_jmax;
    }

    double Estimate(const Splitting &s) const
    {
      double E=1.0-s.m_z;
      return E*m_jmax;
    }

    bool GeneratePoint(Splitting &s) const
    {
      s.m_z=1.0-sqrt(ran->Get())*(1.0-s.m_eta);
      s.m_phi=2.0*M_PI*ran->Get();
      return true;
    }

  };// end of class FFV_IS

  class FVF_IS: public Lorentz_IS {
  private:

    double m_jmax;

  public:

    inline FVF_IS(const Kernel_Key &key):
      Lorentz_IS(key), m_jmax(5.0) {}

    double Value(const Splitting &s) const
    {
      double B=2.0*(1.0-s.m_zi)/s.m_zi+s.m_zi;
      return B;
    }

    double Integral(const Splitting &s) const
    {
      double I=2.0*log(1.0/s.m_eta);
      return I*m_jmax*PDFEstimate(s);
    }

    double Estimate(const Splitting &s) const
    {
      double E=2.0/s.m_z;
      return E*m_jmax*PDFEstimate(s);
    }

    bool GeneratePoint(Splitting &s) const
    {
      s.m_z=pow(s.m_eta,ran->Get());
      s.m_phi=2.0*M_PI*ran->Get();
      return true;
    }

  };// end of class FVF_IS

  class VFF_IS: public Lorentz_IS {
  private:

    double m_jmax;

  public:

    inline VFF_IS(const Kernel_Key &key):
      Lorentz_IS(key), m_jmax(5.0) {}

    double Value(const Splitting &s) const
    {
      double B=1.0-2.0*s.m_zi*(1.0-s.m_zi);
      if (s.m_mec&1) B=1.0-2.0*s.m_x*(1.0-s.m_x)
		       +s.m_y*(s.m_y+2.*s.m_x);
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

  };// end of class VFF_IS

}// end of namespace ALARIC

using namespace ALARIC;

DECLARE_GETTER(FFV_IS,"IS_FFV",Lorentz,Kernel_Key);

Lorentz *ATOOLS::Getter<Lorentz,Kernel_Key,FFV_IS>::
operator()(const Parameter_Type &args) const
{
  if (args.m_type!=1) return NULL;
  if (args.m_swap) return NULL;
  if ((args.m_mode==0 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[1].IntSpin()==1 &&
       args.p_v->in[2].IntSpin()==2) ||
      (args.m_mode==1 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[2].IntSpin()==1 &&
       args.p_v->in[1].IntSpin()==2)) {
    return new FFV_IS(args);
  }
  if ((args.m_mode==0 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[1].IntSpin()==2 &&
       args.p_v->in[2].IntSpin()==1) ||
      (args.m_mode==1 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[2].IntSpin()==2 &&
       args.p_v->in[1].IntSpin()==1)) {
    return new VFF_IS(args);
  }
  if (args.p_v->in[0].IntSpin()==2 &&
      args.p_v->in[1].IntSpin()==1 &&
      args.p_v->in[2].IntSpin()==1) {
    return new FVF_IS(args);
  }
  return NULL;
}

void ATOOLS::Getter<Lorentz,Kernel_Key,FFV_IS>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFV Lorentz Function";
}
