#include "CSSHOWER++/Showers/Splitting_Function_Base.H"

#define USING_DIS_MEC 

namespace CSSHOWER {
  
  class LF_VVH_FF: public SF_Lorentz {
  public:

    inline LF_VVH_FF(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

  class LF_VVH_FI: public SF_Lorentz {
  protected:

    double m_Jmax;

  public:

    inline LF_VVH_FI(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

  class LF_VHV_FF: public SF_Lorentz {
  public:

    inline LF_VHV_FF(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

  class LF_VHV_FI: public SF_Lorentz {
  protected:

    double m_Jmax;

  public:

    inline LF_VHV_FI(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

}

#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Math/Random.H"

using namespace CSSHOWER;
using namespace ATOOLS;

double LF_VVH_FF::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2)
{
  double mj2 = sqr(p_ms->Mass(m_flavs[2]));
  return 2.0 * p_cf->Coupling(scale,0) * mj2 * JFF(y,0.0,0.0,0.0,0.0);
}

double LF_VVH_FF::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  double mj2 = sqr(p_ms->Mass(m_flavs[2]));
  return 2.0*p_cf->MaxCoupling(0)*mj2;
}

double LF_VVH_FF::OverEstimated(const double z,const double y)
{
  double mj2 = sqr(p_ms->Mass(m_flavs[2]));
  return 2.0*p_cf->MaxCoupling(0)*mj2;
}

double LF_VVH_FF::Z()
{
  return m_zmin+(m_zmax-m_zmin)*ATOOLS::ran->Get();
}

double LF_VVH_FI::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2)
{  
  double mj2 = sqr(p_ms->Mass(m_flavs[2]));
  return 2.0 * p_cf->Coupling(scale,0) * mj2 * JFF(y,0.0,0.0,0.0,0.0);
}

double LF_VVH_FI::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  double mj2 = sqr(p_ms->Mass(m_flavs[2]));
  return 2.0*p_cf->MaxCoupling(0)*mj2;
}

double LF_VVH_FI::OverEstimated(const double z,const double y)
{
  double mj2 = sqr(p_ms->Mass(m_flavs[2]));
  return 2.0*p_cf->MaxCoupling(0)*mj2;
}

double LF_VVH_FI::Z()
{
  return m_zmin+(m_zmax-m_zmin)*ATOOLS::ran->Get();
}

double LF_VHV_FF::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2)
{
  double mi2 = sqr(p_ms->Mass(m_flavs[1]));
  return 2.0 * p_cf->Coupling(scale,0) * mi2 * JFF(y,0.0,0.0,0.0,0.0);
}

double LF_VHV_FF::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  double mi2 = sqr(p_ms->Mass(m_flavs[1]));
  return 2.0*p_cf->MaxCoupling(0)*mi2;
}

double LF_VHV_FF::OverEstimated(const double z,const double y)
{
  double mi2 = sqr(p_ms->Mass(m_flavs[1]));
  return 2.0*p_cf->MaxCoupling(0)*mi2;
}

double LF_VHV_FF::Z()
{
  return m_zmin+(m_zmax-m_zmin)*ATOOLS::ran->Get();
}

double LF_VHV_FI::operator()
  (const double z,const double y,
   const double eta, const double scale,const double Q2)
{
  double mi2 = sqr(p_ms->Mass(m_flavs[1]));
  return 2.0 * p_cf->Coupling(scale,0) * mi2 * JFF(y,0.0,0.0,0.0,0.0);
}
double LF_VHV_FI::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  double mi2 = sqr(p_ms->Mass(m_flavs[1]));
  return 2.0*p_cf->MaxCoupling(0)*mi2;
}

double LF_VHV_FI::OverEstimated(const double z,const double y)
{
  double mi2 = sqr(p_ms->Mass(m_flavs[1]));
  return 2.0*p_cf->MaxCoupling(0)*mi2;
}

double LF_VHV_FI::Z()
{
  return m_zmin+(m_zmax-m_zmin)*ATOOLS::ran->Get();
}

DECLARE_GETTER(LF_VVH_FF,"Triangle",SF_Lorentz,SF_Key);

SF_Lorentz *ATOOLS::Getter<SF_Lorentz,SF_Key,LF_VVH_FF>::
operator()(const Parameter_Type &args) const
{
  if (args.m_col<0) return NULL;
  if ((args.m_mode==0 &&
       args.p_v->in[0].IntSpin()==2 &&
       args.p_v->in[1].IntSpin()==2 &&
       args.p_v->in[2].IntSpin()==0) ||
      (args.m_mode==1 &&
       args.p_v->in[0].IntSpin()==2 &&
       args.p_v->in[2].IntSpin()==2 &&
       args.p_v->in[1].IntSpin()==0)) {
    switch (args.m_type) {
    case cstp::FF: return new LF_VVH_FF(args);
    case cstp::FI: return new LF_VVH_FI(args);
    case cstp::IF: return NULL;
    case cstp::II: return NULL;
    case cstp::none: break;
    }
  }
  if ((args.m_mode==0 &&
       args.p_v->in[0].IntSpin()==2 &&
       args.p_v->in[1].IntSpin()==0 &&
       args.p_v->in[2].IntSpin()==2) ||
      (args.m_mode==1 &&
       args.p_v->in[0].IntSpin()==2 &&
       args.p_v->in[2].IntSpin()==0 &&
       args.p_v->in[1].IntSpin()==2)) {
    switch (args.m_type) {
    case cstp::FF: return new LF_VHV_FF(args);
    case cstp::FI: return new LF_VHV_FI(args);
    case cstp::IF: return NULL;
    case cstp::II: return NULL;
    case cstp::none: break;
    }
  }
  return NULL;
}

void ATOOLS::Getter<SF_Lorentz,SF_Key,LF_VVH_FF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"triangle lorentz functions";
}
