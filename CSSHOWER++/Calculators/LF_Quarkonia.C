#include "CSSHOWER++/Showers/Splitting_Function_Base.H"

namespace CSSHOWER {

    class LF_FFV_Quarkonia_FF: public SF_Lorentz {
    public:
    
      inline LF_FFV_Quarkonia_FF(const SF_Key &key): SF_Lorentz(key) {}

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

double LF_FFV_Quarkonia_FF::operator()
  (const double zz, const double y, const double eta,
   const double scale, const double Q2)
{
  const double z  = 1 - zz;
  double muij2 = sqr(p_ms->Mass(m_flavs[0]))/Q2;
  double mi2   = sqr(p_ms->Mass(m_flavs[1]));
  double mj2   = sqr(p_ms->Mass(m_flavs[2]));
  double mk2   = sqr(p_ms->Mass(m_flspec));
  double mui2  = mi2/Q2, muj2 = mj2/Q2, muk2 = mk2/Q2;
  //the massless case
  if (muij2==0. && mui2==0. && muk2==0.) {
    msg_Error() << "Cannot make massless quarkonia emission" << std::endl;
    exit(EXIT_FAILURE);
  }
  else {
    //the massive case
    double vtijk = Lambda(1.,muij2,muk2), vijk = sqr(2.*muk2+(1.-mui2-muk2)*(1.-y))-4.*muk2;
    double vkji = 1.-4.*mui2*((1.-mui2-muk2)*(1.-y)*(1.-z)+muk2)/sqr((1.-mui2-muk2)*(y+(1.-y)*z));
    if (vtijk<0.0 || vijk<0.0 || vkji<0.0) return 0.0;
    vtijk = sqrt(vtijk)/(1.-muij2-muk2);
    vijk  = sqrt(vijk)/((1.-mui2-muk2)*(1.-y));
    vkji  = sqrt(vkji);
    double pipj  = Q2*(1.0-mui2-muk2)*y/2.0;
    double pkpj  = Q2*(1.0-mui2-muk2)*(1.-y)*(1.-z)/2.0;
    double massive = 1./sqr(Q2)/sqr(sqr(1-muij2))*((1-2*muij2-47*sqr(muij2)) - z*(1-muij2)*(1-sqr(sqrt(mui2)+sqrt(muj2))) + 4*(z*(1-z))/(2-z)*(1-mui2) 
                     - 4*(8-7*z-5*z*z)/(2-z)*muij2*(1-muij2) + 12*(z*z*(1-z))/sqr(2-z)*sqr(1-muij2));
    massive *= 1./((1.-mui2-muk2)+1./y*(mui2-muij2));
    msg_Out() << METHOD << ", massive: " << massive << std::endl;
    double longpol = 0.5 * ( 1. - z );
    double value = 8/27/M_PI/(Q2*mui2)*p_cf->Coupling(scale,1) * massive + p_cf->Coupling(scale,0) * longpol;
    msg_Out() << METHOD << ", value: " << value << std::endl;
    msg_Out() << METHOD << ", JFF: " << JFF(y,mui2,muj2,muk2,muij2) << std::endl;
    return value * JFF(y,mui2,muj2,muk2,muij2);
  } 
}

double LF_FFV_Quarkonia_FF::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  return (4.0*p_cf->MaxCoupling(0) + 0.5*p_cf->MaxCoupling(1)) *log((1.-zmin)/(1.-zmax));
}

double LF_FFV_Quarkonia_FF::OverEstimated(const double z,const double y)
{
  return (4.0*p_cf->MaxCoupling(0) + 0.5*p_cf->MaxCoupling(1))/(1.-z);
}

double LF_FFV_Quarkonia_FF::Z()
{
  return 1.-(1.-m_zmin)*pow((1.-m_zmax)/(1.-m_zmin),ATOOLS::ran->Get());
}

DECLARE_GETTER(LF_FFV_Quarkonia_FF,"FFV_Quarkonia",SF_Lorentz,SF_Key);

SF_Lorentz *ATOOLS::Getter<SF_Lorentz,SF_Key,LF_FFV_Quarkonia_FF>::
operator()(const Parameter_Type &args) const
{
  if (args.m_col<0) return NULL;
  if ((args.m_mode==0 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[1].IntSpin()==1 &&
       args.p_v->in[2].IntSpin()==2) ||
      (args.m_mode==1 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[2].IntSpin()==1 &&
       args.p_v->in[1].IntSpin()==2)) {
    switch (args.m_type) {
    case cstp::FF: return new LF_FFV_Quarkonia_FF(args);
    // case cstp::FI: return new LF_FFV_FI(args);
    // case cstp::IF: return new LF_FFV_IF(args);
    // case cstp::II: return new LF_FFV_II(args);
    case cstp::none: break;
    }
  }
  // if ((args.m_mode==0 &&
  //      args.p_v->in[0].IntSpin()==1 &&
  //      args.p_v->in[1].IntSpin()==2 &&
  //      args.p_v->in[2].IntSpin()==1) ||
  //     (args.m_mode==1 &&
  //      args.p_v->in[0].IntSpin()==1 &&
  //      args.p_v->in[2].IntSpin()==2 &&
  //      args.p_v->in[1].IntSpin()==1)) {
  //   switch (args.m_type) {
  //   case cstp::FF: return new LF_FVF_Quarkonia_FF(args);
  //   // case cstp::FI: return new LF_FVF_FI(args);
  //   // case cstp::IF: return new LF_FVF_IF(args);
  //   // case cstp::II: return new LF_FVF_II(args);
  //   case cstp::none: break;
  //   }
  // }
  // if (args.p_v->in[0].IntSpin()==2 &&
  //     args.p_v->in[1].IntSpin()==1 &&
  //     args.p_v->in[2].IntSpin()==1) {
  //   switch (args.m_type) {
  //   case cstp::FF: return new LF_VFF_Quarkonia_FF(args);
  //   // case cstp::FI: return new LF_VFF_FI(args);
  //   // case cstp::IF: return new LF_VFF_IF(args);
  //   // case cstp::II: return new LF_VFF_II(args);
  //   case cstp::none: break;
  //   }
  // }
  return NULL;
}

void ATOOLS::Getter<SF_Lorentz,SF_Key,LF_FFV_Quarkonia_FF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"ffv lorentz functions";
}