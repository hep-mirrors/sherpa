#include "ALARIC++/Shower/Lorentz_FS.H"

#include "MODEL/Main/Single_Vertex.H"
#include "ALARIC++/Shower/Kernel.H"
#include "ALARIC++/Tools/Amplitude.H"
#include "PHASIC++/Channels/Transverse_Kinematics.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace ALARIC {
  
  class Soft_FS: public Lorentz_FS_Rad {
  private:

    int m_id;

    struct ID_Params {
      double m_Q2, m_kap, m_I, m_zp;
      ATOOLS::Vec4D m_qa, m_pk, m_Kt;
      ID_Params(const Splitting &s):
	m_qa(s.p_c->Mom()), m_pk(s.p_s->Mom()) {
	m_Kt=s.m_Kt;
	m_Q2=2.*m_qa*m_Kt;
	m_kap=m_Kt.Abs2()/m_Q2;
	m_I=1.0;
	m_zp=1.0/(1.0+s.m_t0/m_Kt.Abs2());
      }
    };//end of struct ID_Params
    
  public:

    inline Soft_FS(const Kernel_Key &key):
      Lorentz_FS_Rad(key), m_id(key.p_v->in[1]==key.p_v->in[2]) {}

    double Value(const Splitting &s) const
    {
      Vec4D pi(s.m_pi), pk(s.m_pk), pj(s.m_pj);
      if (pk[0]<0.0) pk=-pk;
      Vec4D n(-s.m_Kt+s.p_c->Mom()-s.m_pi);
      double sij(pi*pj), sik(pi*pk), skj(pj*pk);
      double D(sij*(pk*n)+skj*(pi*n));
      if (D==0.0) return 0.0;
      double A(2*sik/(sij*skj)
	       -pi.Abs2()/sqr(sij)
	       -pk.Abs2()/sqr(skj));
      A*=sij*skj*(pi*n)/D;
      double sf(1.0);
      if (m_id) sf=p_sk->Mode()?1.0-s.m_z:s.m_z;
#ifdef DEBUG__Kinematics
      if (s.m_p.size()) {
	Vec4D K(s.m_p[0]+s.m_p[1]), pij(s.p_c->Mom());
	Vec4D n(K+pj), nb(n-n.Abs2()/(2.*pij*K)*pij);
	double y=(s.m_pi*s.m_pj)/(s.m_pi*n);
	double x=(s.m_pi*n)/((s.m_pi+s.m_pj)*n);
	double cp(CosPhi(pk,pj,n,pi));
	DEBUG_VAR(-y<<" "<<s.m_y<<" "<<-y/s.m_y-1.);
	DEBUG_VAR(x<<" "<<s.m_x<<" "<<x/s.m_x-1.);
	DEBUG_VAR(cp<<" "<<cos(s.m_phi)<<" "<<cp/cos(s.m_phi)-1.);
      }
#endif
      return sf*A*(1.0+p_sk->GF()->K(s)+p_sk->GF()->RenCT(s));
    }

    double Integral(const Splitting &s) const
    {
      ID_Params ip(s);
      double I(ip.m_I*4.0*log(1./(1.-ip.m_zp)));
      return I*(1.0+p_sk->GF()->KMax(s));
    }

    double Estimate(const Splitting &s) const
    {
      ID_Params ip(s);
      double E(ip.m_I*4.0/(1.0-s.m_z));
      return E*(1.0+p_sk->GF()->KMax(s));
    }

    bool GeneratePoint(Splitting &s) const
    {
      ID_Params ip(s);
      s.m_z=1.0-pow(1.0-ip.m_zp,ran->Get());
      s.m_phi=2.0*M_PI*ran->Get();
      s.m_soft=1;
      return true;
    }

  };// end of class Soft_FS

}// end of namespace ALARIC

using namespace ALARIC;

DECLARE_GETTER(Soft_FS,"FS_Soft",Lorentz,Kernel_Key);

Lorentz *ATOOLS::Getter<Lorentz,Kernel_Key,Soft_FS>::
operator()(const Parameter_Type &args) const
{
  if (args.m_type!=0 || args.m_swap) return NULL;
  if (args.p_v->in[0]==args.p_v->in[1+args.m_mode].Bar() &&
      args.p_v->in[2-args.m_mode].IntSpin()==2) {
    return new Soft_FS(args);
  }
  return NULL;
}

void ATOOLS::Getter<Lorentz,Kernel_Key,Soft_FS>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Soft Lorentz Function";
}
