#include "PHASIC++/Channels/Single_Channel.H"
#include "PHASIC++/Channels/Channel_Generator.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Channels/Vegas.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;

namespace PHASIC {

  class Chili_Channel : public Single_Channel {

    Vegas *p_vegas;
    int m_type;
    double m_alpha, m_ptmin, m_flmass, m_flwidth;

  public :

    Chili_Channel(int _nin,int _nout,Flavour *fl):
      Single_Channel(_nin,_nout,fl)
    {
      m_name="Chili Channel";
      m_incisr=1;
      m_rannum=3*m_nout-4+2;
      p_rans = new double[m_rannum];
      m_type=1;
      if (m_nout>1) {
	if (!fl[2].Strong() && !fl[3].Strong()) {
	  if (fl[3]==fl[2].Bar()) {
	    m_flmass=Flavour(kf_Z).Mass();
	    m_flwidth=Flavour(kf_Z).Width();
	  }
	  else {
	    m_flmass=Flavour(kf_Wplus).Mass();
	    m_flwidth=Flavour(kf_Wplus).Width();
	  }
	  m_type=2;
	}
      }
      Scoped_Settings s{ Settings::GetMainSettings()["PSI"] };
      m_ptmin=s["Chili_JET_PTMIN"].SetDefault(30).Get<double>();
      m_alpha=s["Chili_TEXP"].SetDefault(2).Get<double>();
      m_flmass=s["Chili_RESONANCE_MASS"].ResetDefault().
	SetDefault(m_flmass).Get<double>();
      m_flwidth=s["Chili_RESONANCE_WIDTH"].ResetDefault().
	SetDefault(m_flwidth).Get<double>();
      p_vegas = new Vegas(m_rannum,100,m_name);
    }

    ~Chili_Channel() { delete p_vegas; }

    void GeneratePoint(Vec4D *p,Cut_Data *cuts,double *_ran)
    {
      DEBUG_FUNC("");
      int irn(0);
      double *ran(p_vegas->GeneratePoint(_ran));
      double S((rpa->gen.PBeam(0)+rpa->gen.PBeam(1)).Abs2());
      double rtS(sqrt(S)), pt2max(S/4), nnoj(m_type);
      m_weight=2*M_PI/S;
      Vec4D psum(0.,0.,0.,0.);
      for(size_t j(0);j<2;++j) {
	Vec4D psum2(psum.Perp()/nnoj);
      for(size_t i(m_nin+m_type);i<m_nin+m_nout;++i) {
	double pt2min(sqr(cuts->etmin[i])-p_ms[i]), pt2;
	if (pt2min) {
	  if (j==1) continue;
	  pt2=PeakedDist
	    (0.,m_alpha,pt2min,pt2max,1,ran[irn++]);
	  m_weight*=pow(pt2,m_alpha)*
	    PeakedWeight(0.,m_alpha,pt2min,pt2max,pt2,1,p_rans[irn-1]);
	}
	else {
	  if (j==0) {
	    ++nnoj;
	    continue;
	  }
	  double xr(ran[irn++]);
	  double ptmin(std::max(m_ptmin,sqrt(p_ms[i]))), ptmax(sqrt(pt2max));
	  double pt(2.*ptmin*ptmax*xr/(2.*ptmin+ptmax*(1.-xr)));
	  pt2=sqr(pt);
	  m_weight*=2.*pt*ptmax/2./ptmin/(2.*ptmin+ptmax)*sqr(2.*ptmin+pt);
	  p_rans[irn-1]=(pt*2.*ptmin+pt*ptmax)/(2.*ptmin*ptmax+pt*ptmax);
	}
        m_weight/=32.*pow(M_PI,3);
	double yb(log(sqrt(S/4/pt2)+sqrt(S/4/pt2-1.)));
	double ymax(yb), ymin(-yb);
        double y(ymin+ran[irn++]*(ymax-ymin));
        m_weight*=ymax-ymin;
	p_rans[irn-1]=(y-ymin)/(ymax-ymin);
        double phi(2*M_PI*ran[irn++]);
	m_weight*=2.*M_PI;
        p_rans[irn-1]=phi/(2*M_PI);
        double sinhy(sinh(y)), coshy(sqrt(1+sinhy*sinhy));
        double pt(sqrt(pt2));
	pt2=(Vec4D(0.,pt*cos(phi),pt*sin(phi),0.)-psum2).PPerp2();
	double mt(sqrt(pt2+p_ms[i]));
	psum+=p[i]=Vec4D(mt*coshy,pt*cos(phi),pt*sin(phi),mt*sinhy)-psum2;
	DEBUG_VAR(i<<" "<<ran[irn-3]<<" "<<ran[irn-2]<<" "<<ran[irn-1]);
	DEBUG_VAR(i<<" "<<p_rans[irn-3]<<" "<<p_rans[irn-2]<<" "<<p_rans[irn-1]);
      }
      }
      double m2(p_ms[m_nin]);
      if (m_type==2) {
	m2=CE.MassivePropMomenta
	  (m_flmass,m_flwidth,
	   cuts->scut[2][3],sqr(rtS-psum.Mass()),ran[irn++]);
	m_weight/=2.*M_PI*CE.MassivePropWeight
	  (m_flmass,m_flwidth,
	   cuts->scut[2][3],sqr(rtS-psum.Mass()),m2,p_rans[irn-1]);
	DEBUG_VAR(m2<<" "<<ran[irn-1]);
	DEBUG_VAR(m2<<" "<<p_rans[irn-1]);
      }
      double mj2(psum.Abs2()), ptj2(psum.PPerp2());
      double yj(m_nout>m_type?psum.Y():0);
      double Qt(sqrt(m2+ptj2)), mt(sqrt(mj2+ptj2));
      double ymin(-log(rtS/Qt*(1.-mt/rtS*exp(-yj))));
      double ymax(log(rtS/Qt*(1.-mt/rtS*exp(yj))));
      double yv(ymin+ran[irn++]*(ymax-ymin));
      m_weight*=ymax-ymin;
      p_rans[irn-1]=(yv-ymin)/(ymax-ymin);
      DEBUG_VAR(yv<<" "<<ran[irn-1]);
      DEBUG_VAR(yv<<" "<<p_rans[irn-1]);
      double sinhy(sinh(yv)), coshy(sqrt(1+sinhy*sinhy));
      p[m_nin]=Vec4D(Qt*coshy,-psum[1],-psum[2],Qt*sinhy);
      double pp((p[m_nin]+psum).PPlus());
      double pm((p[m_nin]+psum).PMinus());
      p[0]=Vec4D(pp/2,0,0,pp/2);
      p[1]=Vec4D(pm/2,0,0,-pm/2);
      m_status=pp>0.&&pp<=rpa->gen.PBeam(0).PPlus();
      m_status&=pm>0.&&pm<=rpa->gen.PBeam(1).PMinus();
      if (m_type==1) {
	if (Qt<cuts->etmin[m_nin]) m_status=0;
      }
      if (m_status && m_type==2) {
	Vec4D pv(p[m_nin]);
	CE.Isotropic2Momenta
	  (pv,p_ms[2],p_ms[3],p[2],p[3],
	   ran[irn++],ran[irn++],-1,1);
	DEBUG_VAR("v "<<ran[irn-2]<<" "<<ran[irn-1]);
	m_weight/=sqr(2*M_PI)*CE.Isotropic2Weight
	  (p[2],p[3],p_rans[irn-1],p_rans[irn-2],-1,1);
	DEBUG_VAR("v "<<p_rans[irn-2]<<" "<<p_rans[irn-1]);
      }
      m_weight*=p_vegas->GenerateWeight(p_rans);
      m_sprime=(p[0]+p[1]).Abs2();
      m_y=(p[0]+p[1]).Y();
      DEBUG_VAR(m_weight);
    }

    void GenerateWeight(Vec4D *p,Cut_Data *cuts,bool recompute)
    {
      if (!recompute) return;
      DEBUG_FUNC("");
      int irn(0);
      double S((rpa->gen.PBeam(0)+rpa->gen.PBeam(1)).Abs2());
      double rtS(sqrt(S)), pt2max(S/4), nnoj(m_type);
      m_weight=2*M_PI/S;
      Vec4D psum(0.,0.,0.,0.);
      for(size_t j(0);j<2;++j) {
	Vec4D psum2(psum.Perp()/nnoj);
      for(size_t i(m_nin+m_type);i<m_nin+m_nout;++i) {
	double pt2min(sqr(cuts->etmin[i])-p_ms[i]);
        double pt2((p[i]+psum2).PPerp2());
	if (pt2min) {
	  if (j==1) continue;
	  m_weight*=pow(pt2,m_alpha)*
	    PeakedWeight(0.,m_alpha,pt2min,pt2max,pt2,1,p_rans[irn++]);
	}
	else {
	  if (j==0) {
	    ++nnoj;
	    continue;
	  }
	  double ptmin(std::max(m_ptmin,sqrt(p_ms[i])));
	  double ptmax(sqrt(pt2max)), pt(sqrt(pt2));
	  m_weight*=2.*pt*ptmax/2./ptmin/(2.*ptmin+ptmax)*sqr(2.*ptmin+pt);
	  p_rans[irn++]=(pt*2.*ptmin+pt*ptmax)/(2.*ptmin*ptmax+pt*ptmax);
	}
        m_weight/=32.*pow(M_PI,3);
	double yb(log(sqrt(S/4/pt2)+sqrt(S/4/pt2-1.)));
	double ymax(yb), ymin(-yb);
        double y(p[i].Y());
        m_weight*=ymax-ymin;
	p_rans[irn++]=(y-ymin)/(ymax-ymin);
        double phi((p[i]+psum2).Phi());
	if (phi<0) phi+=2*M_PI;
	m_weight*=2.*M_PI;
        p_rans[irn++]=phi/(2*M_PI);
        psum+=p[i];
	DEBUG_VAR(i<<" "<<p_rans[irn-3]<<" "<<p_rans[irn-2]<<" "<<p_rans[irn-1]);
      }
      }
      double m2(p_ms[m_nin]), yv(p[m_nin].Y());
      if (m_type==2) {
	m2=(p[2]+p[3]).Abs2();
	yv=(p[2]+p[3]).Y();
	m_weight/=2.*M_PI*CE.MassivePropWeight
	  (m_flmass,m_flwidth,
	   cuts->scut[2][3],sqr(rtS-psum.Mass()),m2,p_rans[irn++]);
	DEBUG_VAR(m2<<" "<<p_rans[irn-1]);
      }
      double mj2(psum.Abs2()), ptj2(psum.PPerp2());
      double yj(m_nout>m_type?psum.Y():0);
      double Qt(sqrt(m2+ptj2)), mt(sqrt(mj2+ptj2));
      double ymin(-log(rtS/Qt*(1.-mt/rtS*exp(-yj))));
      double ymax(log(rtS/Qt*(1.-mt/rtS*exp(yj))));
      m_weight*=ymax-ymin;
      p_rans[irn++]=(yv-ymin)/(ymax-ymin);
      DEBUG_VAR(yv<<" "<<p_rans[irn-1]);
      if (m_type==2) {
	m_weight/=sqr(2*M_PI)*CE.Isotropic2Weight
	  (p[2],p[3],p_rans[irn++],p_rans[irn++],-1,1);
	DEBUG_VAR("v "<<p_rans[irn-2]<<" "<<p_rans[irn-1]);
      }
      m_weight*=p_vegas->GenerateWeight(p_rans);
      DEBUG_VAR(m_weight);
    }
    
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId);   }

    void AddPoint(double w) { p_vegas->AddPoint(w,p_rans); }

    void MPISync() { p_vegas->MPISync(); }
    void Optimize() { p_vegas->Optimize(); }
    void EndOptimize() { p_vegas->EndOptimize(); }

    void ISRInfo(int &type,double &mass,double &width)
    { type=-2; mass=0.; width=0.; }

    std::string ChID() { return "Chili_Channel"; }

  };// end of class Chili_Channel

  class Chili_Channel_Generator: public Channel_Generator {
  public:
    Chili_Channel_Generator(const Channel_Generator_Key &key):
    Channel_Generator(key) {}
    int GenerateChannels()
    {
      p_proc->FillIntegrator(&*p_proc->Integrator()->PSHandler());
      p_mc->DropAllChannels();
      p_mc->Add(new Chili_Channel
		(p_proc->NIn(),p_proc->NOut(),
		 (Flavour*)&p_proc->Flavours().front()));
      return 0;
    }
  };// end of class Chili_Channel_Generator

}// end of namespace PHASIC

using namespace PHASIC;

DECLARE_GETTER(Chili_Channel_Generator,"Chili",
	       Channel_Generator,Channel_Generator_Key);
Channel_Generator *Getter
<Channel_Generator,Channel_Generator_Key,Chili_Channel_Generator>::
operator()(const Channel_Generator_Key &args) const
{
  return new Chili_Channel_Generator(args);
}
void Getter<Channel_Generator,Channel_Generator_Key,
		    Chili_Channel_Generator>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Chili integrator";
}
