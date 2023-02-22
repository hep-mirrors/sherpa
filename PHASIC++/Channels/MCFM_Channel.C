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

  class MCFM_Channel : public Single_Channel {

    Vegas *p_vegas;
    int m_type;
    Flavour m_flv;
    double m_alpha, m_ptmin;

  public :

    MCFM_Channel(int _nin,int _nout,Flavour *fl):
      Single_Channel(_nin,_nout,fl)
    {
      name="MCFM_Channel";
      m_incisr=1;
      rannum=3*nout-4+2;
      rans = new double[rannum];
      m_type=1;
      if (nout>1) {
	if (!fl[2].Strong() && !fl[3].Strong()) {
	  if (fl[3]==fl[2].Bar()) m_flv=Flavour(kf_Z);
	  else m_flv=Flavour(kf_Wplus);
	  m_type=2;
	}
      }
      Data_Reader read(" ",";","!","=");
      m_ptmin=read.GetValue<double>("MCFM_JET_PTMIN",30);
      m_alpha=read.GetValue<double>("MCFM_TEXP",2);
      p_vegas = new Vegas(rannum,100,name,0);
    }

    ~MCFM_Channel() { delete p_vegas; }

    void GeneratePoint(Vec4D *p,Cut_Data *cuts,double *_ran)
    {
      DEBUG_FUNC("");
      int irn(0);
      double *ran(p_vegas->GeneratePoint(_ran));
      double S((rpa->gen.PBeam(0)+rpa->gen.PBeam(1)).Abs2());
      double rtS(sqrt(S)), pt2max(S/4);
      weight=2*M_PI/S;
      Vec4D psum(0.,0.,0.,0.);
      for(size_t i(nin+m_type);i<nin+nout;++i) {
        weight/=32.*pow(M_PI,3);
	double pt2min(sqr(cuts->etmin[i])-ms[i]), pt2;
	if (pt2min) {
	  pt2=Channel_Basics::PeakedDist
	    (0.,m_alpha,pt2min,pt2max,1,ran[irn++]);
	  weight*=pow(pt2,m_alpha)*Channel_Basics::
	    PeakedWeight(0.,m_alpha,pt2min,pt2max,pt2,1,rans[irn-1]);
	}
	else {
	  double xr(ran[irn++]);
	  double ptmin(m_ptmin), ptmax(sqrt(pt2max));
	  double pt(2.*ptmin*ptmax*xr/(2.*ptmin+ptmax*(1.-xr)));
	  pt2=sqr(pt);
	  weight*=2.*pt*ptmax/2./ptmin/(2.*ptmin+ptmax)*sqr(2.*ptmin+pt);
	  rans[irn-1]=(pt*2.*ptmin+pt*ptmax)/(2.*ptmin*ptmax+pt*ptmax);
	}
	double yb(log(sqrt(S/4/pt2)+sqrt(S/4/pt2-1.)));
	double ymax(yb), ymin(-yb);
	if (cuts->cosmax[0][i]<1) ymax=std::min(ymax,atanh(cuts->cosmax[0][i]));
	if (cuts->cosmax[1][i]<1) ymin=std::max(ymin,-atanh(cuts->cosmax[1][i]));
        double y(ymin+ran[irn++]*(ymax-ymin));
        weight*=ymax-ymin;
	rans[irn-1]=(y-ymin)/(ymax-ymin);
        double phi(2*M_PI*ran[irn++]);
	weight*=2.*M_PI;
        rans[irn-1]=phi/(2*M_PI);
        double sinhy(sinh(y)), coshy(sqrt(1+sinhy*sinhy));
        double pt(sqrt(pt2)), mt(sqrt(pt2+ms[i]));
        psum+=p[i]=Vec4D(mt*coshy,pt*cos(phi),pt*sin(phi),mt*sinhy);
	DEBUG_VAR(i<<" "<<ran[irn-3]<<" "<<ran[irn-2]<<" "<<ran[irn-1]);
	DEBUG_VAR(i<<" "<<rans[irn-3]<<" "<<rans[irn-2]<<" "<<rans[irn-1]);
      }
      double m2(ms[nin]);
      if (m_type==2) {
	m2=CE.MassivePropMomenta
	  (m_flv.Mass(),m_flv.Width(),1,
	   cuts->scut[2][3],sqr(rtS-psum.Mass()),ran[irn++]);
	weight/=2.*M_PI*CE.MassivePropWeight
	  (m_flv.Mass(),m_flv.Width(),1,
	   cuts->scut[2][3],sqr(rtS-psum.Mass()),m2,rans[irn-1]);
	DEBUG_VAR(m2<<" "<<ran[irn-1]);
	DEBUG_VAR(m2<<" "<<rans[irn-1]);
      }
      double mj2(psum.Abs2()), ptj2(psum.PPerp2());
      double yj(nout>m_type?psum.Y():0);
      double Qt(sqrt(m2+ptj2)), mt(sqrt(mj2+ptj2));
      double ymin(-log(rtS/Qt*(1.-mt/rtS*exp(-yj))));
      double ymax(log(rtS/Qt*(1.-mt/rtS*exp(yj))));
      double yv(ymin+ran[irn++]*(ymax-ymin));
      weight*=ymax-ymin;
      rans[irn-1]=(yv-ymin)/(ymax-ymin);
      DEBUG_VAR(yv<<" "<<ran[irn-1]);
      DEBUG_VAR(yv<<" "<<rans[irn-1]);
      double sinhy(sinh(yv)), coshy(sqrt(1+sinhy*sinhy));
      p[nin]=Vec4D(Qt*coshy,-psum[1],-psum[2],Qt*sinhy);
      double pp((p[nin]+psum).PPlus());
      double pm((p[nin]+psum).PMinus());
      p[0]=Vec4D(pp/2,0,0,pp/2);
      p[1]=Vec4D(pm/2,0,0,-pm/2);
      m_status=pp>0.&&pp<=rpa->gen.PBeam(0).PPlus();
      m_status&=pm>0.&&pm<=rpa->gen.PBeam(1).PMinus();
      if (m_type==1) {
	if (Qt<cuts->etmin[nin]) m_status=0;
      }
      if (m_status && m_type==2) {
	Vec4D pv(p[nin]);
	CE.Isotropic2Momenta
	  (pv,ms[2],ms[3],p[2],p[3],
	   ran[irn++],ran[irn++],-1,1);
	DEBUG_VAR("v "<<ran[irn-2]<<" "<<ran[irn-1]);
	weight/=sqr(2*M_PI)*CE.Isotropic2Weight
	  (p[2],p[3],rans[irn-1],rans[irn-2],-1,1);
	DEBUG_VAR("v "<<rans[irn-2]<<" "<<rans[irn-1]);
      }
      weight*=p_vegas->GenerateWeight(rans);
      sprime=(p[0]+p[1]).Abs2();
      y=(p[0]+p[1]).Y();
      DEBUG_VAR(weight);
    }

    void GenerateWeight(Vec4D *p,Cut_Data *cuts,bool recompute)
    {
      if (!recompute) return;
      DEBUG_FUNC("");
      int irn(0);
      double S((rpa->gen.PBeam(0)+rpa->gen.PBeam(1)).Abs2());
      double rtS(sqrt(S)), pt2max(S/4);
      weight=2*M_PI/S;
      Vec4D psum(0.,0.,0.,0.);
      for(size_t i(nin+m_type);i<nin+nout;++i) {
        weight/=32.*pow(M_PI,3);
	double pt2min(sqr(cuts->etmin[i])-ms[i]);
        double pt2(p[i].PPerp2());
	if (pt2min) {
	  weight*=pow(pt2,m_alpha)*Channel_Basics::
	    PeakedWeight(0.,m_alpha,pt2min,pt2max,pt2,1,rans[irn++]);
	}
	else {
	  double ptmin(m_ptmin), ptmax(sqrt(pt2max)), pt(sqrt(pt2));
	  weight*=2.*pt*ptmax/2./ptmin/(2.*ptmin+ptmax)*sqr(2.*ptmin+pt);
	  rans[irn++]=(pt*2.*ptmin+pt*ptmax)/(2.*ptmin*ptmax+pt*ptmax);
	}
	double yb(log(sqrt(S/4/pt2)+sqrt(S/4/pt2-1.)));
	double ymax(yb), ymin(-yb);
	if (cuts->cosmax[0][i]<1) ymax=std::min(ymax,atanh(cuts->cosmax[0][i]));
	if (cuts->cosmax[1][i]<1) ymin=std::max(ymin,-atanh(cuts->cosmax[1][i]));
        double y(p[i].Y());
        weight*=ymax-ymin;
	rans[irn++]=(y-ymin)/(ymax-ymin);
        double phi(p[i].Phi());
	if (phi<0) phi+=2*M_PI;
	weight*=2.*M_PI;
        rans[irn++]=phi/(2*M_PI);
        psum+=p[i];
	DEBUG_VAR(i<<" "<<rans[irn-3]<<" "<<rans[irn-2]<<" "<<rans[irn-1]);
      }
      double m2(ms[nin]), yv(p[nin].Y());
      if (m_type==2) {
	m2=(p[2]+p[3]).Abs2();
	yv=(p[2]+p[3]).Y();
	weight/=2.*M_PI*CE.MassivePropWeight
	  (m_flv.Mass(),m_flv.Width(),1,
	   cuts->scut[2][3],sqr(rtS-psum.Mass()),m2,rans[irn++]);
	DEBUG_VAR(m2<<" "<<rans[irn-1]);
      }
      double mj2(psum.Abs2()), ptj2(psum.PPerp2());
      double yj(nout>m_type?psum.Y():0);
      double Qt(sqrt(m2+ptj2)), mt(sqrt(mj2+ptj2));
      double ymin(-log(rtS/Qt*(1.-mt/rtS*exp(-yj))));
      double ymax(log(rtS/Qt*(1.-mt/rtS*exp(yj))));
      weight*=ymax-ymin;
      rans[irn++]=(yv-ymin)/(ymax-ymin);
      DEBUG_VAR(yv<<" "<<rans[irn-1]);
      if (m_type==2) {
	weight/=sqr(2*M_PI)*CE.Isotropic2Weight
	  (p[2],p[3],rans[irn++],rans[irn++],-1,1);
	DEBUG_VAR("v "<<rans[irn-2]<<" "<<rans[irn-1]);
      }
      weight*=p_vegas->GenerateWeight(rans);
      DEBUG_VAR(weight);
    }
    
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId);   }

    void AddPoint(double w) { p_vegas->AddPoint(w,rans); }

    void MPISync() { p_vegas->MPISync(); }
    void Optimize() { p_vegas->Optimize(); }
    void EndOptimize() { p_vegas->EndOptimize(); }

    void ISRInfo(int &type,double &mass,double &width)
    { type=-2; mass=0.; width=0.; }

    std::string ChID() { return "MCFM_Channel"; }

  };// end of class MCFM_Channel

  class MCFM_Channel_Generator: public Channel_Generator {
  public:
    MCFM_Channel_Generator(const Channel_Generator_Key &key):
    Channel_Generator(key) {}
    int GenerateChannels()
    {
      p_proc->FillIntegrator(&*p_proc->Integrator()->PSHandler());
      p_mc->DropAllChannels();
      p_mc->Add(new MCFM_Channel
		(p_proc->NIn(),p_mc->Nout(),
		 (Flavour*)&p_proc->Flavours().front()));
      return 0;
    }
  };// end of class MCFM_Channel_Generator

}// end of namespace PHASIC

using namespace PHASIC;

DECLARE_GETTER(MCFM_Channel_Generator,"MCFM",
	       Channel_Generator,Channel_Generator_Key);
Channel_Generator *Getter
<Channel_Generator,Channel_Generator_Key,MCFM_Channel_Generator>::
operator()(const Channel_Generator_Key &args) const
{
  return new MCFM_Channel_Generator(args);
}
void Getter<Channel_Generator,Channel_Generator_Key,
		    MCFM_Channel_Generator>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"MCFM integrator";
}
