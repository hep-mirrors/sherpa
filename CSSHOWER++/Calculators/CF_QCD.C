#include "CSSHOWER++/Showers/Splitting_Function_Base.H"

#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

#include <algorithm>

namespace CSSHOWER {
  
  const double s_Nc = 3.;
  const double s_CF = (s_Nc*s_Nc-1.)/(2.*s_Nc);
  const double s_CA = s_Nc;
  const double s_TR = 1./2.;

  class CF_QCD: public SF_Coupling {
  protected:

    //! Main underlying coupling (set by SetCoupling)
    MODEL::Running_AlphaS *p_cpl;

    //! (Temporary) alternative coupling (set by SetAlternativeUnderlyingCoupling)
    MODEL::One_Running_AlphaS *p_altcpl;

    //! Alternative renormalisation scale factor (set by SetAlternativeUnderlyingCoupling)
    double m_altrsf;

    //! Buffer of max alphas values to avoid re-calculations
    std::map<MODEL::One_Running_AlphaS *, double> m_altcplmax;

    double m_q, m_rsf, m_k0sq;
    int m_scvmode;

  public:

    inline CF_QCD(const SF_Key &key):
      SF_Coupling(key), p_altcpl(NULL), m_altcplmax(), m_q(0.), m_rsf(1.), m_k0sq(0.0), m_scvmode(0)
    {
      if (key.p_v->in[0].StrongCharge()==8 &&
	  key.p_v->in[1].StrongCharge()==8 &&
	  key.p_v->in[2].StrongCharge()==8) m_q=s_CA;
      else m_q=(key.p_v->in[0].StrongCharge()==8)?s_TR:s_CF;
      if (key.m_type==cstp::FF || key.m_type==cstp::FI) {
	if (key.p_v->in[0].StrongCharge()==8) m_q/=2.0;
      }
      else {
	if (key.m_mode==0) {
	  if (key.p_v->in[1].StrongCharge()==8) m_q/=2.0;
	}
	else {
	  if (key.p_v->in[2].StrongCharge()==8) m_q/=2.0;
	}
      }
    }

    bool SetCoupling(MODEL::Model_Base *md,
		     const double &k0sqi,const double &k0sqf,
		     const double &isfac,const double &fsfac);
    template<class T>
    double CplMax(T * as, double rsf) const;
    double Coupling(const double &scale,const int pol);
    bool AllowSpec(const ATOOLS::Flavour &fl);
    double CplFac(const double &scale) const;

    bool AllowsAlternativeCouplingUsage() const { return true; }
    void SetAlternativeUnderlyingCoupling(void *, double sf=1.0);
  };

}

using namespace CSSHOWER;
using namespace MODEL;
using namespace ATOOLS;

bool CF_QCD::SetCoupling(MODEL::Model_Base *md,
			 const double &k0sqi,const double &k0sqf,
			 const double &isfac,const double &fsfac)
{
  p_cpl=(MODEL::Running_AlphaS*)md->GetScalarFunction("alpha_S");
  SetAlternativeUnderlyingCoupling(NULL);
  m_altcplmax.clear(); // buffered values are not valid anymore
  m_rsf=ToType<double>(rpa->gen.Variable("RENORMALIZATION_SCALE_FACTOR"))
        *ToType<double>(rpa->gen.Variable("CSS_SCALE_FACTOR"));
  m_scvmode=ToType<int>(rpa->gen.Variable("CSS_SCALE_VARIATION_SCHEME"));
  m_cplfac=((m_type/10==1)?fsfac:isfac);
  m_k0sq=(m_type/10==1)?k0sqf:k0sqi;
  m_cplmax.push_back(CplMax(p_cpl, m_rsf));
  m_cplmax.push_back(0.0);
  return true;
}

void CF_QCD::SetAlternativeUnderlyingCoupling(void *cpl, double sf)
{
  m_altrsf = sf;
  if (cpl == NULL) {
    p_altcpl = NULL;
    return;
  } else {
    p_altcpl = static_cast<MODEL::One_Running_AlphaS *>(cpl);
    if (m_altcplmax.find(p_altcpl) == m_altcplmax.end()) {
      m_altcplmax[p_altcpl] = CplMax(p_altcpl, m_altrsf);
    }
  }
}

template<class T>
double CF_QCD::CplMax(T * as, double rsf) const
{
  double minscale = CplFac(m_k0sq)*m_k0sq;
  double ct(0.);
  const double boundedscale(Max(as->ShowerCutQ2(), minscale));
  if (m_rsf>1.) // only for f>1 cpl gets larger
    ct=-(*as)[boundedscale]/M_PI*as->Beta0(0.)*log(m_rsf);
  return (*as)[boundedscale]*(1.-ct)*m_q;
}

double CF_QCD::Coupling(const double &scale,const int pol)
{
  if (pol!=0) return 0.0; // we do not update m_last when polarized
  
  // use nominal or alternative coupling and scale factor
  One_Running_AlphaS * const as = (p_altcpl) ? p_altcpl : p_cpl->GetAs();
  const double rsf = (p_altcpl) ? m_altrsf : m_rsf;

  if (scale<0.0) return m_last = (*as)(sqr(rpa->gen.Ecms()))*m_q;
  double t(CplFac(scale)*scale), scl(CplFac(scale)*scale*rsf);
  if (t < as->ShowerCutQ2()) return m_last = 0.0;
  double cpl=(*as)[t];
  if (!IsEqual(scl,t)) {
    msg_Debugging()<<"as(\\mu_R^2)="<<cpl<<std::endl;
    std::vector<double> ths(as->Thresholds(t,scl));
    ths.push_back((scl>t)?scl:t);
    ths.insert(ths.begin(),(scl>t)?t:scl);
    if (t<scl) std::reverse(ths.begin(),ths.end());
    msg_Debugging()<<"thresholds: "<<ths<<std::endl;
    double fac(1.),ct(0.);
    // Beta0 from One_Running_AlphaS contains extra factor 1/2
    switch (m_scvmode) {
    case 1:
      // replace as(t) -> as(t)*prod[1-as/2pi*beta(nf)*log(th[i]/th[i-1])]
      for (size_t i(1);i<ths.size();++i) {
        ct=cpl/M_PI*as->Beta0((ths[i]+ths[i-1])/2.0)*log(ths[i]/ths[i-1]);
        fac*=1.0-ct;
      }
      break;
    case 2:
      // replace as(t) -> as(t)*[1-sum as/2pi*beta(nf)*log(th[i]/th[i-1])]
      for (size_t i(1);i<ths.size();++i)
        ct+=cpl/M_PI*as->Beta0((ths[i]+ths[i-1])/2.0)*log(ths[i]/ths[i-1]);
      fac=1.-ct;
      break;
    default:
      fac=1.;
      break;
    }
    msg_Debugging()<<"ct="<<ct<<std::endl;
    if (fac<0.) {
      msg_Tracking()<<METHOD<<"(): Renormalisation term too large. Remove."
                    <<std::endl;
      fac=1.;
    }
    cpl*=fac;
    msg_Debugging()<<"as(\\mu_R^2)*(1-ct)="<<cpl<<std::endl;
  }
  cpl*=m_q;
  const double cplmax = (p_altcpl) ? m_altcplmax[p_altcpl] : m_cplmax.front();
  if (cpl>cplmax) {
    msg_Error()<<METHOD<<"(): Value exceeds maximum at t = "
               <<sqrt(t)<<" -> \\mu_R = "<<sqrt(scl)
               <<", qmin = "<<sqrt(as->CutQ2())<<std::endl;
    return m_last = cplmax;
  }
#ifdef DEBUG__Trial_Weight
  msg_Debugging()<<"as weight kt = "<<sqrt(CplFac(scale))<<" * "
		 <<sqrt(scale)<<", \\alpha_s("<<sqrt(scl)<<") = "
		 <<(*as)[scl]<<", m_q = "<<m_q<<"\n";
#endif
  return m_last = cpl;
}

bool CF_QCD::AllowSpec(const ATOOLS::Flavour &fl) 
{
  if (m_type==cstp::FF && !fl.Strong()) return true;
  if (abs(fl.StrongCharge())==3) {
    switch (m_type) {
    case cstp::FF: 
      if (abs(p_lf->FlA().StrongCharge())==3)
	return p_lf->FlA().StrongCharge()==-fl.StrongCharge();
      break;
    case cstp::FI: 
      if (abs(p_lf->FlA().StrongCharge())==3)
	return p_lf->FlA().StrongCharge()==fl.StrongCharge();
      break;
    case cstp::IF: 
      if (abs(p_lf->FlB().StrongCharge())==3)
	return p_lf->FlB().StrongCharge()==fl.StrongCharge();
      break;
    case cstp::II: 
      if (abs(p_lf->FlB().StrongCharge())==3)
	return p_lf->FlB().StrongCharge()==-fl.StrongCharge();
      break;
    case cstp::none: abort();
    }
  }
  return fl.Strong();
}

double CF_QCD::CplFac(const double &scale) const
{
  if (m_kfmode==0) return m_cplfac;
  One_Running_AlphaS * const as = (p_altcpl) ? p_altcpl : p_cpl->GetAs();
  double nf=as->Nf(scale);
  double kfac=exp(-(67.0-3.0*sqr(M_PI)-10.0/3.0*nf)/(33.0-2.0*nf));
  return m_cplfac*kfac;
}

namespace CSSHOWER {

DECLARE_CPL_GETTER(CF_QCD_Getter);

SF_Coupling *CF_QCD_Getter::operator()
  (const Parameter_Type &args) const
{
  return new CF_QCD(args);
}

void CF_QCD_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"strong coupling";
}

}

DECLARE_GETTER(CF_QCD_Getter,"SF_QCD_Fill",
	       void,SFC_Filler_Key);

void *ATOOLS::Getter<void,SFC_Filler_Key,CF_QCD_Getter>::
operator()(const SFC_Filler_Key &key) const
{
  DEBUG_FUNC("model = "<<key.p_md->Name());
  const Vertex_Table *vtab(key.p_md->VertexTable());
  for (Vertex_Table::const_iterator
	 vlit=vtab->begin();vlit!=vtab->end();++vlit) {
    for (Vertex_List::const_iterator 
	   vit=vlit->second.begin();vit!=vlit->second.end();++vit) {
      Single_Vertex *v(*vit);
      if (v->NLegs()>3) continue;
      if (v->Color.front().Type()!=cf::T &&
	  v->Color.front().Type()!=cf::F) continue;
      msg_Debugging()<<"Add "<<v->in[0].Bar()<<" -> "<<v->in[1]<<" "<<v->in[2]<<" {\n";
      std::string atag("{"+v->in[0].Bar().IDName()+"}");
      std::string btag("{"+v->in[1].IDName()+"}");
      std::string ctag("{"+v->in[2].IDName()+"}");
      key.p_gets->push_back(new CF_QCD_Getter(atag+btag+ctag));
    }
  }
  return NULL;
}

void ATOOLS::Getter<void,SFC_Filler_Key,CF_QCD_Getter>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"qcd coupling filler";
}

